#
# Reverse Engineering 3D Printer Files to CAD models.
#
# RepRap Ltd
# https://reprapltd.com
#
# Written by: Adrian Bowyer
#
# 3 September 2021
#
# Licence: GPL
#
# ply2set.py is a program to convert PLY triangulation files to set-theoretic or CSG solid
# models. There is a companion program - set2OpenSCAD.py that converts the results
# to openSCAD models, which can be used in OpenSCAD (clearly) and also FreeCAD.
#
# PLY files are easily created from STL files, which are universally used as input
# to 3D printers. For example the Java program Triangulation.java in this
# repository will do that conversion. Thus this program allows the reverse engineering
# of a 3D Printer file back to a CAD file.
#
# The program uses a modification of Tony Woo's alternating sum of volumes algorithm, which is described
# in his two papers in the Documentation directory of this repository. At the moment it is
# a partial implementation; there are some pathological cases with which it won't deal.
#
# Here are some useful and explanatory web pages:
#
# Ply files: https://en.wikipedia.org/wiki/PLY_(file_format)
# Ply files in Python: https://github.com/dranjan/python-plyfile
# Convex hulls: https://en.wikipedia.org/wiki/Convex_hull
# Convex hulls in Python: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
# Set-theoretic or CSG solid models: https://en.wikipedia.org/wiki/Constructive_solid_geometry
# Boolean logic in Python: https://booleanpy.readthedocs.io/en/latest/index.html
# OpenSCAD: https://openscad.org/
# FreeCAD: https://www.freecadweb.org/
# STL to PLY conversion free online: https://products.aspose.app/3d/conversion/stl-to-ply
# Infix to postfix/reverse polish conversion: http://csis.pace.edu/~wolf/CS122/infix-postfix.htm
#

import copy
import sys
import random
import re
import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
import boolean
from plyfile import PlyData
from scipy.spatial import ConvexHull

# Fudge factor
small = 0.000001

# What (if anything) to plot
#
# 0 - no graphics
# 1 - input triangles at each stage
# 2 - input triangles at each stage plus the CH
# 3 - input and output triangles at each stage plus the convex hulls
graphics = 0

# True to print debugging information
debug = False

# Set operators
leaf = 0
union = 1
intersection = 2
subtraction = 3
empty = 4
universal = 5

# The bounding box and a point near the middle
positiveCorner = [-sys.float_info.max, -sys.float_info.max, -sys.float_info.max]
negativeCorner = [sys.float_info.max, sys.float_info.max, sys.float_info.max]
middle = [0, 0, 0]

# We want there to be just one global instance of this
halfSpaces = None

# How deep to run the recursion before bailing out because
# we may have a pathological infinite-recursion case
maximumLevel = 10

# The Set class uses this to represent the set theoretic expression
booleanAlgebra = boolean.BooleanAlgebra()
TRUE, FALSE, NOT, AND, OR, symbol = booleanAlgebra.definition()

#*******************************************************************************************************

# Small holding class to load the PLY file and provide simple access to its data.

class TriangleFileData:
 def __init__(self, fileName):
  self.plyFileData = PlyData.read(fileName)
  # Ply files can hold polygons other than triangles. Check...
  for face in self.plyFileData['face']:
   if len(face[0]) != 3:
    error = "\nA polygon in ply file " + fileName + " has " + str(len(face[0])) + " vertices (i.e. it's not a triangle).\n"
    error += "ply2set.py can only work with triangles.\n"
    sys.exit(error)

 def VertexCount(self):
  return len(self.plyFileData['vertex'])

 def Vertex(self, index):
  return list(self.plyFileData['vertex'][index])

 def TriangleCount(self):
  return len(self.plyFileData['face'])

 def Triangle(self, index):
  return self.plyFileData['face'][index][0]


#************************************************************************************************************

# Class to hold a single triangle. The three vertices are indices into the list of vertices
# that have been previously loaded from a PLY file into ply.

class Triangle:
 def __init__(self, vertices, colour, ply):
  self.vertices = vertices
  self.ply = ply
  self.points = self.Points()
  self.normal = np.cross( np.subtract(self.points[1], self.points[0]), np.subtract(self.points[2], self.points[0]) )
  s2 = maths.sqrt(np.dot(self.normal, self.normal))
  if s2 < small:
   print("Triangle: vertices are collinear!")
  self.normal = np.multiply(self.normal, 1.0/s2)
  self.colour = RandomShade(colour)
  self.centroid = np.multiply(np.add(np.add(self.points[0], self.points[1]), self.points[2]), 1.0/3.0)
  # This will be an index into halfSpaces, not the halfspace itself
  self.halfSpace = -1

 # The three actual space (x, y, z) coordinates of the vertices
 def Points(self):
  return [self.ply.Vertex(self.vertices[0]), self.ply.Vertex(self.vertices[1]), self.ply.Vertex(self.vertices[2])]

 # The list of vertex indices.
 def Vertices(self):
  return self.vertices

 def SetColour(self, colour):
  self.colour = RandomShade(colour)

 # Change the sense of the triangle. The normal should point from solid to air.
 def Invert(self):
  temp = self.vertices[0]
  self.vertices[0] = self.vertices[1]
  self.vertices[1] = temp
  self.halfSpace = halfSpaces.Get(self.halfSpace).complement
  self.normal = [-self.normal[0], -self.normal[1], -self.normal[2]]
  self.points = self.Points()

 # Is a point inside the triangle?
 # Note - this is doing an essentially 2D calculation in 3D space.
 # To see how this works write it out as vector algebra.
 def ContainsPoint(self, point):
  for p0 in range(3):
   p1 = (p0+1)%3
   p2 = (p0+2)%3
   p10 = np.subtract(self.points[p1], self.points[p0])
   p20 = np.subtract(self.points[p2], self.points[p0])
   d12 = np.cross(p10, p20)
   pp0 = np.subtract(point, self.points[p0])
   dp0 = np.cross(p10, pp0)
   if np.dot(d12, dp0) < 0:
    return False
  return True

 def __str__(self):
  corners = self.Points()
  result = "<\n"
  for c in corners:
   result += " " + str(c) + "\n"
  result += ">\n"
  return result

# Take a base RGB colour and perturb it a bit. Used to allow coplanar triangles
# to be distinguished.
def RandomShade(baseColour):
 colour = [random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1)]
 for c in range(3):
  colour[c] += baseColour[c]
  if colour[c] < 0:
   colour[c] = 0
  if colour[c] > 1:
   colour[c] = 1
 return colour

#*************************************************************************************************************

# A planar half space of the form Ax + By + Cz + D <= 0, where (A, B, C) is the
# plane's normal vector and -D is the distance from it to the origin.
#
# It is the plane through a triangle.

class HalfSpace:
 def __init__(self, triangle = None):
  self.complement = -1
  self.flag = False
  if triangle is None:
   return
  self.normal = triangle.normal
  self.d = -np.dot(self.normal, triangle.centroid)

# Note: two coincident half spaces of opposite sense are not considered equal.
 def __eq__(self, halfSpace):
  if 1.0 - np.dot(halfSpace.normal, self.normal) > small:
   return False
  if abs(self.d - halfSpace.d) > small:
   return False
  return True

# Is a halfspace my complement?
 def Opposite(self, halfSpace):
  if 1.0 + np.dot(halfSpace.normal, self.normal) > small:
   return False
  if abs(self.d + halfSpace.d) > small:
   return False
  return True

 # Work out the complement
 def Complement(self):
  result = HalfSpace()
  result.normal = [-self.normal[0], -self.normal[1], -self.normal[2]]
  result.d = -self.d
  result.complement = -1
  return result

 # Used for output so that only halfSpaces that are actually referred to
 # in a set-theoretic expression are saved.
 def SetFlag(self, f):
  self.flag = f

 # Get the coeficients as an array.
 def Coeficients(self):
  return np.array([self.normal[0], self.normal[1], self.normal[2], self.d])

 def __str__(self):
  result = "{" + str(self.normal[0]) + "x + " + str(self.normal[1]) + "y + " +\
   str(self.normal[2]) + "z + " + str(self.d) + " <= 0}"
  result = result.replace("+ -", "- ") #Sigh!
  return result


#********************************************************************************************************

# A list of half spaces

class HalfSpaceList:
 def __init__(self):
  self.halfSpaceList = []

 # Is halfSpace in the list? If so return an index to it. If not, return -1
 def LookUp(self, halfSpace):
  for hs in range(len(self.halfSpaceList)):
   listHalfSpace = self.halfSpaceList[hs]
   if halfSpace == listHalfSpace:
    return hs
  return -1

 # Add the halfspace to the list, unless that half space is
 # already in the list. Return the index of the half space in the list.
 # Also add its complement, so we can flip between them.
 def AddSpace(self, halfSpace):
  last = len(self.halfSpaceList)
  inList = self.LookUp(halfSpace)
  if inList < 0:
   self.halfSpaceList.append(halfSpace)
   self.halfSpaceList.append(halfSpace.Complement())
   halfSpace.complement = last + 1
   self.halfSpaceList[last+1].complement = last
   return last
  else:
   return inList

 # Add the halfspace corresponding to a triangle using the function above.
 def AddTriangle(self, triangle):
  halfSpace = HalfSpace(triangle)
  index = self.AddSpace(halfSpace)
  triangle.halfSpace = index
  return index

 # Set all the flags false.
 def ResetFlags(self):
  for hs in self.halfSpaceList:
   hs.SetFlag(False)

 # Get an actual half space from its index.
 def Get(self, index):
  return self.halfSpaceList[index]

 def __len__(self):
  return len(self.halfSpaceList)

#********************************************************************************************************

# Unions and intersections of half-spaces.
# Set() gives the empty set
# Set(-n) gives the universal set

class Set:

 # Create a set, optionally corresponding to a halfspace.
 def __init__(self, halfSpaceIndex = None):
  if halfSpaceIndex is None:
   self.expression = FALSE
   return
  if halfSpaceIndex < 0:
   self.expression = TRUE
   return
  self.expression = booleanAlgebra.Symbol(halfSpaceIndex)

 # I hope the following are self-explanatory.

 def Unite(self, set2):
  result = Set()
  result.expression = OR(self.expression, set2.expression)
  return result

 def Intersect(self, set2):
  result = Set()
  result.expression = AND(self.expression, set2.expression)
  return result

 def Subtract(self, set2):
  result = Set()
  result.expression = AND(self.expression, NOT(set2.expression))
  return result

 def Simplify(self):
  result = Set()
  result.expression = self.expression.simplify()
  return result

 def __str__(self):
  return str(self.expression)

 # Convert the set's expression to reverse polish. It's best (except maybe for debugging)
 # To use the result of Simplify() on a set for this.
 def toRP(self):
  # Start by splitting the expression into unique tags, which are either brackets, operators, or symbols.
  # The symbols are the half-space indices (i.e. positive integers).
  s = str(self.expression)
  s = re.split(r'(\d+)', s)

  # This splits chains like ')&~(' into their individual elements
  expression = []
  for e in s:
   if len(e) <= 0:
    pass
   elif len(e) == 1:
    expression.append(e)
   elif e[0].isdigit():
    expression.append(e)
   else:
    for d in e:
     expression.append(d)

  # Now remove complement operators and replace the halfspaces that they complement
  # with the index of that halfspace's actual complement.
  s = expression
  expression = []
  e = 0
  while e < len(s):
   tag = s[e]
   if tag == '~':
    e += 1
    tag = s[e]
    hs = halfSpaces.Get(int(tag)).complement
    expression.append(str(hs))
   else:
    expression.append(tag)
   e += 1

  # Now run through the expression converting it to RP
  rpExpression = []
  stack = []
  e = 0
  while e < len(expression):
   tag = expression[e]
   if tag[0].isdigit():
    rpExpression.append(tag)
   else:
    if len(stack) == 0 or stack[-1] == '(':
     stack.append(tag)
    elif tag == '(':
     stack.append(tag)
    elif tag == ')':
     s = stack.pop()
     while s != '(':
      rpExpression.append(s)
      s = stack.pop()
    else:
     stack.append(tag)
   e += 1
  while len(stack) > 0:
   rpExpression.append(stack.pop())
  return rpExpression


# The end result. Start with the universal set, as the first
# operation that will be done on it is intersections.

finalSet = Set(-1)


#*******************************************************************************************************

# Graphics
#----------

# A model consists of a list of triangles for display. The pattern can be rotated on the screen
# and zoomed using the mouse so it can be easily seen.

class Model:
 def __init__(self, triangles):
  self.triangles = triangles
  self.xangle = 0
  self.yangle = 0
  diag = np.subtract(positiveCorner[0], negativeCorner[0])
  self.distance = -3.0*(maths.sqrt(np.dot(diag, diag)))

 def MouseDrag(self, dx, dy):
  self.xangle += dx
  self.xangle %= 360
  self.yangle += dy
  self.yangle %= 360

 def MouseScroll(self, dy):
  self.distance += dy

 def Draw(self):
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  gluPerspective(60, 1, 0.1, 1000)
  glMatrixMode(GL_MODELVIEW)

  glEnable(GL_DEPTH_TEST)
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
  glEnable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_COLOR_MATERIAL)
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE )

  glLoadIdentity()
  glTranslatef(0, 0, self.distance)
  glRotatef(self.xangle, 0, 1, 0)
  glRotatef(self.yangle, -1, 0, 0)



  for triangle in self.triangles:
   glBegin(GL_POLYGON)
   for point in triangle.Points():
    glColor3f(triangle.colour[0], triangle.colour[1], triangle.colour[2])
    glNormal3f(triangle.normal[0],triangle.normal[1],triangle.normal[2])
    glVertex3f(point[0], point[1], point[2])
   glEnd()

  glFlush()


# The World is the container for models (see above) for display.

class World:
 def __init__(self):
  self.element = []

 def AddModel(self, model):
  self.element.append(model)

 def Draw(self):
  for obj in self.element:
   obj.Draw()

 def Setup(self):
  glEnable(GL_DEPTH_TEST)

 def MouseDrag(self, dx, dy):
  for obj in self.element:
   obj.MouseDrag(dx, dy)

 def MouseScroll(self, scroll_y):
  for obj in self.element:
   obj.MouseScroll(scroll_y)

def PutWorldInWindow(world, title):
 win = window.Window(fullscreen=False, vsync=True, resizable=False, height=600, width=600, caption = title)
 range = np.multiply(np.subtract(positiveCorner, negativeCorner), [2, 2, 5])
 negP = np.subtract(middle, range)
 posP = np.add(middle, range)

 @win.event
 def on_resize(width, height):
  glViewport(0, 0, width, height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(negP[0], posP[0], negP[1], posP[1], negP[2], posP[2])
  glMatrixMode(GL_MODELVIEW)
  return pyglet.event.EVENT_HANDLED

 @win.event
 def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
  world.MouseDrag(dx, dy)

 @win.event
 def on_mouse_scroll(x, y, scroll_x, scroll_y):
  world.MouseScroll(scroll_y)

 @win.event
 def on_draw():
  glClearColor(0.2, 0.2, 0.2, 0.8)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  world.Draw()

 world.Setup()


# Set up a window for a list of triangles
def MakeTriangleWindow(triangles, title):
 if graphics == 0:
  return;
 # They may be given a different colour after this function is called.
 triangles = copy.deepcopy(triangles)
 world = World()
 m = Model(triangles)
 world.AddModel(m)
 PutWorldInWindow(world, title)

#*************************************************************************************************

# The array triangles is assumed to represent the complete surface of a polyhedron.
# This uses the Jordan Curve Theorem to count intersections between the polyhedron
# and a straight line to the centroid of each triangle to decide if the triangle's
# normal is pointing from inside to outside. If it isn't the triangle is inverted.
# The classification is done by counting the number of times the line crosses the polyhedron's
# surface (not including each triangle of interest). If that's odd the line points
# roughly from inside to outside at the triangle. If then the inner product of the line's
# direction and the triangle's normal is negative, the triangle needs to be inverted.
# Similarly, if the count is even and the inner product is positive, the triangle also needs
# to be inverted. Otherwise it is left alone.

def AdjustAttitudes(triangles, hspaces):
 farPoint = np.multiply(positiveCorner, [random.uniform(20, 30), random.uniform(20, 30), random.uniform(20, 30)])
 for triangle in triangles:
  gradient = np.subtract(triangle.centroid, farPoint)
  crossCount = 0
  for t in triangles:
   if t is not triangle:
    hs = hspaces.Get(t.halfSpace)
    parameter = -(np.dot(farPoint, hs.normal) + hs.d)/np.dot(gradient, hs.normal)
    if parameter > 0 and parameter < 1:
     point = np.add(farPoint, np.multiply(gradient, parameter))
     if t.ContainsPoint(point):
      crossCount += 1
  if crossCount%2 == 1:
   if np.dot(triangle.normal, gradient) < 0:
    triangle.Invert()
  else:
   if np.dot(triangle.normal, gradient) > 0:
    triangle.Invert()

#*************************************************************************************************

# Recursive procedure for Woo's alternating sum of volumes algorithm.

def WooStep(trianglesAndPly):
 global finalSet

 # The triangles for which we want the next convex hull
 # Anything to do?
 triangles = trianglesAndPly[0]
 if len(triangles) <= 0:
  return

 # The original data with coordinates in etc.
 triangleFileData = trianglesAndPly[1]

 # The level of recursion
 level = trianglesAndPly[2]
 if level > maximumLevel:
  exit("WooStep(): maximum recursion exceeded.")

 title = "Level: " + str(level)
 if graphics > 0:
  MakeTriangleWindow(triangles, title + " original triangles")

 # The vertices of those triangles with no duplicates. These are indices of the
 # points in the original ply file.
 pointIndices = []

 # Add a halfspace for each triangle as long as it is not already in the halfspace list
 # Also gather all the vertices
 for triangle in triangles:
  halfSpaces.AddTriangle(triangle)
  for v in triangle.Vertices():
   pointIndices.append(v)

 # Remove duplicate point indices
 pointIndices = list(dict.fromkeys(pointIndices))

 # Give the triangles an attitude adjustment
 # TODO Is this correct? triangles may not be a closed surface
 AdjustAttitudes(triangles, halfSpaces)

 # The actual (x, y, z) coordinates of the points of which to compute the convex hull
 points = []
 for p in pointIndices:
  points.append(triangleFileData.Vertex(p))

 # Find the convex hull
 hull = ConvexHull(points)

 # Make a list of the hull triangles and their corresponding half spaces
 hullTriangles = []
 hullHalfSpaces = HalfSpaceList()
 for face in hull.simplices:
  vertices = []
  for f in face:
   # We need to record using the global indexing system
   vertices.append(pointIndices[f])
  triangle = Triangle(vertices, [0.5, 0.5, 1.0], triangleFileData)
  hullTriangles.append(triangle)
  hullHalfSpaces.AddTriangle(triangle)

 # hullTriangles is a complete polygon.
 AdjustAttitudes(hullTriangles, hullHalfSpaces)

 if graphics > 1:
  MakeTriangleWindow(hullTriangles, title + " CH")

 # Classify each triangle we came in with as being on
 # the hull or not. On means it coincides with a hull halfspace; not, not.
 newTriangles = []
 for triangle in triangles:
  halfSpace = HalfSpace(triangle)
  hs = hullHalfSpaces.LookUp(halfSpace)
  if hs < 0:
   triangle.SetColour([1, 0.5, 0.5])
   newTriangles.append(triangle)

 if graphics > 2:
  MakeTriangleWindow(newTriangles, title + " Inside triangles")

 #TODO
 #Split newTriangles into collections of connected triangles and call WooStep for each of them below

 # In what follows, why isn't uniqueHullHalfSpaces just hullHalfSpaces??

 #Add the hull to, or subtract it from, the model.

 uniqueHullHalfSpaces = []
 for t in hullTriangles:
  uniqueHullHalfSpaces.append(t.halfSpace)

 # Remove duplicates
 uniqueHullHalfSpaces = list(dict.fromkeys(uniqueHullHalfSpaces))

 hullSet = Set(-1)
 for h in uniqueHullHalfSpaces:
  bigListIndex = halfSpaces.AddSpace(hullHalfSpaces.Get(h))
  hullSet = hullSet.Intersect(Set(bigListIndex))
 if level == 0:
  finalSet = hullSet
 elif level%2 == 1:
  finalSet = finalSet.Subtract(hullSet)
 else:
  finalSet = finalSet.Unite(hullSet)

 # Simplfying as we go along increases execution time. But it also allows the simplification to
 # achieve greater reduction.
 finalSet = finalSet.Simplify()

 if len(newTriangles) > 0:
  # Carry on down...
  trianglesAndPly = (newTriangles, triangleFileData, level + 1)
  WooStep(trianglesAndPly)

#**************************************************************************************************

# Create the .set file

def ToFile(fileName, halfSpaces, set):
 file = open(fileName, "w")

 # The first things in the file are the bounding box
 file.write(str(negativeCorner) + "\n")
 file.write(str(positiveCorner) + "\n")

 # Get the RP expression to write out at the end
 rp = set.toRP()

 # Count the number of unique halfSpaces in the RP expression and set their flags
 count = 0
 for item in rp:
  if item[0].isdigit():
   hs = int(item)
   halfSpace = halfSpaces.Get(hs)
   if not halfSpace.flag:
    halfSpace.SetFlag(True)
    count += 1

 # Write out the count of the halfSpaces followed by the coefficients of each one
 file.write(str(count) + "\n")
 for c in range(len(halfSpaces)):
  hs = halfSpaces.Get(c)
  if hs.flag:
   file.write(str(c) + " " + str(hs.Coeficients()) + "\n")

 # Write out the reverse polish expression
 for item in rp:
  file.write(item + " ")
 file.write("\n")

 # Tidy up
 file.close()
 halfSpaces.ResetFlags()

#***********************************************************************************************************

# Run the conversion
#-------------------

halfSpaces = HalfSpaceList()

model = "STL2CSG-test-objects-woo-2"
place = "../../"
#fileName = '../../cube.ply'
#fileName = '../../two-disjoint-cubes.ply'
#fileName = '../../two-overlapping-cubes.ply'
#fileName = '../../hole-enclosed-in-cylinder.ply'
#fileName = '../../two-nonmanifold-cubes.ply'
#fileName = '../../two-nasty-nonmanifold-cubes.ply'
#fileName = '../../554.2-extruder-drive-pneumatic.ply'
#fileName = '../../cube-1-cylinder-1.ply'
#fileName = '../../STL2CSG-test-objects-woo-1.ply'
#fileName = '../../STL2CSG-test-objects-woo-2.ply'
#fileName = '../../STL2CSG-test-objects-cube-cylinder.ply'
#fileName = '../../STL2CSG-test-objects-cubePlusCylinder.ply'

# Load up the .ply file of triangles
triangleFileData = TriangleFileData(place+model+".ply")

# Find the bounding box
for v in range(triangleFileData.VertexCount()):
 vertex = triangleFileData.Vertex(v)
 middle = np.add(middle, vertex)
 positiveCorner = np.maximum(positiveCorner, vertex)
 negativeCorner = np.minimum(negativeCorner, vertex)
middle = np.multiply(middle, 1.0 / triangleFileData.VertexCount())

# Make a list of all the original triangles
originalTriangles = []
for t in range(triangleFileData.TriangleCount()):
 triangle = Triangle(triangleFileData.Triangle(t), [0.5, 1.0, 0.5], triangleFileData)
 originalTriangles.append(triangle)

# The recursive alternating sum of volumes function needs
# the list of triangles to work with, the .ply file data
# to get vertex coordinates from, and the current level
# or recursion (here at the start, 0).
trianglesAndPly = (originalTriangles, triangleFileData, 0)

# Run the algorithm
WooStep(trianglesAndPly)

# Save and maybe plot the results
finalSet = finalSet.Simplify()
ToFile(model+".set", halfSpaces, finalSet)
if graphics > 0:
 pyglet.app.run()

