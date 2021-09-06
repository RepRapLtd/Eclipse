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
# OpenSCAD: https://openscad.org/
# FreeCAD: https://www.freecadweb.org/
# STL to PLY conversion free online: https://products.aspose.app/3d/conversion/stl-to-ply
#

import copy
import sys
import random
import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
from plyfile import PlyData
from scipy.spatial import ConvexHull

# Fudge factor
small = 0.000001

# What (if anything) to plot
#
# 0 - no graphics
# 1 - input triangles at each stage
# 2 - input and output triangles at each stage
# 3 - input and output triangles at each stage plus the convex hulls
graphics = 1

# Set operators
leaf = 0
union = 1
intersection = 2
complement = 3

# The bounding box and a point near the middle
positiveCorner = [-sys.float_info.max, -sys.float_info.max, -sys.float_info.max]
negativeCorner = [sys.float_info.max, sys.float_info.max, sys.float_info.max]
middle = [0, 0, 0]


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
  points = self.Points()
  self.normal = np.cross( np.subtract(points[1], points[0]), np.subtract(points[2], points[0]) )
  s2 = maths.sqrt(np.dot(self.normal, self.normal))
  if s2 < small:
   print("Triangle: vertices are collinear!")
  self.normal = np.multiply(self.normal, 1.0/s2)
  self.colour = RandomShade(colour)
  self.centroid = np.multiply(np.add(np.add(points[0],points[1]), points[2]), 1.0/3.0)
  self.halfSpace = (-1, True)

 # The three actual space (x, y, z) coordinates of the vertices
 def Points(self):
  return [self.ply.Vertex(self.vertices[0]), self.ply.Vertex(self.vertices[1]), self.ply.Vertex(self.vertices[2])]

 # The list of vertex indices.
 def Vertices(self):
  return self.vertices

 def SetColour(self, colour):
  self.colour = RandomShade(colour)

#*************************************************************************************************************

# A planar half space of the form Ax + By + Cz + D <= 0, where (A, B, C) is the
# plane's normal vector and D is the distance from it to the origin.
#
# It is the plane through a triangle, and it keeps a list of all the triangles
# that (nearly - global variable small) lie in it. The list is tuples (triangle, same)
# where same  is True means the normals coincide, False means they are exactly opposite.

class HalfSpace:
 def __init__(self, triangle):
  self.triangles = []
  self.triangles.append((triangle, True))
  self.normal = triangle.normal
  self.d = -np.dot(self.normal, triangle.centroid)

# Note: two coincident half spaces of opposite sense are not considered equal.
# To find such call Opposite below.

 def __eq__(self, halfSpace):
  if 1.0 - np.dot(halfSpace.normal, self.normal) > small:
   return False
  if abs(self.d - halfSpace.d) > small:
   return False
  return True

 def Opposite(self, halfSpace):
  if 1.0 + np.dot(halfSpace.normal, self.normal) > small:
   return False
  if abs(self.d + halfSpace.d) > small:
   return False
  return True

 # Keep a list of all triangles that lie in this halfspace
 def AddTriangle(self, triangle, same):
  self.triangles.append((triangle, same))

 def __str__(self):
  result = "{" + str(self.normal[0]) + "x + " + str(self.normal[1]) + "y + " +\
   str(self.normal[2]) + "z + " + str(self.d) + " <= 0}"
  result = result.replace("+ -", "- ") #Sigh!
  return result

#********************************************************************************************************

# Unions and intersections of half-spaces.

class Set:

 def __init__(self, halfSpaceIndex = None):
  if halfSpaceIndex is None:
   return
  self.o1 = halfSpaceIndex
  self.op = leaf

 def Union(self, set2):
  result = Set()
  result.o1 = self
  result.o2 = set2
  result.op = union
  return result

 def Intersection(self, set2):
  result = Set()
  result.o1 = self
  result.o2 = set2
  result.op = intersection
  return result

 def Complement(self):
  result = Set()
  result.o1 = self
  result.op = complement
  return result

 def Subtract(self, set2):
  result = set2.Complement()
  result = self.Intersection(result)
  return result

 # Convert a Set to a string, representing it in reverse polish
 def __str__(self):
  if self.op == leaf:
   return str(self.o1)
  elif self.op == union:
   return str(self.o1) + " " + str(self.o2) + " u"
  elif self.op == intersection:
   return str(self.o1) + " " + str(self.o2) + " ^"
  elif self.op == complement:
   return str(self.o1) + " -"
  sys.exit("Set(): operator not defined!")



#********************************************************************************************************

# A list of half spaces

class HalfSpaceList:
 def __init__(self):
  self.halfSpaceList = []

 # Is halfSpace in the list? If so return an index to it, together with a logical flag
 # indicating if it is the same sense or opposite. If not, return (-1, ...)
 def LookUp(self, halfSpace):
  for hs in range(len(self.halfSpaceList)):
   listHalfSpace = self.halfSpaceList[hs]
   if halfSpace == listHalfSpace:
    return (hs, True)
   if halfSpace.Opposite(listHalfSpace):
    return (hs, False)
  return (-1, True)

 # Add the halfspace corresponding to triangle to the list, unless that half space is
 # already in the list. Return the index of the half space in the list and a same/opposite flag.
 # Maintain the triangle's pointer to its corresponding half space, and the half space's list
 # of triangles in it.

 def Add(self, triangle):
  last = len(self.halfSpaceList)
  halfSpace = HalfSpace(triangle)
  inList = self.LookUp(halfSpace)
  hs = inList[0]
  if hs < 0:
   self.halfSpaceList.append(halfSpace)
   halfSpace.AddTriangle(triangle, True)
   triangle.halfSpace = (last, True)
   return (last, True)

  listHalfSpace = self.halfSpaceList[hs]
  if inList[1]:
   listHalfSpace.AddTriangle(triangle, True)
   triangle.halfSpace = (hs, True)
  else:
   listHalfSpace.AddTriangle(triangle, False)
   triangle.halfSpace = (hs, False)
  return inList

 def Get(self, index):
  return self.halfSpaceList[index]

# We want there to be just one global instance of this

halfSpaces = HalfSpaceList()

#*******************************************************************************************************

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
  glDisable(GL_LIGHTING)
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

#*************************************************************************************************************

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

# Recursive procedure for my variation on Woo's alternating sum of volumes algorithm.

def WooStep(pointsTrianglesAndPly):

 # The triangles for which we want the next convex hull
 # Anything to do?
 triangles = pointsTrianglesAndPly[1]
 if len(triangles) <= 0:
  return

 # The vertices of those triangles with no duplicates. These are indices of the
 # points in the original ply file.
 pointIndices = pointsTrianglesAndPly[0]

 # The original data with coordinates in etc.
 ply = pointsTrianglesAndPly[2]

 # The level of recursion
 level = pointsTrianglesAndPly[3]

 if graphics > 0:
  title = "Level: " + str(level)
  MakeTriangleWindow(triangles, title + " original triangles")

 # Add a halfspace for each triangle as long as it is not already in the halfspace list
 for triangle in triangles:
  halfSpaces.Add(triangle)

 # The actual (x, y, z) coordinates of the points to compute the convex hull of
 points = []
 for p in pointIndices:
  points.append(ply.Vertex(p))

 # Find the convex hull
 hull = ConvexHull(points)

 # Make a list of the hull triangles and their corresponding halfspaces
 # Note that any halfspaces that are already in the main list will not be duplicated
 # because of the way the halfSpaces.Add() function works.
 hullTriangles = []
 hullHalfSpaces = []
 for face in hull.simplices:
  vertices = []
  for f in face:
   vertices.append(pointIndices[f])
  triangle = Triangle(vertices, [0.5, 1.0, 0.5], ply)
  hullTriangles.append(triangle)
  hullHalfSpaces.append(halfSpaces.Add(triangle)[0])

 if graphics > 2:
  MakeTriangleWindow(hullTriangles, title + " CH")

 # Classify each triangle we came in with as being on
 # the hull or not. On means it coincides with a hull halfspace; not, not.
 newTriangles = []
 newPointIndices = []
 for triangle in triangles:
  hs = triangle.halfSpace[0]
  if hs < 0:
   sys.exit("WooStep(): Original triangle with no halfspace! Error!")
  else:
   alreadyThere = False
   for hhs in hullHalfSpaces:
    if hhs == hs:
     alreadyThere = True
   if not alreadyThere:
    # This triangle is not on the hull. So it needs to be considered at the next level of the recursion.
    triangle.SetColour([1, 0.5, 0.5])
    newTriangles.append(triangle)
    for v in triangle.Vertices():
     newPointIndices.append(v)

 # remove duplicate vertices
 newPointIndices = list(dict.fromkeys(newPointIndices))

 if graphics > 1:
  MakeTriangleWindow(newTriangles, title + " Inside triangles")

 if len(newTriangles) <= 0:
  print("Nothing left at recursion level " + str(level))
 else:
  # Carry on down...
  WooStep((newPointIndices, newTriangles, ply, level + 1))

#**************************************************************************************************

# Run the conversion

#fileName = '../../cube.ply'
#fileName = '../../two-disjoint-cubes.ply'
#fileName = '../../two-overlapping-cubes.ply'
#fileName = '../../hole-enclosed-in-cylinder.ply'
#fileName = '../../two-nonmanifold-cubes.ply'
#fileName = '../../two-nasty-nonmanifold-cubes.ply'
#fileName = '../../554.2-extruder-drive-pneumatic.ply'
#fileName = '../../cube-1-cylinder-1.ply'
fileName = '../../STL2CSG-test-objects-woo-1.ply'
#fileName = '../../STL2CSG-test-objects-woo-2.ply'
#fileName = '../../STL2CSG-test-objects-cube-cylinder.ply'
#fileName = '../../STL2CSG-test-objects-cubePlusCylinder.ply'
triangleFileData = TriangleFileData(fileName)

originalPointIndices = []
originalTriangles = []

for v in range(triangleFileData.VertexCount()):
 originalPointIndices.append(v)
 vertex = triangleFileData.Vertex(v)
 middle = np.add(middle, vertex)
 positiveCorner = np.maximum(positiveCorner, vertex)
 negativeCorner = np.minimum(negativeCorner, vertex)

middle = np.multiply(middle, 1.0 / triangleFileData.VertexCount())

print(negativeCorner, positiveCorner, middle)

for t in range(triangleFileData.TriangleCount()):
 triangle = Triangle(triangleFileData.Triangle(t), [0.5, 1.0, 0.5], triangleFileData)
 originalTriangles.append(triangle)

# Testing the Set class

s = Set(5).Subtract(Set(0).Union(Set(1)).Intersection(Set(2)))
print(str(s))

pointsTrianglesAndPly = (originalPointIndices, originalTriangles, triangleFileData, 0)

WooStep(pointsTrianglesAndPly)

if graphics > 0:
 pyglet.app.run()

