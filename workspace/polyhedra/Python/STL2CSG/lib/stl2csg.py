import sys
import random
import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
from plyfile import PlyData, PlyElement
from scipy.spatial import ConvexHull

small = 0.000000001

def RandomShade(baseColour):
 colour = [random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1)]
 for c in range(3):
  colour[c] += baseColour[c]
  if colour[c] < 0:
   colour[c] = 0
  if colour[c] > 1:
   colour[c] = 1
 return colour


class Triangle:
 def __init__(self, points, colour):
  self.points = points
  self.normal = np.cross( np.subtract(points[1], points[0]), np.subtract(points[2], points[0]) )
  s2 = maths.sqrt(np.dot(self.normal, self.normal))
  if s2 < small:
   print("Triangle: corners are collinear.")
  self.normal = np.multiply(self.normal, 1.0/s2)
  self.colour = RandomShade(colour)
  self.centroid = np.multiply(np.add(np.add(points[0],points[1]), points[2]), 1.0/3.0)
  self.halfSpace = (-1, True)

#*************************************************************************************************************

class HalfSpace:
 def __init__(self, triangle):
  self.triangles = []
  self.triangles.append(triangle)
  self.normal = triangle.normal
  self.d = np.dot(self.normal, triangle.centroid)

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

 def AddTriangle(self, triangle, same):
  self.triangles.append((triangle, same))

 def __str__(self):
  return "{" + str(self.normal[0]) + "X + " + str(self.normal[1]) + "Y + " + \
   str(self.normal[2]) + "Z + " + str(self.d) + " <= 0}"

#********************************************************************************************************

class HalfSpaceList:
 def __init__(self):
  self.halfSpaceList = []

 def LookUp(self, halfSpace):
  for hs in range(len(self.halfSpaceList)):
   listHalfSpace = self.halfSpaceList[hs]
   if halfSpace == listHalfSpace:
    return (hs, True)
   if halfSpace.Opposite(listHalfSpace):
    return (hs, False)
  return (-1, True)

 def add(self, triangle):
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




#*******************************************************************************************************

class Model:
 def __init__(self, triangles):
  self.triangles = triangles
  self.angle = 0

 def Update(self):
  self.angle += 1
  self.angle %= 360


 def Draw(self):
  glMatrixMode(GL_MODELVIEW)
  glEnable(GL_DEPTH_TEST)
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
  glDisable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_COLOR_MATERIAL)
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE )

  glLoadIdentity()
  glRotatef(self.angle, 0.8, 0.6, 0)


  for triangle in self.triangles:
   glBegin(GL_POLYGON)
   for point in triangle.points:
    glColor3f(triangle.colour[0], triangle.colour[1], triangle.colour[2])
    glNormal3f(triangle.normal[0],triangle.normal[1],triangle.normal[2])
    glVertex3f(point[0], point[1], point[2])
   glEnd()

  glFlush()



#*************************************************************************************************************

class World:
 def __init__(self):
  self.element = []

 def Update(self, dt):
  for obj in self.element:
   obj.Update()

 def AddModel(self, model):
  self.element.append(model)

 def Draw(self):
  for obj in self.element:
   obj.Draw()

 def Setup(self):
  glEnable(GL_DEPTH_TEST)

def Run(world, negPoint, posPoint, centre):
 win = window.Window(fullscreen=False, vsync=True, resizable=True, height=600, width=600)
 range = np.multiply(np.subtract(posPoint, negPoint), [2, 2, 5])
 negP = np.subtract(centre, range)
 posP = np.add(centre, range)

 @win.event
 def on_resize(width, height):
  glViewport(0, 0, width, height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(negP[0], posP[0], negP[1], posP[1], negP[2], posP[2])
  glMatrixMode(GL_MODELVIEW)
  return pyglet.event.EVENT_HANDLED

 @win.event
 def on_draw():
  glClearColor(0.2, 0.2, 0.2, 0.8)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  world.Draw()

 pyglet.clock.schedule(world.Update)
 world.Setup()
 pyglet.app.run()

#*************************************************************************************************


world = World()

allPoints = []
hullTriangles = []
originalTriangles = []
halfSpaces = HalfSpaceList()

#plyFileData = PlyData.read('../../../cube.ply')
#plyFileData = PlyData.read('../../../two-disjoint-cubes.ply')
plyFileData = PlyData.read('../../../two-overlapping-cubes.ply')
#plyFileData = PlyData.read('../../../hole-enclosed-in-cylinder.ply')
#plyFileData = PlyData.read('../../../two-nonmanifold-cubes.ply')
#plyFileData = PlyData.read('../../../two-nasty-nonmanifold-cubes.ply')
#plyFileData = PlyData.read('../../../554.2-extruder-drive-pneumatic.ply')
#plyFileData = PlyData.read('../../../cube-1-cylinder-1.ply')

positiveCorner = [-sys.float_info.max, -sys.float_info.max, -sys.float_info.max]
negativeCorner = [sys.float_info.max, sys.float_info.max, sys.float_info.max]
centroid = [0, 0, 0]

for v in plyFileData['vertex']:
 coords = []
 for c in v:
  coords.append(c)
 centroid = np.add(centroid, coords)
 positiveCorner = np.maximum(positiveCorner, coords)
 negativeCorner = np.minimum(negativeCorner, coords)
 allPoints.append(coords)

centroid = np.multiply(centroid, 1.0/len(allPoints))

print(negativeCorner, positiveCorner, centroid)

for f in plyFileData['face']:
 points = []
 for v in f[0]:
  coords = []
  for c in plyFileData['vertex'][v]:
   coords.append(c)
  points.append(coords)
 triangle = Triangle(points, [0.5, 1.0, 0.5])
 originalTriangles.append(triangle)
 halfSpaces.add(triangle)

for h in halfSpaces.halfSpaceList:
 print(str(h))


hull = ConvexHull(allPoints)
for simplex in hull.simplices:
 points = []
 points.append(allPoints[simplex[0]])
 points.append(allPoints[simplex[1]])
 points.append(allPoints[simplex[2]])
 hullTriangles.append(Triangle(points, [1, 0.5, 0.5]))

for triangle in hullTriangles:
 hs = HalfSpace(triangle)
 print(halfSpaces.LookUp(hs))


#m = Model(hullTriangles)
m = Model(originalTriangles)
world.AddModel(m)
Run(world, negativeCorner, positiveCorner, centroid)

