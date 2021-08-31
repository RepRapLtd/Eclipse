import sys
import random
import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
from plyfile import PlyData, PlyElement
from scipy.spatial import ConvexHull

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
  self.normal = np.multiply(self.normal, 1.0/s2)
  self.colour = RandomShade(colour)

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

def Run(world):
 win = window.Window(fullscreen=False, vsync=True, resizable=True, height=600, width=600)

 @win.event
 def on_resize(width, height):
  glViewport(0, 0, width, height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(-10, 20, -10, 20, -50, 50)
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

plyFileData = PlyData.read('../../../cube.ply')

positiveCorner = [-sys.float_info.max, -sys.float_info.max, -sys.float_info.max]
negativeCorner = [sys.float_info.max, sys.float_info.max, sys.float_info.max]

for v in plyFileData['vertex']:
 coords = []
 for c in v:
  coords.append(c)
 positiveCorner = np.maximum(positiveCorner, coords)
 negativeCorner = np.minimum(negativeCorner, coords)
 allPoints.append(coords)

print(positiveCorner, negativeCorner)

for f in plyFileData['face']:
 points = []
 for v in f[0]:
  coords = []
  for c in plyFileData['vertex'][v]:
   coords.append(c)
  points.append(coords)
 originalTriangles.append(Triangle(points, [0.5, 1.0, 0.5]))



hull = ConvexHull(allPoints)
for simplex in hull.simplices:
 points = []
 points.append(allPoints[simplex[0]])
 points.append(allPoints[simplex[1]])
 points.append(allPoints[simplex[2]])
 hullTriangles.append(Triangle(points, [1, 0.5, 0.5]))


#m = Model(hullTriangles)
m = Model(originalTriangles)
world.AddModel(m)
Run(world)

