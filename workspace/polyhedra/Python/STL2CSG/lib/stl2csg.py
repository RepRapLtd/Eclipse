import random

import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
from plyfile import PlyData, PlyElement
from scipy.spatial import ConvexHull

class Triangle:
 def __init__(self, points):
  self.points = points
  self.normal = np.cross( np.subtract(points[1], points[0]), np.subtract(points[2], points[0]) )
  s2 = maths.sqrt(np.dot(self.normal, self.normal))
  self.normal = np.multiply(self.normal, 1.0/s2)
  self.colour = (random.random(), random.random(), random.random())



def vector(type, *args):
    return (type*len(args))(*args)

class model:
    def __init__(self, triangles):
        self.triangles = triangles
        self.angle = 0

    def update(self):
        self.angle += 1
        self.angle %= 360


    def draw(self):
        glMatrixMode(GL_MODELVIEW)

        glEnable(GL_DEPTH_TEST)

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

        glEnable(GL_LIGHTING)
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





class world:
    def __init__(self):
        self.element = []

    def update(self, dt):
        for obj in self.element:
            obj.update()

    def addModel(self, model):
        self.element.append(model)

    def draw(self):
        for obj in self.element:
            obj.draw()


def setup():
    # look for GL_DEPTH_BUFFER_BIT
    glEnable(GL_DEPTH_TEST)




def Run(mWorld):
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
    mWorld.draw()

 pyglet.clock.schedule(mWorld.update)
 setup()
 pyglet.app.run()

#*************************************************************************************************


mWorld = world()


triangles = []

'''
rng = np.random.default_rng()
allPoints = rng.random((4, 3))
hull = ConvexHull(allPoints)
for simplex in hull.simplices:
 points = []
 points.append(allPoints[simplex, 0])
 points.append(allPoints[simplex, 1])
 points.append(allPoints[simplex, 2])
 triangles.append(Triangle(points))
'''

ply = PlyData.read('../../../cube.ply')
for f in ply['face']:
 points = []
 for v in f[0]:
  coords = []
  for c in ply['vertex'][v]:
   coords.append(c)
  points.append(coords)
 triangles.append(Triangle(points))

m = model(triangles)
mWorld.addModel(m)
Run(mWorld)

