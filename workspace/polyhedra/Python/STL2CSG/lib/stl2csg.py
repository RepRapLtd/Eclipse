

import pyglet
from pyglet.gl import *
from pyglet import window
import numpy as np
from plyfile import PlyData, PlyElement
from scipy.spatial import ConvexHull




def vector(type, *args):
    return (type*len(args))(*args)

class model:
    def __init__(self, hull, points): #, colorMatrix, index):
        self.hull = hull
        self.points = points
        #self.colourMatrix = colorMatrix
        #self.index = index
        self.angle = 0

    def update(self):
        self.angle += 1
        self.angle %= 360


    def draw(self):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(self.angle, 1, 1, 0)

        glEnableClientState(GL_VERTEX_ARRAY)
        #glEnableClientState(GL_COLOR_ARRAY)

        for simplex in self.hull.simplices:
         glBegin(GL_POLYGON)
         for s in range(3):
          glVertex3f(self.points[simplex, s][0], self.points[simplex, s][1], self.points[simplex, s][2])
         glEnd()

        glFlush()

        #glVertexPointer(3, GL_FLOAT, 0, self.vertices)
        #glColorPointer(3, GL_FLOAT, 0, self.colourMatrix)
        #glDrawElements(GL_QUADS, len(self.index), GL_UNSIGNED_INT, self.index)

        glDisableClientState(GL_COLOR_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)



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




def Run(win, mWorld):

 @win.event
 def on_resize(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(-2, 2, -2, 2, -2, 2)
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



win = window.Window(fullscreen=False, vsync=True, resizable=True, height=600, width=600)
mWorld = world()


vertices = []
rng = np.random.default_rng()
points = rng.random((4, 3))
hull = ConvexHull(points)
m = model(hull, points)
mWorld.addModel(m)
Run(win, mWorld)


'''
import matplotlib.pyplot as plt
rng = np.random.default_rng()
points = rng.random((5, 3))
hull = ConvexHull(points)
for simplex in hull.simplices:
 print(points[simplex, 0], points[simplex, 1], points[simplex, 2])




  self.ply = PlyData.read(fileName)
  for f in self.ply['face']:
   if len(f[0]) != 3:
    print("Non-trangular face in " + fileName + "has " + str(len(f[0])) + " vertices.")
    '''
