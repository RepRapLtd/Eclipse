
import numpy as np
from plyfile import PlyData, PlyElement
'''
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
rng = np.random.default_rng()
points = rng.random((30, 2))
hull = ConvexHull(points)
plt.plot(points[:,0], points[:,1], 'o')
for simplex in hull.simplices:
 plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.show()
'''

def ReadPLY(fileName):
 ply = PlyData.read(fileName)
 for f in ply['face']:
  if len(f[0]) != 3:
   print("Non-trangular face in " + fileName + "has " + str(len(f[0])) + " edges.")
 return ply

ply = ReadPLY('../../../cube.ply')
