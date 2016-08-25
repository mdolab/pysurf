from __future__ import division
import pysurf
import numpy as np


coor = np.array([[0.0, 0.0, 0.0],
                 [-1, -1, -1],
                 [-1, -1, -1],
                 [1.0, 0.0, 0.0],
                 [2.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0]], order='F').T

barsConn = np.array([[1,4],
                     [6,5],
                     [4,5],
                     [5,4],
                     [5,4],
                     [5,4]], order='F').T

name = 'test_curve'

curve = pysurf.tsurf_geometry.Curve(coor, barsConn, name)

print curve.coor.T
print curve.barsConn.T
