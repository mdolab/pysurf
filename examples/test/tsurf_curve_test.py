from __future__ import division
import pysurf
import numpy as np


coor = np.array([[0.0, 0.0, 0.0],
                 [-1, -1, -1],
                 [-1, -1, -1],
                 [1.0, 0.0, 0.0],
                 [2.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0]]).T

barsConn = np.array([[1,4],
                     [6,5],
                     [4,5]]).T

name = 'test_curve'

curve = pysurf.TSurfCurve(coor, barsConn, name)

print curve.coor
print curve.barsConn

pts = np.array([[3.0, 0.0, 0.0]])
projPts = pts.copy()


curve.project(pts, xyzProj=projPts)
print projPts
