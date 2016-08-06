# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np

# TESTING FUNCTION

cube = pysurf.TSurfComponent('../inputs/cube.cgns', MPI.COMM_WORLD)

pts = np.array([[.6, .5, 0.1]], order='F')

xyzProj, normProj = cube.project_on_surface(pts)

print 'Original cube projection:'
print xyzProj
print normProj

xyzProj, normProj = cube.project_on_curve(pts)

print 'Original cube edge projection:'
print xyzProj
print normProj

# Change cube and update points
coor = cube.coor.copy()
coor[2,:] = coor[2,:] + 3.0
cube.update(coor)

# Repeat projections

xyzProj, normProj = cube.project_on_surface(pts)

print 'Modified cube projection:'
print xyzProj
print normProj

xyzProj, normProj = cube.project_on_curve(pts)

print 'Modified cube edge projection:'
print xyzProj
print normProj
