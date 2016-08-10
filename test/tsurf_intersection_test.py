# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import copy

# TESTING FUNCTION

# Set communicator
comm = MPI.COMM_WORLD

# Define components
# cube = pysurf.TSurfComponent('../inputs/cube.cgns', comm)
# cylinder = pysurf.TSurfComponent('../inputs/cylinder.cgns', comm)

bigCube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)
smallCube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)

# cube = pysurf.TSurfComponent('../inputs/cubeAndCylinder.cgns', ['geom'], comm)
# cylinder = pysurf.TSurfComponent('../inputs/cubeAndCylinder.cgns', ['cylinder'], comm)


smallCube.scale(1.5)
smallCube.translate(.25, 0., 0.)
smallCube.rotate(30, 0)
smallCube.rotate(30, 1)
smallCube.rotate(30, 2)

# coor = cylinder.coor * .5
# coor[0,:] = coor[0,:] - 0.25
# coor[1,:] = coor[1,:] + 0.5
# coor[2,:] = coor[2,:] + 0.5
# cylinder.update(coor)


# Call intersection function
Intersections = pysurf.compute_intersections([bigCube, smallCube])
# Intersections = pysurf.compute_intersections([cube, cylinder])

# print 'FIRST'
# print Intersections[0].barsConn
# print 'second'
# print Intersections[1].barsConn
# print 'THIRD'
# print Intersections[2].barsConn

print ' Number of intersections found:\n          ', len(Intersections)

curve = pysurf.plot3d_interface.Curve(Intersections[0].coor, Intersections[0].barsConn)
curve.export_plot3d('curve0')

# curve = pysurf.plot3d_interface.Curve(Intersections[1].coor, Intersections[1].barsConn)
# curve.export_plot3d('curve1')

'''
# Move the cube far away
cube.update(cube.coor + 5.0)

# Recompute intersection
pysurf.compute_intersections([cube, cylinder])
'''
