# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import copy

# TESTING FUNCTION

# Set communicator
comm = MPI.COMM_WORLD

# Define components
#cube = pysurf.TSurfComponent('../inputs/cube.cgns', comm)
#cylinder = pysurf.TSurfComponent('../inputs/cylinder.cgns', comm)

bigCube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)
smallCube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)

#cube = pysurf.TSurfComponent('../inputs/cubeAndCylinder.cgns', ['geom'], comm)
#cylinder = pysurf.TSurfComponent('../inputs/cubeAndCylinder.cgns', ['cylinder'], comm)

# Reduce cube size
smallCube.update(smallCube.coor*0.5)

# Translate cube
coor = smallCube.coor
coor[0,:] = coor[0,:] + 0.77
coor[1,:] = coor[1,:] + 0.25
coor[2,:] = coor[2,:] + 0.25
smallCube.update(coor)

# Call intersection function
Intersections = pysurf.compute_intersections([bigCube, smallCube])

print Intersections[0].barsConn

'''
# Move the cube far away
cube.update(cube.coor + 5.0)

# Recompute intersection
pysurf.compute_intersections([cube, cylinder])
'''
