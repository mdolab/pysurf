# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import copy

# TESTING FUNCTION

# Set communicator
comm = MPI.COMM_WORLD

# Define components
cube = pysurf.ADTComponent('../inputs/cube.cgns', comm)
cylinder = pysurf.ADTComponent('../inputs/cubeAndCylinder.cgns', comm)

# Call intersection function
pysurf.compute_intersections([cube, cylinder])

# Move the cube far away
cube.update(cube.coor + 5.0)

# Recompute intersection
pysurf.compute_intersections([cube, cylinder])
