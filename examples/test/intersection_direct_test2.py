# This script is used to test intersection derivatives.
# I had to do it because the derivatives of the interior points of the intersection
# curve were twice as expected. I suspect this is due to the double propagation of derivative seeds.

# IMPORTS
from __future__ import division
import numpy as np
import pysurf
from mpi4py import MPI

# USER-DEFINED VARIABLES
pertNode = 1 # Node from component A to be perturbed
pertCoor = 2 # Coordinate to be perturbed [0,1,2]
stepSize = 1e-7

# EXECUTION

# Define communicator
comm = MPI.COMM_WORLD

# Define distance tolerance
distTol = 1e-7

# Load dummy components. Its data will be replaced later
TSurfGeometryA = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)
TSurfGeometryB = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)

# Define component A
TSurfGeometryA.coor = np.array([[0.0, 0.0, 0.0],
                                [3.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [3.0, 1.0, 0.0]],order='F').T
TSurfGeometryA.triaConn = np.array([[1,2,3],
                                    [2,4,3]],order='F').T
TSurfGeometryA.quadsConn = np.zeros((4,0),order='F')

# Define component B
TSurfGeometryB.coor = np.array([[0.0, 0.2, -0.2],
                                [3.0, 0.2, -0.2],
                                [0.0, 0.2, 0.2],
                                [3.0, 0.2, 0.2]],order='F').T
TSurfGeometryB.triaConn = np.array([[1,2,3],
                                    [2,4,3]],order='F').T
TSurfGeometryB.quadsConn = np.zeros((4,0),order='F')

# FINITE DIFFERENCING

# Perturb one of the variables to do finite differencing
TSurfGeometryA.coor[pertCoor,pertNode] = TSurfGeometryA.coor[pertCoor,pertNode] + stepSize

# Call Fortran code to find perturbed intersections
Intersections = pysurf.tsurf_tools._compute_pair_intersection(TSurfGeometryA, TSurfGeometryB, distTol, comm)

# Store perturbed coordinates
coorIntPert = np.array(Intersections[0].coor)

# Remove perturbation
TSurfGeometryA.coor[pertCoor,pertNode] = TSurfGeometryA.coor[pertCoor,pertNode] - stepSize

# Call Fortran code to find unperturbed intersections
Intersections = pysurf.tsurf_tools._compute_pair_intersection(TSurfGeometryA, TSurfGeometryB, distTol, comm)

# Store initial point
coorInt0 = np.array(Intersections[0].coor)

# Compute derivatives
coorIntdFD = (coorIntPert - coorInt0)/stepSize

# DERIVATIVES

# Set seeds for forward mode
coorAd = np.zeros(TSurfGeometryA.coor.shape,order='F')
coorAd[pertCoor,pertNode] = 1.0

coorBd = np.zeros(TSurfGeometryB.coor.shape,order='F')

# Compute derivatives in forward mode
coorIntd = pysurf.tsurf_tools._compute_pair_intersection_d(TSurfGeometryA, TSurfGeometryB, Intersections[0], coorAd, coorBd, distTol)

# Print results
print ''
print 'Number of intersections:', len(Intersections)
print Intersections[0].barsConn
print Intersections[0].coor
print coorIntd
print coorIntdFD
