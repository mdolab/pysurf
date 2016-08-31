# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np

# TESTING FUNCTION

# Define number of nodes in the circle
nNodes = 4

# Initialize ccor array
coor = np.zeros((3,nNodes), order='F')

# Initialize angles
theta = np.linspace(0, 2.0*np.pi, nNodes)

# Compute coordinates
for ii in range(nNodes):
    coor[0,ii] = np.cos(theta[ii])
    coor[1,ii] = np.sin(theta[ii])

# Generate connectivity
nNodes = coor.shape[1]
barsConn = np.zeros((2,nNodes-1))
barsConn[0,:] = range(1,nNodes)
barsConn[1,:] = range(2,nNodes+1)

# Make it periodic
barsConn[1,-1] = barsConn[0,0]

# Create curve object
curve = pysurf.TSurfCurve(coor, barsConn, 'test')

print curve.coor
print curve.barsConn

# Reorder curve
curve.shift_end_nodes(criteria='maxY')

print curve.coor
print curve.barsConn
