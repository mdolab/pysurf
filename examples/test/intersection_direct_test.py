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
comp1 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)
comp2 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)

# Define component A
comp1.coor = np.array([[-0.1, 0.0, 0.0],
                       [3.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [3.0, 1.0, 0.0]],order='F').T
comp1.triaConn = np.array([[1,2,3],
                           [2,4,3]],order='F').T
comp1.quadsConn = np.zeros((4,0),order='F')

# Update ADT
comp1.update()

# Define component B
comp2.coor = np.array([[0.0, 0.2, -0.2],
                       [3.0, 0.2, -0.2],
                       [0.0, 0.2, 0.2],
                       [3.0, 0.2, 0.2]],order='F').T
comp2.triaConn = np.array([[1,2,3],
                           [2,4,3]],order='F').T
comp2.quadsConn = np.zeros((4,0),order='F')

# Update ADT
comp2.update()

#===========================

# FORWARD PASS

# Call intersection code
intNames = comp1.intersect(comp2, distTol)

# Store coordinates of the initial curve
intCoor0 = np.array(comp1.curves[intNames[0]].coor, order='F')

# SETTING DERIVATIVE SEEDS

# Set seeds
comp1.set_randomADSeeds()
comp2.set_randomADSeeds()

# Get seeds for future comparions
coor1d, _ = comp1.get_forwardADSeeds()
coor2d, _ = comp2.get_forwardADSeeds()
intCoorb = comp1.curves[intNames[0]].get_reverseADSeeds(clean=False)

# FORWARD AD

# Compute derivatives in forward mode
comp1.intersect_d(comp2, distTol)

# Get derivative seeds of the intersected curve
intCoord = comp1.curves[intNames[0]].get_forwardADSeeds()

# REVERSE AD

# Make sure the curve has the same seeds on both sides
#comp1.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})
#comp2.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})

# Compute derivatives in forward mode
comp1.intersect_b(comp2, distTol, accumulateSeeds=False)

# Get derivative seeds of the geometry components
coor1b, _ = comp1.get_reverseADSeeds()
coor2b, _ = comp2.get_reverseADSeeds()

# FINITE DIFFERENCES

# Set step size
stepSize = 1e-7

# Change the geometry
comp1.update(comp1.coor + stepSize*coor1d)
comp2.update(comp2.coor + stepSize*coor2d)

# Run the intersection code once again
intNames = comp1.intersect(comp2, distTol)

# Get the new intersection coordinates
intCoor = np.array(comp1.curves[intNames[0]].coor, order='F')

# Compute derivatives with finite differences
intCoord_FD = (intCoor - intCoor0)/stepSize

# FINITE DIFFERENCE TEST

# Compute the maximum difference between derivatives
maxError = np.max(np.abs(intCoord - intCoord_FD))

print ''
print 'FD test (this should be around',stepSize,'):'
print maxError
print ''

# DOT PRODUCT TEST
dotProd = 0.0
dotProd = dotProd + np.sum(coor1d*coor1b)
dotProd = dotProd + np.sum(coor2d*coor2b)
dotProd = dotProd - np.sum(intCoord*intCoorb)

print ''
print 'dot product test (this should be around 1e-14):'
print dotProd
print ''


'''
def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

# BUILDING THE JACOBIAN

nNodes1 = comp1.coor.shape[1]
nNodes2 = comp2.coor.shape[1]
nNodesInt = intCoor0.shape[1]

print 'nNodes1:',nNodes1
print 'nNodes2:',nNodes2
print 'nNodesInt:',nNodesInt

Jac_fwd = np.zeros((3*(nNodes1+nNodes2),3*nNodesInt))
Jac_rev = np.zeros((3*(nNodes1+nNodes2),3*nNodesInt))

# Initialize arrays of derivatives
coor1d = np.zeros(comp1.coor.shape, order='F')
coor2d = np.zeros(comp2.coor.shape, order='F')
intCoorb = np.zeros(comp1.curves[intNames[0]].coor.shape,order='F')

# FORWARD AD

for ii in range(nNodes1):
    for jj in range(3):
        coor1d[:,:] = 0.0
        coor1d[jj,ii] = 1.0
        coor2d[:,:] = 0.0

        comp1.set_forwardADSeeds(coord=coor1d)
        comp2.set_forwardADSeeds(coord=coor2d)

        # Compute derivatives in forward mode
        comp1.intersect_d(comp2, distTol)
        
        # Get derivative seeds of the intersected curve
        intCoord = comp1.curves[intNames[0]].get_forwardADSeeds()

        Jac_fwd[3*ii+jj,::3] = intCoord[0,:]
        Jac_fwd[3*ii+jj,1::3] = intCoord[1,:]
        Jac_fwd[3*ii+jj,2::3] = intCoord[2,:]

for ii in range(nNodes2):
    for jj in range(3):
        coor1d[:,:] = 0.0
        coor2d[:,:] = 0.0
        coor2d[jj,ii] = 1.0

        comp1.set_forwardADSeeds(coord=coor1d)
        comp2.set_forwardADSeeds(coord=coor2d)

        # Compute derivatives in forward mode
        comp1.intersect_d(comp2, distTol)
        
        # Get derivative seeds of the intersected curve
        intCoord = comp1.curves[intNames[0]].get_forwardADSeeds()

        Jac_fwd[3*(ii+nNodes1)+jj,::3] = intCoord[0,:]
        Jac_fwd[3*(ii+nNodes1)+jj,1::3] = intCoord[1,:]
        Jac_fwd[3*(ii+nNodes1)+jj,2::3] = intCoord[2,:]

# REVERSE AD

for ii in range(nNodesInt):
    for jj in range(3):
        intCoorb[:,:] = 0.0
        intCoorb[jj,ii] = 1.0

        comp1.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})
        comp2.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})

        # Compute derivatives in reverse mode
        comp1.intersect_b(comp2, distTol, accumulateSeeds=False)

        # Get derivative seeds
        coor1b,_ = comp1.get_reverseADSeeds()
        coor2b,_ = comp2.get_reverseADSeeds()

        Jac_rev[:3*nNodes1:3,3*ii+jj] = coor1b[0,:]
        Jac_rev[1:3*nNodes1:3,3*ii+jj] = coor1b[1,:]
        Jac_rev[2:3*nNodes1:3,3*ii+jj] = coor1b[2,:]

        Jac_rev[3*nNodes1:3*(nNodes1+nNodes2):3,3*ii+jj] = coor2b[0,:]
        Jac_rev[3*nNodes1+1:3*(nNodes1+nNodes2):3,3*ii+jj] = coor2b[1,:]
        Jac_rev[3*nNodes1+2:3*(nNodes1+nNodes2):3,3*ii+jj] = coor2b[2,:]

print ''
print 'Maximum difference between Jacobians:'
print np.max(np.abs(Jac_fwd-Jac_rev))
print ''

view_mat(Jac_fwd)

view_mat(Jac_rev)

view_mat(Jac_fwd-Jac_rev)
'''
