# This script is used to test intersection derivatives.
# I had to do it because the derivatives of the interior points of the intersection
# curve were twice as expected. I suspect this is due to the double propagation of derivative seeds.

# IMPORTS
import numpy as np
import pysurf
from mpi4py import MPI

# USER-DEFINED VARIABLES
stepSize = 1e-7

# EXECUTION

# Define communicator
comm = MPI.COMM_WORLD

# Define distance tolerance
distTol = 1e-7

# Load dummy components. Its data will be replaced later on
comp1 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)
comp2 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)

# Define component A
comp1.coor = np.array([[-0.1, 0.0, 0.0],
                       [3.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [3.0, 1.0, 0.0]])
comp1.triaConnF = np.array([[1,2,3],
                            [2,4,3]])
comp1.quadsConnF = np.zeros((0,4))

# Update ADT
comp1.update()

# Define component B
comp2.coor = np.array([[0.0, 0.2, -0.2],
                       [3.5, 0.2, -0.2],
                       [0.0, 0.2, 0.2],
                       [3.5, 0.2, 0.2]])
comp2.triaConnF = np.array([[1,2,3],
                            [2,4,3]])
comp2.quadsConnF = np.zeros((0,4))

# Update ADT
comp2.update()

#===========================

# FORWARD PASS

# Call intersection code
intCurves = comp1.intersect(comp2, distTol)
intCurve = intCurves[0]
intName = intCurve.name

# Store coordinates of the initial curve
intCoor0 = intCurve.get_points()

# SETTING DERIVATIVE SEEDS

# Set seeds
coor1d, curveCoor1d = comp1.set_randomADSeeds(mode='forward')
coor2d, curveCoor2d = comp2.set_randomADSeeds(mode='forward')
intCoorb = intCurve.set_randomADSeeds(mode='reverse')

# FORWARD AD

# Compute derivatives in forward mode
comp1.intersect_d(comp2, intCurve, distTol)

# Get derivative seeds of the intersected curve
intCoord = intCurve.get_forwardADSeeds()

# REVERSE AD

# Make sure the curve has the same seeds on both sides
#comp1.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})
#comp2.set_reverseADSeeds(curveCoorb={intNames[0]:intCoorb})

# Compute derivatives in reverse mode
comp1.intersect_b(comp2, intCurve, distTol, accumulateSeeds=False)

# Get derivative seeds of the geometry components
coor1b, curveCoor1b = comp1.get_reverseADSeeds()
coor2b, curveCoor2b = comp2.get_reverseADSeeds()

# FINITE DIFFERENCES

# Set step size
stepSize = 1e-7

# Change the geometry
comp1.update(comp1.coor + stepSize*coor1d)
for curveName in curveCoor1d:
    comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize*curveCoor1d[curveName])
comp2.update(comp2.coor + stepSize*coor2d)
for curveName in curveCoor2d:
    comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize*curveCoor2d[curveName])

# Run the intersection code once again
intCurves = comp1.intersect(comp2, distTol)
intCurve = intCurves[0]

# Get the new intersection coordinates
intCoor = intCurve.get_points()

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
for curveName in curveCoor1d:
    dotProd = dotProd - np.sum(curveCoor1b[curveName]*curveCoor1d[curveName])
for curveName in curveCoor2d:
    dotProd = dotProd - np.sum(curveCoor2b[curveName]*curveCoor2d[curveName])

print ''
print 'dot product test (this should be around 1e-14):'
print dotProd
print ''


def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()
'''
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
intCoorb = np.zeros(intCurve.coor.shape,order='F')

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
