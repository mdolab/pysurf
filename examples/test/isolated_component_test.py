"""
This script tests the projection algorithm for a TSurf geometry object.

First it checks the projection on the original cube, then it translate the cube
3 units in the z-direction and again checks the projection.

"""

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np

# Load geometry
cube = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', MPI.COMM_WORLD)

# Define generic function to compute projections
def computeProjections(xyz, xyzd, coord, xyzProjb, normProjb, coor=None):

    # Replace baseline surface coordinates if the user provided none
    if coor is not None:
        cube.update(coor)

    # Call projection algorithm
    xyzProj, normProj, projDict = cube.project_on_surface(xyz)
    
    # Call derivatives code in forward mode
    xyzProjd, normProjd = cube.project_on_surface_d(xyz, xyzd, xyzProj, normProj, projDict, coord)

    # Call derivatives code in backward mode
    xyzb, coorb = cube.project_on_surface_b(xyz, xyzProj, xyzProjb, normProj, normProjb, projDict)

    # Print results
    print
    print 'Original cube projection:'
    print xyzProj
    print normProj
    print

    # Return results
    return xyzProj, xyzProjd, normProj, normProjd, xyzb, coorb

# BACK TO MAIN PROGRAM

# Define points
xyz = np.array([[.7, .55, 0.1],
                [.7, .55, 0.9]], order='F')

# Define derivatives and normalize them to a given step size
stepSize = 1e-7

xyzd = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
xyzd = xyzd/np.array([np.linalg.norm(xyzd,axis=1)]).T
coord = np.array(np.random.rand(cube.coor.shape[0],cube.coor.shape[1]),order='F')
coord = coord/np.array([np.linalg.norm(coord,axis=1)]).T

xyzProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
xyzProjb = xyzProjb/np.array([np.linalg.norm(xyzProjb,axis=1)]).T
normProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
normProjb = normProjb/np.array([np.linalg.norm(normProjb,axis=1)]).T

# Call projection function at the original point
xyzProj0, xyzProjd_AD, normProj0, normProjd_AD, xyzb_AD, coorb_AD = computeProjections(xyz,
                                                                                       xyzd, coord,
                                                                                       xyzProjb, normProjb)

# Call projection function at the perturbed point
xyzProj, xyzProjd, normProj, normProjd, dummy, dummy2 = computeProjections(xyz+stepSize*xyzd,
                                                                           xyzd, coord,
                                                                           xyzProjb, normProjb,
                                                                           coor=cube.coor+stepSize*coord)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0)/stepSize
normProjd_FD = (normProj - normProj0)/stepSize

# Compute normalize directional derivatives from AD
xyzProjd_AD = xyzProjd_AD
normProjd_AD = normProjd_AD

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd*xyzb_AD)
dotProd = dotProd + np.sum(coord*coorb_AD)
dotProd = dotProd - np.sum(xyzProjd_AD*xyzProjb)
dotProd = dotProd - np.sum(normProjd_AD*normProjb)

# Compare directional derivatives
print ''
print 'xyzProjd_AD'
print xyzProjd_AD
print 'xyzProjd_FD'
print xyzProjd_FD
print 'normProjd_AD'
print normProjd_AD
print 'normProjd_FD'
print normProjd_FD

# Dot product test
print 'dot product (should be 0.0)'
print dotProd
