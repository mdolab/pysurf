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
cube = pysurf.TSurfGeometry('../inputs/cylinder.cgns', MPI.COMM_WORLD)

# Define generic function to compute projections
def computeProjections(pts, ptsd, coord, coor=None):

    # Replace baseline surface coordinates if the user provided none
    if coor is not None:
        cube.update(coor)

    # Call projection algorithm
    xyzProj, normProj, projDict = cube.project_on_surface(pts)
    
    # Call derivatives code in forward mode
    xyzProjd, normProjd = cube.project_on_surface_d(pts, ptsd, xyzProj, normProj, projDict, coord)

    # Print results
    print
    print 'Original cube projection:'
    print xyzProj
    print normProj
    print

    # Return results
    return xyzProj, xyzProjd, normProj, normProjd

# BACK TO MAIN PROGRAM

# Define points
pts = np.array([[.6, .5, 1.0],
                [.6, .5, 0.1]], order='F')

# Define derivatives and normalize them to a given step size
stepSize = 1e-8
ptsd = np.array(np.random.rand(pts.shape[0],pts.shape[1]),order='F')
ptsd = stepSize*ptsd/np.array([np.linalg.norm(ptsd,axis=1)]).T
coord = np.array(np.random.rand(cube.coor.shape[0],cube.coor.shape[1]),order='F')
coord = stepSize*coord/np.array([np.linalg.norm(coord,axis=1)]).T

# Call projection function at the original point
xyzProj0, xyzProjd_AD, normProj0, normProjd_AD = computeProjections(pts, ptsd, coord)

# Call projection function at the perturbed point
xyzProj, xyzProjd, normProj, normProjd = computeProjections(pts+ptsd, ptsd, coord, coor=cube.coor + coord)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0)
normProjd_FD = (normProj - normProj0)

# Compute normalize directional derivatives from AD
xyzProjd_AD = xyzProjd_AD
normProjd_AD = normProjd_AD

# Compare directional derivatives
print ''
print 'xyzProjd'
print xyzProjd_AD
print xyzProjd_FD
print 'normProjd'
print normProjd_AD
print normProjd_FD
