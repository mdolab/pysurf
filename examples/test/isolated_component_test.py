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

    # Replace baseline surface coordinates if the user provided it
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

def computeCurveProjections(xyz, xyzd, allCoord, xyzProjb, allCoor=None):

    # Replace baseline coordinates if the user provided it
    if allCoor is not None:
        for curveName in cube.curves:
            cube.curves[curveName].update(allCoor[curveName])

    # Call projection algorithm
    xyzProj, tanProj, curveProjDict = cube.project_on_curve(xyz)
    
    # Call derivatives code in forward mode
    xyzProjd = cube.project_on_curve_d(xyz, xyzd, allCoord, xyzProj, curveProjDict)

    # Call derivatives code in backward mode
    xyzb, allCoorb = cube.project_on_curve_b(xyz, xyzProj, xyzProjb, curveProjDict)

    # Print results
    print
    print 'Original cube edge projection:'
    print xyzProj
    print normProj
    print

    # Return results
    return xyzProj, xyzProjd, tanProj, xyzb, allCoorb

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

allCoord = {}
for curveName in cube.curves:
    curveCoord = np.array(np.random.rand(cube.curves[curveName].coor.shape[0],cube.curves[curveName].coor.shape[1]),order='F')
    allCoord[curveName] = curveCoord/np.array([np.linalg.norm(curveCoord,axis=1)]).T

xyzProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
xyzProjb = xyzProjb/np.array([np.linalg.norm(xyzProjb,axis=1)]).T
normProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
normProjb = normProjb/np.array([np.linalg.norm(normProjb,axis=1)]).T

#############################
# SURFACE PROJECTION
#############################

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

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd*xyzb_AD)
dotProd = dotProd + np.sum(coord*coorb_AD)
dotProd = dotProd - np.sum(xyzProjd_AD*xyzProjb)
dotProd = dotProd - np.sum(normProjd_AD*normProjb)

# Compare directional derivatives
print 'SURFACE PROJECTION'
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

#############################
# CURVE PROJECTION
#############################

# Call projection function at the original point
xyzProj0, xyzProjd_AD, tanProj, xyzb_AD, allCoorb_AD = computeCurveProjections(xyz,
                                                                               xyzd, allCoord,
                                                                               xyzProjb)

# Create coordinates of perturbed points
allCoor = {}
for curveName in cube.curves:
    coor = cube.curves[curveName].coor
    allCoor[curveName] = coor + allCoord[curveName]*stepSize

# Call projection function at the perturbed point
xyzProj, xyzProjd, tanProj, dummy, dummy2 = computeCurveProjections(xyz+stepSize*xyzd,
                                                                    xyzd, allCoord,
                                                                    xyzProjb, allCoor)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0)/stepSize

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd*xyzb_AD)
for curveName in cube.curves:
    dotProd = dotProd + np.sum(allCoord[curveName]*allCoorb_AD[curveName])
dotProd = dotProd - np.sum(xyzProjd_AD*xyzProjb)

# Compare directional derivatives
print 'CURVE PROJECTION'
print ''
print 'xyzProjd_AD'
print xyzProjd_AD
print 'xyzProjd_FD'
print xyzProjd_FD

# Dot product test
print 'dot product (should be 0.0)'
print dotProd
