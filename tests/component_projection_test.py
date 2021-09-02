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
cube = pysurf.TSurfGeometry("../inputs/cube.cgns", MPI.COMM_WORLD)

# Define generic function to compute projections
def computeProjections(xyz, xyzd, coord, xyzProjb, normProjb, coor=None):

    # Replace baseline surface coordinates if the user provided it
    if coor is not None:
        cube.update(coor)

    # Call projection algorithm
    xyzProj, normProj, projDict = cube.project_on_surface(xyz)

    # Set derivative seeds
    cube.set_forwardADSeeds(coord=coord)

    # Call derivatives code in forward mode
    xyzProjd, normProjd = cube.project_on_surface_d(xyz, xyzd, xyzProj, normProj, projDict)

    # Call derivatives code in backward mode
    xyzb = cube.project_on_surface_b(xyz, xyzProj, xyzProjb, normProj, normProjb, projDict)

    coorb = cube.coorb

    # Print results
    print()
    print("Cube projection:")
    print(xyzProj)
    print(normProj)
    print()

    # Return results
    return xyzProj, xyzProjd, normProj, normProjd, xyzb, coorb


def computeCurveProjections(xyz, xyzd, allCoord, xyzProjb, tanProjb, allCoor=None):

    # Replace baseline coordinates if the user provided it
    if allCoor is not None:
        for curveName in cube.curves:
            cube.curves[curveName].set_points(allCoor[curveName])

    # Call projection algorithm
    xyzProj, tanProj, curveProjDict = cube.project_on_curve(xyz)

    # Set derivative seeds
    cube.set_forwardADSeeds(curveCoord=allCoord)

    # Call derivatives code in forward mode
    xyzProjd, tanProjd = cube.project_on_curve_d(xyz, xyzd, xyzProj, tanProj, curveProjDict)

    # Call derivatives code in backward mode
    xyzb = cube.project_on_curve_b(xyz, xyzProj, xyzProjb, tanProj, tanProjb, curveProjDict)

    allCoorb = {}
    for curve in cube.curves.itervalues():
        allCoorb[curve.name] = curve.get_reverseADSeeds()

    # Print results
    print()
    print("Cube edge projection:")
    print(xyzProj)
    print(tanProj)
    print()

    # Return results
    return xyzProj, xyzProjd, tanProj, tanProjd, xyzb, allCoorb


# BACK TO MAIN PROGRAM

# Define points
xyz = np.array([[0.7, 0.55, 0.1], [0.7, 0.55, 0.9]])

# Define derivatives and normalize them to a given step size
stepSize = 1e-7

xyzd = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))
xyzd = xyzd / np.sqrt(np.sum(xyzd ** 2))

coord = np.array(np.random.rand(cube.coor.shape[0], cube.coor.shape[1]))
coord = coord / np.sqrt(np.sum(coord ** 2))

allCoord = {}
for curveName in cube.curves:
    # Get curve coordinate size
    coor = cube.curves[curveName].get_points()
    curveCoord = np.array(np.random.rand(coor.shape[0], coor.shape[1]))
    allCoord[curveName] = curveCoord / np.sqrt(np.sum(curveCoord ** 2))

xyzProjb = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))

normProjb = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))

tanProjb = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))

#############################
# SURFACE PROJECTION
#############################

# Call projection function at the original point
xyzProj0, xyzProjd_AD, normProj0, normProjd_AD, xyzb_AD, coorb_AD = computeProjections(
    xyz, xyzd, coord, xyzProjb, normProjb
)

# Call projection function at the perturbed point
xyzProj, xyzProjd, normProj, normProjd, dummy, dummy2 = computeProjections(
    xyz + stepSize * xyzd, xyzd, coord, xyzProjb, normProjb, coor=cube.coor + stepSize * coord
)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0) / stepSize
normProjd_FD = (normProj - normProj0) / stepSize

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd * xyzb_AD)
dotProd = dotProd + np.sum(coord * coorb_AD)
dotProd = dotProd - np.sum(xyzProjd_AD * xyzProjb)
dotProd = dotProd - np.sum(normProjd_AD * normProjb)

# Compare directional derivatives
print("SURFACE PROJECTION")
print("")
print("xyzProjd_AD")
print(xyzProjd_AD)
print("xyzProjd_FD")
print(xyzProjd_FD)
print("normProjd_AD")
print(normProjd_AD)
print("normProjd_FD")
print(normProjd_FD)

# Get finite difference error
FD_error_projection = np.max(np.abs(xyzProjd_AD - xyzProjd_FD))
FD_error_normal = np.max(np.abs(normProjd_AD - normProjd_FD))

# Print results
print("")
print("#=================================================================#")
print("Finite difference error for projection (should be around 1e-7)")
print(FD_error_projection)
print("Finite difference error for normals (should be around 1e-7)")
print(FD_error_normal)
print("dot product (should be 0.0)")
print(dotProd)
print("#=================================================================#")
print("")

#############################
# CURVE PROJECTION
#############################

# Call projection function at the original point
xyzProj0, xyzProjd_AD, tanProj0, tanProjd_AD, xyzb_AD, allCoorb_AD = computeCurveProjections(
    xyz, xyzd, allCoord, xyzProjb, tanProjb
)

# Create coordinates of perturbed points
allCoor = {}
for curveName in cube.curves:
    coor = cube.curves[curveName].get_points()
    allCoor[curveName] = coor + allCoord[curveName] * stepSize

# Call projection function at the perturbed point
xyzProj, dummy0, tanProj, dummy1, dummy2, dummy3 = computeCurveProjections(
    xyz + stepSize * xyzd, xyzd, allCoord, xyzProjb, tanProjb, allCoor
)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0) / stepSize
tanProjd_FD = (tanProj - tanProj0) / stepSize

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd * xyzb_AD)
for curveName in cube.curves:
    dotProd = dotProd + np.sum(allCoord[curveName] * allCoorb_AD[curveName])
dotProd = dotProd - np.sum(xyzProjd_AD * xyzProjb)
dotProd = dotProd - np.sum(tanProjd_AD * tanProjb)

# Compare directional derivatives
print("CURVE PROJECTION")
print("")
print("xyzProjd_AD")
print(xyzProjd_AD)
print("xyzProjd_FD")
print(xyzProjd_FD)
print("tanProjd_AD")
print(tanProjd_AD)
print("tanProjd_FD")
print(tanProjd_FD)

# Get finite difference error
FD_error_projection = np.max(np.abs(xyzProjd_AD - xyzProjd_FD))
FD_error_tangent = np.max(np.abs(tanProjd_AD - tanProjd_FD))

# Print results
print("")
print("#=================================================================#")
print("Finite difference error for projection (should be around 1e-7)")
print(FD_error_projection)
print("Finite difference error for normals (should be around 1e-7)")
print(FD_error_tangent)
print("dot product (should be 0.0)")
print(dotProd)
print("#=================================================================#")
print("")
