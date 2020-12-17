# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# MESH PARAMETERS
numSkinNodes = 129
LE_spacing = 0.01
TE_spacing = 0.001

# TESTING FUNCTION

os.system("rm *.plt")

# Load components
comp1 = pysurf.TSurfGeometry("../../inputs/initial_full_wing_crm4.cgns", ["wing", "curve_le"])
comp2 = pysurf.TSurfGeometry("../../inputs/fuselage_crm4.cgns", ["fuse"])

name1 = comp1.name
name2 = comp2.name

# comp1.rename(name1)
# comp2.rename(name2)

# Load TE curves and append them to the wing component
curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves("curve_te_upp.plt_")
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves("curve_te_low.plt_")
curveName = curve_te_upp.keys()[0]
comp1.add_curve(curve_te_upp[curveName])
curveName = curve_te_low.keys()[0]
comp1.add_curve(curve_te_low[curveName])
comp1.curves[curveName].flip()

# Translate the wing
comp1.translate(0.0, 0.0, 0.0001)

# Create manager object and add the geometry objects to it
manager0 = pysurf.Manager()
manager0.add_geometry(comp1)
manager0.add_geometry(comp2)

distTol = 1e-7

# ======================================================
# FORWARD PASS


def forward_pass(manager):

    """
    This function will apply all geometry operations to the given manager.
    """

    # INTERSECT

    # Call intersection function
    intCurveNames = manager.intersect(distTol=distTol)
    intCurveName = intCurveNames[0]

    # SPLIT

    # Split curves based on TE and LE curves
    optionsDict = {
        "splittingCurves": [comp1.curves["curve_le"], comp1.curves["curve_te_upp"], comp1.curves["curve_te_low"]]
    }

    splitCurveNames = manager.split_intCurve(intCurveName, optionsDict, criteria="curve")

    # REMESH

    # Find the highest z-coordinate of the entire intersection (vertical position)
    maxZ = -99999
    for curve in splitCurveNames:
        curr_maxZ = np.max(manager.intCurves[curve].coor[:, 2])
        maxZ = max(maxZ, curr_maxZ)

    # Now we can identify and remesh each curve properly
    for curveName in splitCurveNames:

        # Get pointer to the curve object
        curve = manager.intCurves[curveName]

        # The trailing edge curve will have less nodes than the other ones
        if curve.numNodes < 20:

            # This is the trailing edge curve.
            # Just apply an uniform spacing
            optionsDict = {"nNewNodes": 9}
            TE_curveName = manager.remesh_intCurve(curveName, optionsDict)

        else:

            # We have an upper or lower skin curve.
            # First let's identify if the curve is defined from
            # LE to TE or vice-versa
            curveCoor = curve.get_points()
            deltaX = curveCoor[-1, 0] - curveCoor[0, 0]

            if deltaX > 0:
                LE_to_TE = True
            else:
                LE_to_TE = False

            # Compute the highest vertical coordinate of the curve
            curr_maxZ = np.max(curve.coor[:, 2])

            # Now we can determine if we have upper or lower skin
            if curr_maxZ < maxZ:

                # We are at the lower skin

                if LE_to_TE:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": TE_spacing,
                        "finalSpacing": LE_spacing,
                    }

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": LE_spacing,
                        "finalSpacing": TE_spacing,
                    }

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

            else:

                # We are at the upper skin

                if LE_to_TE:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": TE_spacing,
                        "finalSpacing": LE_spacing,
                    }

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": LE_spacing,
                        "finalSpacing": TE_spacing,
                    }

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

    # MERGE
    curveNames = [TE_curveName, LS_curveName, US_curveName]
    mergedCurveName = "intersection"
    manager.merge_intCurves(curveNames, mergedCurveName)

    # REORDER
    manager.intCurves[mergedCurveName].shift_end_nodes(criteria="maxX")

    # Export final curve
    manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName)

    # Flip the curve for marching if necessary
    mergedCurveCoor = manager.intCurves[mergedCurveName].get_points()
    deltaZ = mergedCurveCoor[1, 2] - mergedCurveCoor[0, 2]

    if deltaZ > 0:
        manager.intCurves[mergedCurveName].flip()

    # MARCH SURFACE MESHES
    meshName = "mesh"

    options_body = {
        "bc1": "continuous",
        "bc2": "continuous",
        "dStart": 0.01,
        "numLayers": 49,
        "extension": 1.3,
        "epsE0": 4.5,
        "theta": -0.5,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "nuArea": 0.16,
        "numAreaPasses": 20,
        "sigmaSplay": 0.3,
        "cMax": 5.0,
        "ratioGuess": 10.0,
    }

    options_wing = {
        "bc1": "curve:curve_te_upp",
        "bc2": "curve:curve_te_upp",
        "dStart": 0.01,
        "numLayers": 49,
        "extension": 1.4,
        "epsE0": 12.5,
        "theta": -0.8,
        "alphaP0": 0.25,
        "numSmoothingPasses": 4,
        "nuArea": 0.16,
        "numAreaPasses": 0,
        "sigmaSplay": 0.3,
        "cMax": 10.0,
        "ratioGuess": 10.0,
        "guideCurves": ["curve_te_low"],
        "remesh": True,
    }

    meshNames = manager.march_intCurve_surfaceMesh(
        mergedCurveName, options0=options_wing, options1=options_body, meshName=meshName
    )

    for meshName in meshNames:
        manager.meshGenerators[meshName].export_plot3d(meshName + ".xyz")

    return mergedCurveName


# END OF forward_pass
# ======================================================

# Call the forward pass function to the original manager
mergedCurveName = forward_pass(manager0)

# Export surface meshes of manager0
manager0.export_meshes("./surf_mesh_output")

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoor1d = manager0.geoms[name1].set_randomADSeeds(mode="forward")
coor2d, curveCoor2d = manager0.geoms[name2].set_randomADSeeds(mode="forward")

meshb = []
for meshGen in manager0.meshGenerators.itervalues():
    meshb.append(meshGen.meshObj.set_randomADSeeds(mode="reverse"))

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
meshd = []
for meshGen in manager0.meshGenerators.itervalues():
    meshd.append(meshGen.meshObj.get_forwardADSeeds())

# REVERSE AD

# Call AD code
manager0.reverseAD()

# Get relevant seeds
coor1b, curveCoor1b = manager0.geoms[name1].get_reverseADSeeds()
coor2b, curveCoor2b = manager0.geoms[name2].get_reverseADSeeds()

# Dot product test
dotProd = 0.0

for ii in range(len(meshd)):
    dotProd = dotProd + np.sum(meshd[ii] * meshb[ii])

dotProd = dotProd - np.sum(coor1b * coor1d)

for curveName in curveCoor1d:
    dotProd = dotProd - np.sum(curveCoor1d[curveName] * curveCoor1b[curveName])

dotProd = dotProd - np.sum(coor2b * coor2d)

for curveName in curveCoor2d:
    dotProd = dotProd - np.sum(curveCoor2d[curveName] * curveCoor2b[curveName])

print("dotProd test (this will be repeated at the end as well)")
print(dotProd)

# FINITE DIFFERENCE
stepSize = 1e-7

# Apply perturbations to the geometries
comp1.update(comp1.coor + stepSize * coor1d)

for curveName in comp1.curves:
    comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize * curveCoor1d[curveName])

comp2.update(comp2.coor + stepSize * coor2d)

for curveName in comp2.curves:
    comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize * curveCoor2d[curveName])

# Create new manager with perturbed components
manager1 = pysurf.Manager()
manager1.add_geometry(comp1)
manager1.add_geometry(comp2)

# Do the forward pass to the perturbed case
mergedCurveName = forward_pass(manager1)

# Get coordinates of the meshes in both cases to compute FD derivatives
meshd_FD = []
for meshName in manager0.meshGenerators:
    meshCoor0 = manager0.meshGenerators[meshName].meshObj.get_points()
    meshCoor1 = manager1.meshGenerators[meshName].meshObj.get_points()

    # Compute derivatives with FD
    meshCoord_FD = (meshCoor1 - meshCoor0) / stepSize

    # Append to the dictionary
    meshd_FD.append(meshCoord_FD)


def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt

    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation="none")
    plt.colorbar(im, orientation="horizontal")
    plt.show()


# Find the largest difference in derivatives
FD_error = 0.0
for ii in range(len(manager0.meshGenerators)):
    numNodes = manager0.meshGenerators[manager0.meshGenerators.keys()[ii]].numNodes
    curr_error = np.abs(meshd[ii] - meshd_FD[ii]).reshape((numNodes, -1, 3), order="F")

    view_mat(np.abs(curr_error[:, :, 0]))
    view_mat(np.abs(curr_error[:, :, 1]))
    view_mat(np.abs(curr_error[:, :, 2]))

    FD_error = max(FD_error, np.max(curr_error))

# Print results
print("dotProd test")
print(dotProd)
print("FD test")
print(FD_error)
