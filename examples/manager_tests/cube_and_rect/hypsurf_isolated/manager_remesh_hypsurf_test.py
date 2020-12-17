# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

"""
This script tests just the surface mesh generation derivatives
starting from a given intersection curve.
This isolate the surface marching process from the geometry-handling steps.
"""

# MESH PARAMETERS
numSkinNodes = 129
LE_spacing = 0.001
TE_spacing = 0.01

# TESTING FUNCTION
os.system("rm *.plt")
os.system("rm *.xyz")

# Load components
comp1 = pysurf.TSurfGeometry("../../../inputs/cube_uns.cgns", ["geom"])
comp2 = pysurf.TSurfGeometry("../../../inputs/rect_uns.cgns", ["geom"])

name1 = comp1.name
name2 = comp2.name

# name1 = 'wing'
# name2 = 'body'

# comp1.rename(name1)
# comp2.rename(name2)

# ADDING GUIDE CURVES
# Create curve dictionary based on imported curves
# !!! Make sure to call `extract_curves.py` before running this script
curves = []
curves.append(pysurf.tsurf_tools.read_tecplot_curves("extracted_curve_000.plt_"))
curves.append(pysurf.tsurf_tools.read_tecplot_curves("extracted_curve_001.plt_"))
curves.append(pysurf.tsurf_tools.read_tecplot_curves("extracted_curve_002.plt_"))
curves.append(pysurf.tsurf_tools.read_tecplot_curves("extracted_curve_003.plt_"))

# Create an empty list in which we'll store the long (longitudinal) edges of
# the rect
long_curves = []
counter = 0

# Examine each of the imported curves
for ext_curve in curves:

    # Split these curves based on sharpness to get all edges of the rect
    split_curve = pysurf.tsurf_tools.split_curve_single(ext_curve[ext_curve.keys()[0]], "int", criteria="sharpness")

    # Loop over these split curves
    for name in split_curve:

        # Give the curves new names so they do not conflict with each other
        split_curve[name].name = "int_" + "{}".format(counter).zfill(3)
        counter += 1

        # Extract the curve points and compute the length
        pts = split_curve[name].get_points()
        length = pts[2, 0] - pts[2, -1]

        # Hardcoded logic here based on the length of edges
        if np.abs(length) > 1:

            # Flip the long curve if it's facing the 'wrong' direction
            if length > 0:
                split_curve[name].flip()

            # Add the long curve to the list
            long_curves.append(split_curve[name])

# Create a list of guideCurves based on the extracted long curves
# Note here we need the strings, not the curve objects.
# We use the same loop to add the guide curves to the rectangle object
guideCurves = []
for ext_curve in long_curves:
    print(ext_curve.name)
    comp2.add_curve(ext_curve)
    if ext_curve.name != "int_011":
        guideCurves.append(ext_curve.name)

# Rotate the rectangle in 5 degrees
comp2.rotate(5, 2)

distTol = 1e-7

# Set up integer to export different meshes
mesh_pass = 0

# Get the intersection curve from a file
curveDict = pysurf.tsurf_tools.read_tecplot_curves("intersection_shifted.plt_")
intCurveName = curveDict.keys()[0]
curve = curveDict[intCurveName]
curve.extra_data["parentGeoms"] = [name1, name2]

# Create manager object and add the geometry objects to it
manager0 = pysurf.Manager()
manager0.add_geometry(comp1)
manager0.add_geometry(comp2)
manager0.add_curve(curve)

distTol = 1e-7


def forward_pass(manager):

    # FORWARD PASS

    # REMESH
    optionsDict = {
        "nNewNodes": 41,
        "spacing": "linear",
        "initialSpacing": 0.005,
        "finalSpacing": 0.005,
    }
    remeshedCurveName = manager.remesh_intCurve(intCurveName, optionsDict)

    # MARCH SURFACE MESHES
    meshName = "mesh"

    options_rect = {
        "bc1": "curve:int_011",
        "bc2": "curve:int_011",
        "dStart": 0.03,
        "numLayers": 17,
        "extension": 3.5,
        "epsE0": 4.5,
        "theta": -0.5,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "nuArea": 0.16,
        "numAreaPasses": 20,
        "sigmaSplay": 0.3,
        "cMax": 10000.0,
        "ratioGuess": 1.5,
        "guideCurves": guideCurves,
    }

    options_cube = {
        "bc1": "continuous",
        "bc2": "continuous",
        "dStart": 0.02,
        "numLayers": 17,
        "extension": 2.5,
        "epsE0": 4.5,
        "theta": -0.5,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "nuArea": 0.16,
        "numAreaPasses": 20,
        "sigmaSplay": 0.3,
        "cMax": 10000.0,
        "ratioGuess": 1.5,
    }

    meshNames = manager.march_intCurve_surfaceMesh(
        intCurveName, options0=options_cube, options1=options_rect, meshName=meshName
    )

    # EXPORT
    for meshName in meshNames:
        manager.meshes[meshName].exportPlot3d(meshName + "_" + str(mesh_pass) + ".xyz")

    return meshNames


# ===============================================

# RUN FORWARD CODE
meshNames = forward_pass(manager0)
mesh_pass = mesh_pass + 1

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoor1d = manager0.geoms[name1].set_randomADSeeds(mode="forward")
coor2d, curveCoor2d = manager0.geoms[name2].set_randomADSeeds(mode="forward")

intCoord = manager0.intCurves[intCurveName].set_randomADSeeds(mode="forward")

meshb = []
for meshName in meshNames:
    meshb.append(manager0.meshes[meshName].set_randomADSeeds(mode="reverse"))

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
meshd = []
for meshName in meshNames:
    meshd.append(manager0.meshes[meshName].get_forwardADSeeds())

# REVERSE AD

# Call AD code
manager0.reverseAD()

# Get relevant seeds
coor1b, curveCoor1b = manager0.geoms[name1].get_reverseADSeeds()
coor2b, curveCoor2b = manager0.geoms[name2].get_reverseADSeeds()

intCoorb = manager0.intCurves[intCurveName].get_reverseADSeeds()

# Dot product test
dotProd = 0.0

print("int")
dotProd = dotProd + np.sum(intCoorb * intCoord)
print(dotProd)

print("coor1")
dotProd = dotProd + np.sum(coor1b * coor1d)
print(dotProd)

print("curves1")
for curveName in curveCoor1d:
    dotProd = dotProd + np.sum(curveCoor1d[curveName] * curveCoor1b[curveName])
    print(dotProd)

print("coor2")
dotProd = dotProd + np.sum(coor2b * coor2d)
print(dotProd)

print("curves2")
for curveName in curveCoor2d:
    dotProd = dotProd + np.sum(curveCoor2d[curveName] * curveCoor2b[curveName])
    print(dotProd)

print("mesh")
for ii in range(len(meshNames)):
    dotProd = dotProd - np.sum(meshb[ii] * meshd[ii])
    print(dotProd)

print("dotProd test")
print(dotProd)

# FINITE DIFFERENCE
stepSize = 1e-7

comp1.update(comp1.coor + stepSize * coor1d)

for curveName in comp1.curves:
    comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize * curveCoor1d[curveName])

comp2.update(comp2.coor + stepSize * coor2d)

for curveName in comp2.curves:
    comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize * curveCoor2d[curveName])

curve.set_points(curve.get_points() + stepSize * intCoord)

# Create new manager with perturbed components
manager2 = pysurf.Manager()
manager2.add_geometry(comp1)
manager2.add_geometry(comp2)
manager2.add_curve(curve)

# Execute code with the perturbed geometry
meshNames = forward_pass(manager2)

# Get coordinates of the mesh nodes
meshd_FD = []
for meshName in meshNames:
    mesh0 = manager0.meshes[meshName].mesh
    mesh = manager2.meshes[meshName].mesh
    curr_meshd = (mesh - mesh0) / stepSize
    meshd_FD.append(curr_meshd)

for ii in range(len(meshNames)):

    # Print results
    print("FD test for mesh", meshNames[ii])
    print(np.max(np.abs(meshd[ii] - meshd_FD[ii])))

print("dotProd test")
print(dotProd)
