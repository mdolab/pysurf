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
comp1 = pysurf.TSurfGeometry("../../../inputs/initial_full_wing_crm4.cgns", ["wing"])
comp2 = pysurf.TSurfGeometry("../../../inputs/fuselage_crm4.cgns", ["fuse"])

# name1 = 'wing'
# name2 = 'body'

# comp1.rename(name1)
# comp2.rename(name2)

name1 = comp1.name
name2 = comp2.name

# Load TE curves and append them to the wing component
curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves("curve_te_upp.plt_")
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves("curve_te_low.plt_")
curveName = curve_te_upp.keys()[0]
comp1.add_curve(curve_te_upp[curveName])
curveName = curve_te_low.keys()[0]
comp1.add_curve(curve_te_low[curveName])
comp1.curves[curveName].flip()

# Move the wing to match intersection point
comp1.translate(0.0, 0.0, 1.5)

# Get the intersection curve from a file
curveDict = pysurf.tsurf_tools.read_tecplot_curves("intersection_000.plt_")
intCurveName = curveDict.keys()[0]
curve = curveDict[intCurveName]
curve.extra_data["parentGeoms"] = [name1, name2]

# Create manager object and add the geometry objects to it
manager0 = pysurf.Manager()
manager0.add_geometry(comp1)
manager0.add_geometry(comp2)
manager0.add_curve(curve)

distTol = 1e-7

numPasses = 0


def forward_pass(manager):

    # FORWARD PASS

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
        intCurveName, options0=options_wing, options1=options_body, meshName=meshName
    )

    for meshName in meshNames:
        manager.meshes[meshName].exportPlot3d(meshName + "_%03d" % numPasses + ".xyz")

    return meshNames


# ===============================================

# RUN FORWARD CODE
meshNames = forward_pass(manager0)
numPasses = numPasses + 1

# DERIVATIVE SEEDS

# Generate random seeds
manager0.geoms[name1].set_randomADSeeds(mode="forward")

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
    dotProd = dotProd + np.sum(curveCoor2d[curveName] * curveCoor2b[curevName])
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
