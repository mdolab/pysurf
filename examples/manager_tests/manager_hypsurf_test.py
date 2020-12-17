# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle
import copy
import unittest


def forward_pass(manager, geom):

    # FORWARD PASS

    # Erase previous data
    manager.clear_all()

    # First we will create a new "intersection" curve based on the curves defined in the cylinder geometry
    ref_curve = geom.curves["source"]
    # curve = copy.deepcopy(ref_curve)
    curve = pysurf.TSurfCurve(ref_curve.name, ref_curve.coor, ref_curve.barsConn)
    curve.extra_data["parentGeoms"] = [geom.name, geom.name]
    manager.add_curve(curve)
    curveName = curve.name

    # Set marching options

    options1 = {
        "bc1": "splay",
        "bc2": "curve:bc1",
        "dStart": 0.01,
        "numLayers": 30,
        "extension": 2.8,
        "epsE0": 5.5,
        "theta": 0.0,
        "alphaP0": 0.25,
        "numSmoothingPasses": 3,
        "nuArea": 0.16,
        "numAreaPasses": 3,
        "sigmaSplay": 0.2,
        "cMax": 1.0,
        "ratioGuess": 10.0,
    }

    options2 = {
        "bc1": "splay",
        "bc2": "splay",
        "dStart": 0.01,
        "numLayers": 30,
        "extension": 2.8,
        "epsE0": 5.5,
        "theta": 0.0,
        "alphaP0": 0.25,
        "numSmoothingPasses": 3,
        "nuArea": 0.16,
        "numAreaPasses": 3,
        "sigmaSplay": 0.2,
        "cMax": 1.0,
        "ratioGuess": 10.0,
    }

    # Call meshing routine
    meshName = "mesh"
    meshNames = manager.march_intCurve_surfaceMesh(curveName, options1, options2, meshName)

    for meshName in meshNames:
        manager.meshGenerators[meshName].export_plot3d(meshName + ".xyz")

    return curveName, meshNames


# ===============================================


class HypsurfTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(HypsurfTest, self).__init__(*args, **kwargs)

    def test_hypsurf(self):
        # TESTING FUNCTION

        os.system("rm *.plt")

        # Load components
        geom = pysurf.TSurfGeometry("../inputs/cylinder.cgns")

        # We need to translate to avoid indetermination in derivative due to an projection exactly on the edge
        geom.curves["bc1"].translate(-0.00001, 0, 0)

        name1 = geom.name

        # Flip BC curve
        geom.curves["source"].flip()

        # Create manager object and add the geometry object to it
        manager0 = pysurf.Manager()
        manager0.add_geometry(geom)

        distTol = 1e-7

        # RUN FORWARD CODE
        curveName, meshNames = forward_pass(manager0, geom)

        # DERIVATIVE SEEDS

        # Generate random seeds
        coor1d, curveCoord = manager0.geoms[name1].set_randomADSeeds(mode="forward")

        intCoord = manager0.intCurves[curveName].set_randomADSeeds(mode="forward")

        meshb = []
        for meshName in meshNames:
            meshb.append(manager0.meshGenerators[meshName].meshObj.set_randomADSeeds(mode="reverse"))

        # FORWARD AD

        # Call AD code
        manager0.forwardAD()

        # Get relevant seeds
        meshd = []
        for meshName in meshNames:
            meshd.append(manager0.meshGenerators[meshName].meshObj.get_forwardADSeeds())

        # REVERSE AD

        # Call AD code
        manager0.reverseAD()

        # Get relevant seeds
        coor1b, curveCoorb = manager0.geoms[name1].get_reverseADSeeds()
        intCoorb = manager0.intCurves[curveName].get_reverseADSeeds()

        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(intCoorb * intCoord)
        dotProd = dotProd + np.sum(coor1b * coor1d)
        for curveName in curveCoord:
            dotProd = dotProd + np.sum(curveCoord[curveName] * curveCoorb[curveName])
        for ii in range(len(meshNames)):
            dotProd = dotProd - np.sum(meshb[ii] * meshd[ii])

        print("dotProd test")
        print(dotProd)
        np.testing.assert_almost_equal(dotProd, 0.0, decimal=13)

        # FINITE DIFFERENCE
        stepSize = 1e-7

        # Store initial coordinates of the seed curve
        seedCurveCoor = geom.curves["source"].get_points()

        # Perturb the geometries
        geom.update(geom.coor + stepSize * coor1d)
        for curveName in curveCoord:
            geom.curves[curveName].set_points(geom.curves[curveName].get_points() + stepSize * curveCoord[curveName])
        geom.curves["source"].set_points(seedCurveCoor + stepSize * intCoord)

        # Create new manager with perturbed components
        manager2 = pysurf.Manager()
        manager2.add_geometry(geom)

        # Execute code with the perturbed geometry
        curveName, meshNames = forward_pass(manager2, geom)

        # Get coordinates of the mesh nodes
        meshd_FD = []
        for meshName in meshNames:
            mesh0 = manager0.meshGenerators[meshName].meshObj.get_points()
            mesh = manager2.meshGenerators[meshName].meshObj.get_points()
            curr_meshd = (mesh - mesh0) / stepSize
            meshd_FD.append(curr_meshd)

        for ii in range(len(meshNames)):

            # Print results
            print("FD test")
            FD_results = np.max(np.abs(meshd[ii] - meshd_FD[ii]))
            print(FD_results)
            np.testing.assert_almost_equal(FD_results, 0.0, decimal=6)

        print("dotProd test")
        print(dotProd)
        np.testing.assert_almost_equal(dotProd, 0.0, decimal=13)


if __name__ == "__main__":
    unittest.main()
