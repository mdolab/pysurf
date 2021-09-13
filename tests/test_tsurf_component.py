"""
This script tests the projection algorithm for a TSurf geometry object.

First it checks the projection on the original cube, then it translates the cube
3 units in the z-direction and again checks the projection.

"""

import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os

# Set random seed for derivative tests
np.random.seed(123)


class TestTSurfProjection(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):
        baseDir = os.path.dirname(os.path.abspath(__file__))
        surfFile = os.path.join(baseDir, "..", "input_files", "cube.cgns")
        self.cube = pysurf.TSurfGeometry(surfFile, MPI.COMM_WORLD)
        self.pts = np.array([[0.6, 0.5, 1.0], [0.6, 0.5, 0.1]], order="F")

    def test_orig_cube_projection(self):
        xyzProj, normProj, projDict = self.cube.project_on_surface(self.pts)

        np.testing.assert_allclose(xyzProj, np.array([[0.6, 0.5, 1.0], [0.6, 0.5, 0.0]]))
        np.testing.assert_allclose(normProj, np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]))

    def test_orig_cube_edge_projection(self):
        xyzProj, normProj, curveProjDict = self.cube.project_on_curve(self.pts)

        np.testing.assert_allclose(xyzProj, np.array([[0.6, 0.70999998, 1.0], [0.55, 0.55, 0.0]]))
        np.testing.assert_allclose(normProj, np.array([[-1.0, 0.0, 0.0], [-0.70710678, -0.70710678, 0.0]]))

    def test_mod_cube_projection(self):
        self.cube.translate(0, 0, 3)
        xyzProj, normProj, projDict = self.cube.project_on_surface(self.pts)
        self.cube.translate(0, 0, -3)

        np.testing.assert_allclose(xyzProj, np.array([[0.6, 0.5, 3.0], [0.6, 0.5, 3.0]]))
        np.testing.assert_allclose(normProj, np.array([[0.0, 0.0, -1.0], [0.0, 0.0, -1.0]]))

    def test_mod_cube_edge_projection(self):
        self.cube.translate(0, 0, 3)
        xyzProj, normProj, curveProjDict = self.cube.project_on_curve(self.pts)
        self.cube.translate(0, 0, -3)

        np.testing.assert_allclose(xyzProj, np.array([[0.55, 0.55, 3.0], [0.55, 0.55, 3.0]]))
        np.testing.assert_allclose(
            normProj, np.array([[-0.707106778, -0.70710678, 0.0], [-0.70710678, -0.70710678, 0.0]])
        )

    # Define helper functions for derivative testing
    def computeProjections(self, xyz, xyzd, coord, xyzProjb, normProjb, coor=None):

        cube = self.cube

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

        # Return results
        return xyzProj, xyzProjd, normProj, normProjd, xyzb, coorb

    def computeCurveProjections(self, xyz, xyzd, allCoord, xyzProjb, tanProjb, allCoor=None):

        cube = self.cube

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
        for curve in cube.curves.values():
            allCoorb[curve.name] = curve.get_reverseADSeeds()

        # Return results
        return xyzProj, xyzProjd, tanProj, tanProjd, xyzb, allCoorb

    def test_projection_deriv(self):

        cube = self.cube

        # Define points
        xyz = np.array([[0.7, 0.55, 0.1], [0.7, 0.55, 0.9]])

        # Define derivatives and normalize them to a given step size
        stepSize = 1e-7

        # Set random AD seeds
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

        # === SURFACE PROJECTION ===

        # Call projection function at the original point
        xyzProj0, xyzProjd_AD, normProj0, normProjd_AD, xyzb_AD, coorb_AD = self.computeProjections(
            xyz, xyzd, coord, xyzProjb, normProjb
        )

        # Call projection function at the perturbed point
        xyzProj, xyzProjd, normProj, normProjd, dummy, dummy2 = self.computeProjections(
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
        np.testing.assert_allclose(dotProd, 0, atol=1e-15)

        # Compare AD to FD
        np.testing.assert_allclose(xyzProjd_AD, xyzProjd_FD, rtol=1e-6)
        np.testing.assert_allclose(normProjd_AD, normProjd_FD, rtol=3e-6)

        # === CURVE PROJECTION ===

        # Call projection function at the original point
        xyzProj0, xyzProjd_AD, tanProj0, tanProjd_AD, xyzb_AD, allCoorb_AD = self.computeCurveProjections(
            xyz, xyzd, allCoord, xyzProjb, tanProjb
        )

        # Create coordinates of perturbed points
        allCoor = {}
        for curveName in cube.curves:
            coor = cube.curves[curveName].get_points()
            allCoor[curveName] = coor + allCoord[curveName] * stepSize

        # Call projection function at the perturbed point
        xyzProj, dummy0, tanProj, dummy1, dummy2, dummy3 = self.computeCurveProjections(
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
        np.testing.assert_allclose(dotProd, 0, atol=1e-15)

        # Compare AD to FD
        np.testing.assert_allclose(xyzProjd_AD, xyzProjd_FD, rtol=3e-6)
        np.testing.assert_allclose(tanProjd_AD, tanProjd_FD, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
