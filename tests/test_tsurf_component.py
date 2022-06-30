"""
This script tests the projection algorithm for a TSurf geometry object.

First it checks the projection on the original cube, then it translates the cube
3 units in the z-direction and again checks the projection.

"""

import os
import numpy as np
import unittest
from mpi4py import MPI
from pysurf import TSurfGeometry

# Set random seed for derivative tests
np.random.seed(123)


class TestTSurfProjection(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):
        baseDir = os.path.dirname(os.path.abspath(__file__))
        surfFile = os.path.join(baseDir, "..", "input_files", "cube.cgns")
        self.cube = TSurfGeometry(surfFile, comm=MPI.COMM_WORLD)
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
    def setUpDerivTest(self):
        cube = self.cube

        # Define points
        xyz = np.array([[0.7, 0.55, 0.1], [0.7, 0.55, 0.9]])

        # Define step sizes
        stepSize_FD = 1e-7
        stepSize_CS = 1e-200

        # Set random AD seeds
        xyz_d = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))
        xyz_d = xyz_d / np.sqrt(np.sum(xyz_d**2))

        coor_d = np.array(np.random.rand(cube.coor.shape[0], cube.coor.shape[1]))
        coor_d = coor_d / np.sqrt(np.sum(coor_d**2))

        allCoord = {}
        for curveName in cube.curves:
            # Get curve coordinate size
            coor = cube.curves[curveName].get_points()
            curveCoord = np.array(np.random.rand(coor.shape[0], coor.shape[1]))
            allCoord[curveName] = curveCoord / np.sqrt(np.sum(curveCoord**2))

        xyzProj_b = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))
        normProj_b = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))
        tanProj_b = np.array(np.random.rand(xyz.shape[0], xyz.shape[1]))

        return xyz, stepSize_FD, stepSize_CS, xyz_d, coor_d, allCoord, xyzProj_b, normProj_b, tanProj_b

    def computeSurfaceProjections_AD(self, xyz, xyz_d, coor_d, xyzProj_b, normProj_b):

        cube = self.cube

        # Call projection algorithm
        xyzProj, normProj, projDict = cube.project_on_surface(xyz)

        # Set derivative seeds
        cube.set_forwardADSeeds(coord=coor_d)

        # Call derivatives code in forward mode
        xyzProj_d, normProj_d = cube.project_on_surface_d(xyz, xyz_d, xyzProj, normProj, projDict)

        # Call derivatives code in backward mode
        xyz_b = cube.project_on_surface_b(xyz, xyzProj, xyzProj_b, normProj, normProj_b, projDict)
        coor_b = cube.coorb

        # Return results
        return xyzProj, xyzProj_d, normProj, normProj_d, xyz_b, coor_b

    def computeSurfaceProjections_pert(self, xyz_pert, coor_pert, CS=False):

        if CS:
            cube = self.cube_CS
        else:
            cube = self.cube

        # Replace the baseline surface coordinates
        cube.update(coor_pert)

        # Call projection algorithm
        xyzProj_pert, normProj_pert = cube.project_on_surface(xyz_pert)[0:2]

        # Return results
        return xyzProj_pert, normProj_pert

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

        # Call AD code in forward mode
        xyzProj_d, tanProj_d = cube.project_on_curve_d(xyz, xyzd, xyzProj, tanProj, curveProjDict)

        # Call AD code in backward mode
        xyz_b = cube.project_on_curve_b(xyz, xyzProj, xyzProjb, tanProj, tanProjb, curveProjDict)

        allCoor_b = {}
        for curve in cube.curves.values():
            allCoor_b[curve.name] = curve.get_reverseADSeeds()

        # Return results
        return xyzProj, xyzProj_d, tanProj, tanProj_d, xyz_b, allCoor_b

    def test_surface_projection_deriv(self):

        cube = self.cube
        xyz, stepSize_FD, stepSize_CS, xyz_d, coor_d, allCoord, xyzProj_b, normProj_b, tanProj_b = self.setUpDerivTest()

        # Compute projection and AD at the original point
        xyzProj0, xyzProj_d, normProj0, normProj_d, xyz_b, coor_b = self.computeSurfaceProjections_AD(
            xyz, xyz_d, coor_d, xyzProj_b, normProj_b
        )

        # Compute projection at the FD perturbed point
        xyz_pert = xyz + stepSize_FD * xyz_d
        coor_pert = cube.coor + stepSize_FD * coor_d
        xyzProj_pert, normProj_pert = self.computeSurfaceProjections_pert(xyz_pert, coor_pert)

        # Compute FD derivatives
        xyzProj_FD = (xyzProj_pert - xyzProj0) / stepSize_FD
        normProj_FD = (normProj_pert - normProj0) / stepSize_FD

        # Compute dot product test
        dotProd_LHS = np.sum(xyz_d * xyz_b) + np.sum(coor_d * coor_b)
        dotProd_RHS = np.sum(xyzProj_d * xyzProj_b) + np.sum(normProj_d * normProj_b)
        np.testing.assert_allclose(dotProd_LHS, dotProd_RHS, rtol=1e-15)

        # Compare AD to FD
        np.testing.assert_allclose(xyzProj_d, xyzProj_FD, rtol=1e-6)
        np.testing.assert_allclose(normProj_d, normProj_FD, rtol=3e-6)

    def test_curve_projection_deriv(self):

        cube = self.cube
        xyz, stepSize_FD, stepSize_CS, xyz_d, coor_d, allCoord, xyzProj_b, normProj_b, tanProj_b = self.setUpDerivTest()

        # Call projection function at the original point
        xyzProj0, xyzProj_d, tanProj0, tanProj_d, xyz_b, allCoor_b = self.computeCurveProjections(
            xyz, xyz_d, allCoord, xyzProj_b, tanProj_b
        )

        # Create coordinates of perturbed points
        allCoor = {}
        for curveName in cube.curves:
            coor = cube.curves[curveName].get_points()
            allCoor[curveName] = coor + allCoord[curveName] * stepSize_FD

        # Call projection function at the perturbed point
        xyzProj_pert, _, tanProj_pert = self.computeCurveProjections(
            xyz + stepSize_FD * xyz_d, xyz_d, allCoord, xyzProj_b, tanProj_b, allCoor
        )[0:3]

        # Compute directional derivatives with finite differencing
        xyzProj_FD = (xyzProj_pert - xyzProj0) / stepSize_FD
        tanProj_FD = (tanProj_pert - tanProj0) / stepSize_FD

        # Compute dot product test
        dotProd_LHS = np.sum(xyz_d * xyz_b)
        for curveName in cube.curves:
            dotProd_LHS += np.sum(allCoord[curveName] * allCoor_b[curveName])
        dotProd_RHS = np.sum(xyzProj_d * xyzProj_b) + np.sum(tanProj_d * tanProj_b)
        np.testing.assert_allclose(dotProd_LHS, dotProd_RHS, rtol=1e-15)

        # Compare AD to FD
        np.testing.assert_allclose(xyzProj_d, xyzProj_FD, rtol=3e-6)
        np.testing.assert_allclose(tanProj_d, tanProj_FD, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
