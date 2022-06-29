import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os


class TestCurveIntersection(unittest.TestCase):

    N_PROCS = 1

    def setUp(self):
        # Set communicator
        comm = MPI.COMM_WORLD

        # Load cube and cylinder files
        baseDir = os.path.dirname(os.path.abspath(__file__))
        cubeFile = os.path.join(baseDir, "..", "input_files", "cube.cgns")
        cylinderFile = os.path.join(baseDir, "..", "input_files", "cylinder.cgns")

        # Create geometry objects
        self.comp1 = pysurf.TSurfGeometry(cubeFile, comm=comm)
        self.comp2 = pysurf.TSurfGeometry(cylinderFile, comm=comm)

    def test_curve_intersection(self):

        comp1 = self.comp1
        comp2 = self.comp2

        # Check the number of curves for each component
        self.assertEqual(len(comp1.curves.values()), 3)
        self.assertEqual(len(comp2.curves.values()), 3)

        # Call intersection function
        intersectionOrig = comp1.intersect(comp2, distTol=1e-7)

        # Check the number of intersection curves
        self.assertEqual(len(intersectionOrig), 27)

        # Split the original curves
        pysurf.tsurf_tools.split_curves(comp1.curves)
        pysurf.tsurf_tools.split_curves(comp2.curves)

        # Check the number of curves for each component after splitting
        self.assertEqual(len(comp1.curves.values()), 14)
        self.assertEqual(len(comp2.curves.values()), 5)

        # Call intersection function with the split curves
        intersectionSplit = comp1.intersect(comp2, distTol=1e-7)

        # Check that the number of intersection curves is the same as before
        self.assertEqual(len(intersectionSplit), 27)

    def test_intersection_deriv(self):

        # Define distance tolerance
        distTol = 1e-7

        # Load dummy components whose data  will be overwritten
        comp1 = self.comp1
        comp2 = self.comp2

        # Define component A
        comp1.coor = np.array([[-0.1, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [3.0, 1.0, 0.0]])
        comp1.triaConnF = np.array([[1, 2, 3], [2, 4, 3]])
        comp1.quadsConnF = np.zeros((0, 4))

        # Update ADT
        comp1.update()

        # Define component B
        comp2.coor = np.array([[0.0, 0.2, -0.2], [3.5, 0.2, -0.2], [0.0, 0.2, 0.2], [3.5, 0.2, 0.2]])
        comp2.triaConnF = np.array([[1, 2, 3], [2, 4, 3]])
        comp2.quadsConnF = np.zeros((0, 4))

        # Update ADT
        comp2.update()

        # === FORWARD PASS ===

        # Call intersection code
        intCurves = comp1.intersect(comp2, distTol)
        intCurve = intCurves[0]

        # Store coordinates of the initial curve
        intCoor0 = intCurve.get_points()

        # === FORWARD AD ===

        # Set random AD seeds
        coor1d, curveCoor1d = comp1.set_randomADSeeds(mode="forward")
        coor2d, curveCoor2d = comp2.set_randomADSeeds(mode="forward")
        intCoorb = intCurve.set_randomADSeeds(mode="reverse")

        # Compute derivatives in forward mode
        comp1.intersect_d(comp2, intCurve, distTol)

        # Get derivative seeds of the intersected curve
        intCoord = intCurve.get_forwardADSeeds()

        # === REVERSE AD ===

        # Compute derivatives in reverse mode
        comp1.intersect_b(comp2, intCurve, distTol, accumulateSeeds=False)

        # Get derivative seeds of the geometry components
        coor1b, curveCoor1b = comp1.get_reverseADSeeds()
        coor2b, curveCoor2b = comp2.get_reverseADSeeds()

        # === FINITE DIFFERENCES ===

        # Set step size
        stepSize = 1e-7

        # Change the geometry
        comp1.update(comp1.coor + stepSize * coor1d)
        for curveName in curveCoor1d:
            comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize * curveCoor1d[curveName])
        comp2.update(comp2.coor + stepSize * coor2d)
        for curveName in curveCoor2d:
            comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize * curveCoor2d[curveName])

        # Run the intersection code once again
        intCurves = comp1.intersect(comp2, distTol)
        intCurve = intCurves[0]

        # Get the new intersection coordinates
        intCoor = intCurve.get_points()

        # Compute derivatives with finite differences
        intCoord_FD = (intCoor - intCoor0) / stepSize

        # === DERIVATIVE TESTS ===

        # Compute dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(coor1d * coor1b)
        dotProd = dotProd + np.sum(coor2d * coor2b)
        dotProd = dotProd - np.sum(intCoord * intCoorb)
        for curveName in curveCoor1d:
            dotProd = dotProd - np.sum(curveCoor1b[curveName] * curveCoor1d[curveName])
        for curveName in curveCoor2d:
            dotProd = dotProd - np.sum(curveCoor2b[curveName] * curveCoor2d[curveName])
        np.testing.assert_allclose(dotProd, 0, atol=3e-16)

        # Compare AD to FD
        np.testing.assert_allclose(intCoord, intCoord_FD, rtol=1e-7)


if __name__ == "__main__":
    unittest.main()
