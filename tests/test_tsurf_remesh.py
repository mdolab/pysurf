import pysurf
import numpy as np
import unittest


class TestRemesh(unittest.TestCase):

    N_PROCS = 2

    def setUp(self):

        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.2, 0.9], [0.1, 0.3, 0.7], [0.5, 0.6, 0.5], [0.8, 0.7, 0.4], [1.0, 0.9, 0.2]])

        # Create curve object
        self.curve = pysurf.tsurf_tools.create_curve_from_points(coor, "test", periodic=False)

        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 1.0, 0.5], [0.8, 0.7, 0.4], [0.0, 0.0, 0.0]])

        # Create curve object
        self.periodic_curve = pysurf.tsurf_tools.create_curve_from_points(coor, "test", periodic=True)

        # Define relative tolerance
        self.rtol = 5e-7

    def test_linear_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh()
        ref_curve = np.array(
            [
                [0.0, 0.1894573, 0.4608581, 0.782837, 1.0],
                [0.2, 0.367093, 0.5706436, 0.694279, 0.9],
                [0.9, 0.6552714, 0.519571, 0.405721, 0.2],
            ]
        ).T
        np.testing.assert_allclose(ref_curve, curve.coor, rtol=self.rtol)

    def test_linear_remesh_newNodes(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(nNewNodes=3)
        ref_curve = np.array([[0.0, 0.4608581, 1.0], [0.2, 0.5706436, 0.9], [0.9, 0.519571, 0.2]]).T
        np.testing.assert_allclose(ref_curve, curve.coor, rtol=self.rtol)

    def test_cosine_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(spacing="cosine")
        ref_curve = np.array(
            [
                [0.0, 0.0873804, 0.4608581, 0.8764255, 1.0],
                [0.2, 0.2873804, 0.5706436, 0.7764255, 0.9],
                [0.9, 0.7252393, 0.519571, 0.3235745, 0.2],
            ]
        ).T
        np.testing.assert_allclose(ref_curve, curve.coor, rtol=self.rtol)

    def test_hypTan_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(spacing="hypTan")
        ref_curve = np.array(
            [
                [0.0, 0.0852627, 0.4608581, 0.8794204, 1.0],
                [0.2, 0.2852627, 0.5706436, 0.7794204, 0.9],
                [0.9, 0.7294747, 0.519571, 0.3205796, 0.2],
            ]
        ).T
        np.testing.assert_allclose(ref_curve, curve.coor, rtol=self.rtol)

    def test_periodic_linear_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.periodic_curve.remesh()
        ref_curve = np.array(
            [
                [0.0, 0.42402, 0.5707677, 0.5172996],
                [0.0, 0.42402, 0.9292323, 0.4526371],
                [0.0, 0.42402, 0.4764108, 0.2586498],
            ]
        ).T
        np.testing.assert_allclose(ref_curve, curve.coor, rtol=self.rtol)

    def test_remesh_deriv(self):
        # Remesh the curve
        optionsDict = {
            "nNewNodes": 13,
            "spacing": "linear",
            "initialSpacing": 0.005,
            "finalSpacing": 0.005,
        }
        remeshedCurve = self.curve.remesh(**optionsDict)

        # Assign derivative seeds
        coord = self.curve.set_randomADSeeds(mode="forward")
        remeshedCoorb = remeshedCurve.set_randomADSeeds(mode="reverse")

        # Adjust seeds
        coord[:, :] = 0.0
        coord[1, 2] = -1.0
        self.curve.set_forwardADSeeds(coord)

        # Forward AD
        self.curve.remesh_d(remeshedCurve, **optionsDict)

        # Store relevant seeds
        remeshedCoord = remeshedCurve.get_forwardADSeeds()

        # Reverse AD
        self.curve.remesh_b(remeshedCurve, **optionsDict)

        # Store relevant seeds
        coorb = self.curve.get_reverseADSeeds()

        # Compute dot product test
        dotProd = np.sum(coord * coorb) - np.sum(remeshedCoord * remeshedCoorb)
        np.testing.assert_allclose(dotProd, 0, atol=1e-15)

        # Define step size for finite difference test
        stepSize = 1e-7

        # Apply perturbations for finite difference tests
        self.curve.set_points(self.curve.get_points() + stepSize * coord)

        # Remesh the curve
        remeshedCurve2 = self.curve.remesh(**optionsDict)

        # Compare coordinates to compute derivatives
        remeshedCoor0 = remeshedCurve.get_points()
        remeshedCoorf = remeshedCurve2.get_points()
        remeshedCoord_FD = (remeshedCoorf - remeshedCoor0) / stepSize

        # Compare AD to FD
        np.testing.assert_allclose(remeshedCoord, remeshedCoord_FD, rtol=3e-5)


if __name__ == "__main__":
    unittest.main()
