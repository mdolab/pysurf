import numpy as np
import unittest
from pysurf.tsurf_tools import create_curve_from_points


class TestRemesh(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.2, 0.9], [0.1, 0.3, 0.7], [0.5, 0.6, 0.5], [0.8, 0.7, 0.4], [1.0, 0.9, 0.2]])

        # Create curve object
        self.curve = create_curve_from_points(coor, "test", periodic=False)

        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 1.0, 0.5], [0.8, 0.7, 0.4], [0.0, 0.0, 0.0]])

        # Create curve object
        self.periodic_curve = create_curve_from_points(coor, "test", periodic=True)

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
        coor_d = self.curve.set_randomADSeeds(mode="forward")
        remeshedCoor_b = remeshedCurve.set_randomADSeeds(mode="reverse")

        # Forward AD
        self.curve.remesh_d(remeshedCurve, **optionsDict)

        # Store relevant seeds
        remeshedCoor_d = remeshedCurve.get_forwardADSeeds()

        # Reverse AD
        self.curve.remesh_b(remeshedCurve, **optionsDict)

        # Store relevant seeds
        coor_b = self.curve.get_reverseADSeeds()

        # Dot product test
        dotProd_LHS = np.sum(coor_d * coor_b)
        dotProd_RHS = np.sum(remeshedCoor_d * remeshedCoor_b)
        np.testing.assert_allclose(dotProd_LHS, dotProd_RHS, rtol=1e-15)

        # Define step sizes
        stepSize_FD = 1e-7
        stepSize_CS = 1e-200

        # Apply FD perturbation
        self.curve.set_points(self.curve.get_points() + stepSize_FD * coor_d)

        # Remesh the curve
        remeshedCurve_pert = self.curve.remesh(**optionsDict)

        # Compute FD derivative
        remeshedCoor0 = remeshedCurve.get_points()
        remeshedCoor_pert = remeshedCurve_pert.get_points()
        remeshedCoor_FD = (remeshedCoor_pert - remeshedCoor0) / stepSize_FD

        # Compare AD to FD
        np.testing.assert_allclose(remeshedCoor_d, remeshedCoor_FD, rtol=1e-7)

        # Create a complex curve object for CS test
        self.curve_CS = create_curve_from_points(self.curve.get_points(), "test", periodic=False, dtype=complex)

        # Apply CS perturbation
        self.curve_CS.set_points(self.curve_CS.get_points() + stepSize_CS * coor_d * 1j)

        # Remesh the complex curve
        remeshedCurve_pert = self.curve_CS.remesh(**optionsDict)

        # Compute CS derivative
        remeshedCoor_pert = remeshedCurve_pert.get_points()
        remeshedCoor_CS = np.imag(remeshedCoor_pert) / stepSize_CS

        # Compare AD to CS
        np.testing.assert_allclose(remeshedCoor_d, remeshedCoor_CS, rtol=1e-7)


if __name__ == "__main__":
    unittest.main()
