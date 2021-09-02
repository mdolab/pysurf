"""
OUTDATED
"""

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import pysurf
import unittest

np.random.seed(123)


class TestRemesh(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestRemesh, self).__init__(*args, **kwargs)

        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.2, 0.9], [0.1, 0.3, 0.7], [0.5, 0.6, 0.5], [0.8, 0.7, 0.4], [1.0, 0.9, 0.2]])

        # Create curve object
        self.curve = pysurf.tsurf_tools.create_curve_from_points(coor, "test", periodic=False)

        # Create simple curve and store within the class
        coor = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 1.0, 0.5], [0.8, 0.7, 0.4], [0.0, 0.0, 0.0]])

        # Create curve object
        self.periodic_curve = pysurf.tsurf_tools.create_curve_from_points(coor, "test", periodic=True)

    def test_linear_remesh_(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh()
        master_curve = np.array(
            [
                [0.0, 0.1894573, 0.4608581, 0.782837, 1.0],
                [0.2, 0.367093, 0.5706436, 0.694279, 0.9],
                [0.9, 0.6552714, 0.519571, 0.405721, 0.2],
            ]
        ).T
        np.testing.assert_almost_equal(master_curve, curve.coor)

    def test_linear_remesh_newNodes(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(nNewNodes=3)
        master_curve = np.array([[0.0, 0.4608581, 1.0], [0.2, 0.5706436, 0.9], [0.9, 0.519571, 0.2]]).T
        np.testing.assert_almost_equal(master_curve, curve.coor)

    def test_cosine_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(spacing="cosine")
        master_curve = np.array(
            [
                [0.0, 0.0873804, 0.4608581, 0.8764255, 1.0],
                [0.2, 0.2873804, 0.5706436, 0.7764255, 0.9],
                [0.9, 0.7252393, 0.519571, 0.3235745, 0.2],
            ]
        ).T
        np.testing.assert_almost_equal(master_curve, curve.coor)

    def test_hypTan_remesh(self):
        # Remesh and ensure values are as expected
        curve = self.curve.remesh(spacing="hypTan")
        master_curve = np.array(
            [
                [0.0, 0.0852627, 0.4608581, 0.8794204, 1.0],
                [0.2, 0.2852627, 0.5706436, 0.7794204, 0.9],
                [0.9, 0.7294747, 0.519571, 0.3205796, 0.2],
            ]
        ).T
        np.testing.assert_almost_equal(master_curve, curve.coor)

    # def test_tangent_remesh(self):
    #     # Remesh and ensure values are as expected
    #     curve = self.curve.remesh(spacing='tangent')
    #     master_curve = np.array([[ 0.       ,  0.1675572,  0.4170479,  0.702797 ,  0.9318736],
    #                              [ 0.2      ,  0.3506679,  0.5377859,  0.667599 ,  0.8318736],
    #                              [ 0.9      ,  0.6662214,  0.5414761,  0.432401 ,  0.2681264]])
    #     np.testing.assert_almost_equal(master_curve, curve.coor)

    def test_periodic_linear_remesh_(self):
        # Remesh and ensure values are as expected
        curve = self.periodic_curve.remesh()
        master_curve = np.array(
            [
                [0.0, 0.42402, 0.5707677, 0.5172996],
                [0.0, 0.42402, 0.9292323, 0.4526371],
                [0.0, 0.42402, 0.4764108, 0.2586498],
            ]
        ).T
        np.testing.assert_almost_equal(master_curve, curve.coor)

    # These functions have changed but remesh is checked by manager_tests
    # def test_linear_remesh_derivs(self):
    #     # Forward mode
    #
    #     # Get random input seeds
    #     coord = np.random.random_sample(self.curve.coor.shape)
    #     coord_copy = coord.copy()
    #
    #     # Get output seeds
    #     newCoord = pysurf.tsurf_tools._remesh_d(self.curve, coord)
    #
    #     # Backward mode
    #
    #     # Get random input seeds
    #     newCoorb = np.random.random_sample(newCoord.shape)
    #     newCoorb_copy = newCoorb.copy()
    #
    #     # Get output seeds
    #     coorb = pysurf.tsurf_tools._remesh_b(self.curve, newCoorb)
    #
    #     # Compute dot product and make sure it's equal to 0
    #     lhs = np.sum(coord_copy * coorb)
    #     rhs = np.sum(newCoorb_copy * newCoord)
    #
    #     np.testing.assert_almost_equal(rhs-lhs, 0.)
    #
    # def test_periodic_linear_remesh_derivs(self):
    #     # Forward mode
    #
    #     # Get random input seeds
    #     coord = np.random.random_sample(self.periodic_curve.coor.shape)
    #     coord_copy = coord.copy()
    #
    #     # Get output seeds
    #     newCoord = pysurf.tsurf_tools._remesh_d(self.periodic_curve, coord)
    #
    #     # Backward mode
    #
    #     # Get random input seeds
    #     newCoorb = np.random.random_sample(newCoord.shape)
    #     newCoorb_copy = newCoorb.copy()
    #
    #     # Get output seeds
    #     coorb = pysurf.tsurf_tools._remesh_b(self.periodic_curve, newCoorb)
    #
    #     # Compute dot product and make sure it's equal to 0
    #     lhs = np.sum(coord_copy * coorb)
    #     rhs = np.sum(newCoorb_copy * newCoord)
    #
    #     np.testing.assert_almost_equal(rhs-lhs, 0.)
    #
    # def test_cosine_remesh_derivs(self):
    #     # Forward mode
    #
    #     # Get random input seeds
    #     coord = np.random.random_sample(self.curve.coor.shape)
    #     coord_copy = coord.copy()
    #
    #     # Get output seeds
    #     newCoord = pysurf.tsurf_tools._remesh_d(self.curve, coord, spacing='cosine')
    #
    #     # Backward mode
    #
    #     # Get random input seeds
    #     newCoorb = np.random.random_sample(newCoord.shape)
    #     newCoorb_copy = newCoorb.copy()
    #
    #     # Get output seeds
    #     coorb = pysurf.tsurf_tools._remesh_b(self.curve, newCoorb, spacing='cosine')
    #
    #     # Compute dot product and make sure it's equal to 0
    #     lhs = np.sum(coord_copy * coorb)
    #     rhs = np.sum(newCoorb_copy * newCoord)
    #
    #     np.testing.assert_almost_equal(rhs-lhs, 0.)
    #
    # def test_hypTan_remesh_derivs(self):
    #     # Forward mode
    #
    #     # Get random input seeds
    #     coord = np.random.random_sample(self.curve.coor.shape)
    #     coord_copy = coord.copy()
    #
    #     # Get output seeds
    #     newCoord = pysurf.tsurf_tools._remesh_d(self.curve, coord, spacing='hypTan')
    #
    #     # Backward mode
    #
    #     # Get random input seeds
    #     newCoorb = np.random.random_sample(newCoord.shape)
    #     newCoorb_copy = newCoorb.copy()
    #
    #     # Get output seeds
    #     coorb = pysurf.tsurf_tools._remesh_b(self.curve, newCoorb, spacing='hypTan')
    #
    #     # Compute dot product and make sure it's equal to 0
    #     lhs = np.sum(coord_copy * coorb)
    #     rhs = np.sum(newCoorb_copy * newCoord)
    #
    #     np.testing.assert_almost_equal(rhs-lhs, 0.)
    #
    # def test_periodic_hypTan_remesh_derivs(self):
    #     # Forward mode
    #
    #     # Get random input seeds
    #     coord = np.random.random_sample(self.periodic_curve.coor.shape)
    #     coord_copy = coord.copy()
    #
    #     # Get output seeds
    #     newCoord = pysurf.tsurf_tools._remesh_d(self.periodic_curve, coord, spacing='hypTan')
    #
    #     # Backward mode
    #
    #     # Get random input seeds
    #     newCoorb = np.random.random_sample(newCoord.shape)
    #     newCoorb_copy = newCoorb.copy()
    #
    #     # Get output seeds
    #     coorb = pysurf.tsurf_tools._remesh_b(self.periodic_curve, newCoorb, spacing='hypTan')
    #
    #     # Compute dot product and make sure it's equal to 0
    #     lhs = np.sum(coord_copy * coorb)
    #     rhs = np.sum(newCoorb_copy * newCoord)
    #
    #     np.testing.assert_almost_equal(rhs-lhs, 0.)


if __name__ == "__main__":
    unittest.main()
