"""
This script tests the projection algorithm for a TSurf geometry object.

First it checks the projection on the original cube, then it translate the cube
3 units in the z-direction and again checks the projection.

"""

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest

# TESTING FUNCTION
class TestTSurfProjection(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestTSurfProjection, self).__init__(*args, **kwargs)
        self.cube = pysurf.TSurfGeometry('../inputs/cube.cgns', MPI.COMM_WORLD)
        self.pts = np.array([[.6, .5, 1.0],
                             [.6, .5, 0.1]], order='F')

    def test_orig_cube_projection(self):
        xyzProj, normProj, projDict = self.cube.project_on_surface(self.pts)

        print()
        print('Original cube projection:')
        print(xyzProj)
        print(normProj)
        print()

        np.testing.assert_almost_equal(xyzProj, np.array([[ 0.6, 0.5, 1.],
                                                          [ 0.6, 0.5, 0.]]))
        np.testing.assert_almost_equal(normProj, np.array([[ 0., 0.,  1.],
                                                           [ 0., 0., -1.]]))

        # Call derivative code
        ptsd = np.zeros(self.pts.shape)
        ptsd[0,0] = 1.0
        #xyzProjd, normProjd = self.cube.project_on_surface_d(self.pts, ptsd, xyzProj, normProj, projDict, np.zeros(self.cube.coor.shape))

    def test_orig_cube_edge_projection(self):
        xyzProj, normProj, curveProjDict = self.cube.project_on_curve(self.pts)

        print()
        print('Original cube edge projection:')
        print(xyzProj)
        print(normProj)
        print()

        np.testing.assert_almost_equal(xyzProj, np.array([[ 0.6, 0.70999998, 1.],
                                                          [ 0.55, 0.55, 0.]]))
        np.testing.assert_almost_equal(normProj, np.array([[ -1., 0.,  0.],
                                                           [ -.70710678, -.70710678, 0.]]))
    def test_mod_cube_projection(self):
        self.cube.translate(0, 0, 3)
        xyzProj, normProj, projDict = self.cube.project_on_surface(self.pts)
        self.cube.translate(0, 0, -3)

        print()
        print('Modified cube projection:')
        print(xyzProj)
        print(normProj)
        print()

        np.testing.assert_almost_equal(xyzProj, np.array([[ 0.6, 0.5, 3.],
                                                          [ 0.6, 0.5, 3.]]))
        np.testing.assert_almost_equal(normProj, np.array([[ 0., 0.,  -1.],
                                                           [ 0., 0., -1.]]))

    def test_mod_cube_edge_projection(self):
        self.cube.translate(0, 0, 3)
        xyzProj, normProj, curveProjDict = self.cube.project_on_curve(self.pts)
        self.cube.translate(0, 0, -3)

        print()
        print('Modified cube edge projection:')
        print(xyzProj)
        print(normProj)
        print()

        np.testing.assert_almost_equal(xyzProj, np.array([[ 0.55, 0.55, 3.],
                                                          [ 0.55, 0.55, 3.]]))
        np.testing.assert_almost_equal(normProj, np.array([[ -.707106778, -.70710678,  0.],
                                                           [ -.70710678, -.70710678, 0.]]))

if __name__ == "__main__":
    unittest.main()
