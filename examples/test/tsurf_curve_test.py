"""
This script tests the projection algorithm for a TSurf curve.

"""

import pysurf
import numpy as np
import unittest

# TESTING FUNCTION
class TestCurveProjection(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestCurveProjection, self).__init__(*args, **kwargs)
        coor = np.array([[0.0, 0.0, 0.0],
                         [-1, -1, -1],
                         [-1, -1, -1],
                         [1.0, 0.0, 0.0],
                         [2.0, 0.0, 0.0],
                         [1.0, 1.0, 0.0]])

        self.curve = pysurf.tsurf_tools.create_curve_from_points(coor, 'test_curve')

    def test_orig_curve_projection(self):
        pts = np.array([[3.0, 0.0, 0.0],
                        [2.0, 1.0, 0.5]])
        projPts = pts.copy()

        self.curve.project(pts, xyzProj=projPts)

        np.testing.assert_almost_equal(projPts, np.array([[ 2., 0., 0. ],
                                                          [ 1.5, .5, 0.]]))

if __name__ == "__main__":
    unittest.main()
