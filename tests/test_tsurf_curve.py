"""
This script tests the projection algorithm for a TSurf curve.

"""

import pysurf
import numpy as np
import unittest


class TestCurveProjection(unittest.TestCase):
    N_PROCS = 1

    def test_orig_curve_projection(self):
        coor = np.array(
            [[0.0, 0.0, 0.0], [-1, -1, -1], [-1, -1, -1], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
        )

        self.curve = pysurf.tsurf_tools.create_curve_from_points(coor, "test_curve")

        pts = np.array([[3.0, 0.0, 0.0], [2.0, 1.0, 0.5]])
        projPts = pts.copy()

        self.curve.project(pts, xyzProj=projPts)

        np.testing.assert_allclose(projPts, np.array([[2.0, 0.0, 0.0], [1.5, 0.5, 0.0]]))


if __name__ == "__main__":
    unittest.main()
