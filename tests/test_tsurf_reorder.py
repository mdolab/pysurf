import pysurf
import numpy as np
import unittest


class TestReorder(unittest.TestCase):

    N_PROCS = 2

    def test_reorder(self):

        # Define number of nodes in the circle
        nNodes = 5

        # Initialize coor array
        coor = np.zeros((nNodes, 3))

        # Initialize angles
        theta = np.linspace(0, 2.0 * np.pi, nNodes + 1)

        # Compute coordinates
        for ii in range(nNodes):
            coor[ii, 0] = np.cos(theta[ii])
            coor[ii, 1] = np.sin(theta[ii])

        # Generate connectivity
        nNodes = coor.shape[0]
        barsConn = np.zeros((nNodes, 2))
        barsConn[:, 0] = range(0, nNodes)
        barsConn[:, 1] = range(1, nNodes + 1)

        # Make it periodic
        barsConn[-1, 1] = barsConn[0, 0]

        # Create curve object
        curve = pysurf.TSurfCurve(coor, barsConn, "test")

        # Reorder curve
        curve.shift_end_nodes(criteria="maxY")

        np.testing.assert_allclose(curve.barsConn, np.array([[1, 2], [2, 3], [3, 4], [4, 0], [0, 1]]))


if __name__ == "__main__":
    unittest.main()
