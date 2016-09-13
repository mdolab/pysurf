# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import pysurf
import unittest

class TestRemesh(unittest.TestCase):

    def test_remesh(self):
        coor = np.array([[0.0, 0.0, 0.9],
                         [0.1, 0.0, 0.7],
                         [0.5, 0.0, 0.5],
                         [0.9, 0.0, 0.4],
                         [1.0, 0.0, 0.2]],order='F').T

        # Generate connectivity
        nNodes = coor.shape[1]
        barsConn = np.zeros((2,nNodes-1))
        barsConn[0,:] = range(1,nNodes)
        barsConn[1,:] = range(2,nNodes+1)

        # Create curve object
        curve = pysurf.TSurfCurve(coor, barsConn, 'test')

        # Remesh curve
        curve.remesh(spacing='cosine')
        master_curve = np.array([[ 0., 0.0855821, 0.4843909, 0.9144179, 1.],
                                 [ 0., 0., 0., 0., 0.],
                                 [ 0.9, .72883587, .50780456, .37116413, .2]])
        np.testing.assert_almost_equal(master_curve, curve.coor)

if __name__ == "__main__":
    unittest.main()
