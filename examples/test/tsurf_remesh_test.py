# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import pysurf
import unittest

np.random.seed(1222223)

class TestRemesh(unittest.TestCase):

    def test_remesh(self):
        coor = np.array([[0.0, 0.0, 0.9],
                         [0.1, 0.0, 0.7],
                         [0.5, 0.0, 0.5],
                         [0.8, 0.0, 0.4],
                         [1.0, 0.0, 0.2]],order='F').T

        # Generate connectivity
        nNodes = coor.shape[1]
        barsConn = np.zeros((2,nNodes-1))
        barsConn[0,:] = range(1,nNodes)
        barsConn[1,:] = range(2,nNodes+1)

        # Create curve object
        curve = pysurf.TSurfCurve(coor, barsConn, 'test')

        ###### copying from fortran
        # coord = np.random.rand(curve.coor.shape[0], curve.coor.shape[1])

        coord = coor * coor
        coord_copy = coord.copy()

        newCoord = curve._remesh_d(coord)

        # newCoorb = np.random.rand(curve.coor.shape[0], curve.coor.shape[1])

        newCoorb = coor * .5 + 1.
        newCoorb_copy = newCoorb.copy()
        coorb = curve._remesh_b(newCoorb)

        lhs = np.sum(coord_copy * coorb)
        rhs = np.sum(newCoorb_copy * newCoord)

        print
        print 'THIS SHOULD BE ZERO:'
        print lhs-rhs
        curve.remesh(spacing='linear')


        master_curve = np.array([[ 0., 0.1839562, 0.4679125, 0.7671471, 1.],
                                 [0., 0., 0., 0., 0.],
                                 [0.9, 0.6580219, 0.5160438, 0.410951, 0.2]])
        np.testing.assert_almost_equal(master_curve, curve.coor)

if __name__ == "__main__":
    unittest.main()
