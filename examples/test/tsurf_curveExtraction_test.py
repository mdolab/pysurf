# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import hashlib
import unittest
import os
import pickle


class TestCurveExtraction(unittest.TestCase):
    def test_curve_extraction(self):

        # TESTING FUNCTION

        os.system("rm curve*")

        cube = pysurf.TSurfGeometry("../inputs/wingBody.cgns", ["wing"], MPI.COMM_WORLD)
        # cube = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', MPI.COMM_WORLD)
        # cube = pysurf.TSurfGeometry('../inputs/crm.cgns', ['w_upp','w_low','w_ted'], MPI.COMM_WORLD)

        cube.extract_curves()

        with open("extraction_dict.pickle", "r") as f:
            master_dict = pickle.load(f)

        curveID = 0
        # master_dict = {}
        for curve in cube.curves.itervalues():
            curveID = curveID + 1
            # curve.export_plot3d('curve_%03d'%curveID)

            # master_dict.update({curveID : curve.coor})
            coor_master = master_dict[curveID]
            np.testing.assert_almost_equal(coor_master, curve.coor)

        # Save the extracted curves for comparison
        # with open('extraction_dict.pickle', 'w') as f:
        #     pickle.dump(master_dict, f)

        # os.system('tec360 layout_cube.lay')


if __name__ == "__main__":
    unittest.main()
