# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION
'''
# Big cube and small cube example
comp1 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)
comp2 = pysurf.TSurfGeometry('../inputs/simpleCube.cgns', comm)
comp2.translate(.25, 0., 0.)
comp2.rotate(30, 0)
comp2.rotate(30, 1)
comp2.rotate(30, 2)
comp2.scale(1.5)
'''

'''
# Wing-body example
comp1 = pysurf.TSurfGeometry('../inputs/crm.cgns', ['w_upp','w_low','w_ted'])
comp2 = pysurf.TSurfGeometry('../inputs/crm.cgns', ['b_fwd','b_cnt','b_rrf'])
comp1.translate(0.0, -10.0, 0.0) # Move wing inboard so we have a well-defined intersection
'''

class TestCurveIntersection(unittest.TestCase):

    def test_curve_intersection(self):

        os.system('rm curve*')

        # Set communicator
        comm = MPI.COMM_WORLD

        # Cube and cylinder example
        comp1 = pysurf.TSurfGeometry('../inputs/cube.cgns', comm)
        comp2 = pysurf.TSurfGeometry('../inputs/cylinder.cgns', comm)

        # Call intersection function
        comp1.intersect(comp2)

        # Split intersections
        pysurf.tsurf_tools.split_curves(comp2.curves)

        with open('intersection_dict.pickle', 'r') as f:
            master_dict = pickle.load(f)

        # Export curves in plot3d format
        curveID = 0
        # master_dict = {}
        for curve in comp2.curves.itervalues():
            # Save only intersection curves
            if 'int' in curve.name:
                curveID = curveID + 1
                # curve.export_plot3d('curve_%03d'%curveID)

                # master_dict.update({curveID : curve.coor})
                coor_master = master_dict[curveID]
                np.testing.assert_almost_equal(coor_master, curve.coor)

        # Save the intersection curves for comparison
        # with open('intersection_dict.pickle', 'w') as f:
        #     pickle.dump(master_dict, f)


        # Testing derivatives
        
        # Get intersection curve
        for curve in comp1.curves:
            if 'int' in curve:
                intCurve = comp1.curves[curve]
                coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                                 comp2,
                                                                                 intCurve,
                                                                                 np.ones(intCurve.coor.shape))
                print 'boo'
                print np.min(coorAb)


if __name__ == "__main__":
    unittest.main()
