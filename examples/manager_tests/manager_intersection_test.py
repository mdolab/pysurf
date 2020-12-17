# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle
import unittest

class IntersectionTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(IntersectionTest, self).__init__(*args, **kwargs)

    def test_intersection(self):

        # TESTING FUNCTION

        os.system('rm *.plt')

        # Load components
        comp1 = pysurf.TSurfGeometry('../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
        comp2 = pysurf.TSurfGeometry('../inputs/fuselage_crm4.cgns',['fuse'])
        '''
        #comp1 = pysurf.TSurfGeometry('../inputs/crm.cgns',['w_upp','w_low'])
        #comp2 = pysurf.TSurfGeometry('../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])
        #comp1.translate(0.0,-1.0,0.0)
        '''
        name1 = comp1.name
        name2 = comp2.name

        #comp1.rename(name1)
        #comp2.rename(name2)

        # Create manager object and add the geometry objects to it
        manager = pysurf.Manager()
        manager.add_geometry(comp1)
        manager.add_geometry(comp2)

        distTol = 1e-7

        # FORWARD PASS

        # Call intersection function
        manager.intersect(distTol=distTol)

        curveName = manager.intCurves.keys()[0]
        manager.intCurves[curveName].export_tecplot(curveName)

        # DERIVATIVE SEEDS

        # Generate random seeds
        coor1d, curveCoor1d = manager.geoms[name1].set_randomADSeeds(mode='forward')
        coor2d, curveCoor2d = manager.geoms[name2].set_randomADSeeds(mode='forward')
        intCoorb = manager.intCurves[curveName].set_randomADSeeds(mode='reverse')

        # FORWARD AD

        # Call AD code
        manager.forwardAD()

        # Get relevant seeds
        intCoord = manager.intCurves[curveName].get_forwardADSeeds()

        # REVERSE AD

        # Call AD code
        manager.reverseAD()

        # Get relevant seeds
        coor1b, curveCoor1b = manager.geoms[name1].get_reverseADSeeds()
        coor2b, curveCoor2b = manager.geoms[name2].get_reverseADSeeds()

        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(intCoorb*intCoord)
        dotProd = dotProd - np.sum(coor1b*coor1d)
        dotProd = dotProd - np.sum(coor2b*coor2d)
        for curveName in curveCoor1d:
            dotProd = dotProd - np.sum(curveCoor1b[curveName]*curveCoor1d[curveName])
        for curveName in curveCoor2d:
            dotProd = dotProd - np.sum(curveCoor2b[curveName]*curveCoor2d[curveName])

        print('dotProd test')
        print(dotProd)
        np.testing.assert_almost_equal(dotProd, 0., decimal=11)

        # FINITE DIFFERENCE
        # We can't do a FD test here because the number of points will change

if __name__ == "__main__":
    unittest.main()
