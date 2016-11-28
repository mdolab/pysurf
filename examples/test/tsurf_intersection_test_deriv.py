# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

os.system('rm curve_*')

# Load components
wing = pysurf.TSurfGeometry('../inputs/crm.cgns',['w_upp','w_low'])
body = pysurf.TSurfGeometry('../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

comp1 = wing
comp2 = body

distTol = 1e-7

#========================================================

def curve_intersection(deltaZ,ii):

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    wing.translate(0,-100.0,deltaZ)

    # Call intersection function
    wing.intersect(body)

    # Testing derivatives

    # Get intersection curve
    for curve in comp1.curves:
        if 'int' in curve:
            intCurve = comp1.curves[curve]

            # Running forward mode

            coorAd = np.array(np.random.rand(comp1.coor.shape[0],comp1.coor.shape[1]),order='F')
            coorBd = np.array(np.random.rand(comp2.coor.shape[0],comp2.coor.shape[1]),order='F')

            coorIntd = pysurf.tsurf_tools._compute_pair_intersection_d(comp1,
                                                                       comp2,
                                                                       intCurve,
                                                                       coorAd,
                                                                       coorBd,
                                                                       distTol)

            newCoorIntd = pysurf.tsurf_tools._remesh_d(intCurve, coorIntd)


            # Running backward mode

            newCoorIntb = np.random.rand(intCurve.coor.shape[0],intCurve.coor.shape[1])

            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)

            # Dot product test
            dotProd = 0.0
            dotProd = dotProd + np.sum(newCoorIntb*newCoorIntd)
            dotProd = dotProd - np.sum(coorAb*coorAd)
            dotProd = dotProd - np.sum(coorBb*coorBd)

            # Print results; should be 0
            print 'dotProd test'
            print dotProd
            np.testing.assert_almost_equal(dotProd, 0.)

            # Remesh the curve in a linear spacing
            newIntCurve = intCurve.remesh()

            # Save the curve
            newIntCurve.export_tecplot('curve_%03d'%ii)

            # Find the upper skin trailing edge point
            pt0 = newIntCurve.barsConn[0,0]
            pt1 = newIntCurve.barsConn[-1,-1]
            Z0 = newIntCurve.coor[2,pt0-1]
            Z1 = newIntCurve.coor[2,pt1-1]
            if Z0 > Z1:
                pointID = pt0
            else:
                pointID = pt1

            # Now compute the derivative for the Y coordinate of the first point of the intersection
            X = newIntCurve.coor[0,pointID-1]
            Y = newIntCurve.coor[1,pointID-1]
            Z = newIntCurve.coor[2,pointID-1]
            newCoorIntb[:,:] = 0.0
            newCoorIntb[1,pointID-1] = 1.0

            # Get nodes of the parent triangles
            parentTriaA = newIntCurve

            # Compute the remesh derivatives
            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)
            dYdZ = np.sum(coorAb[2,:])

    # Remove translation
    comp1.translate(0.0,+100.0,-deltaZ)

    # Print results
    print 'results'
    print 'Y:',Y
    print 'dYdZ:',dYdZ
    print ''

    # Return results
    return Y, dYdZ

#========================================================

# TESTING FUNCTION
class TestIntersectionDerivs(unittest.TestCase):

    def test_intersection_derivative(self):
        curve_intersection(50, 0)

# MAIN PROGRAM
if __name__ == "__main__":
    unittest.main()

    # # MAIN PROGRAM
    # nStates = 15
    # Z = np.linspace(0.0, 140.0, nStates)
    #
    # Y = np.zeros(len(Z))
    # dYdZ = np.zeros(len(Z))
    #
    # for ii in range(len(Z)):
    #     print ''
    #     print 'translation'
    #     print Z[ii]
    #     Y[ii], dYdZ[ii] = curve_intersection(Z[ii],ii)
    #     print 'results'
    #     print Y[ii]
    #     print dYdZ[ii]
    #     print ''
    #
    # results = np.vstack([Z,Y,dYdZ])
    # with open('results.pickle','w') as fid:
    #     pickle.dump(results,fid)
