# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

os.system('rm curve_*')

#comp1.curves['geom_1d'].export_plot3d('plate_outline')
#comp2.curves['cylinder_1d'].export_plot3d('cyl_outline')

distTol = 1e-7

#========================================================

def curve_intersection(comp1,comp2,deltaZ,ii):

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    comp1.translate(0.0,0.0,deltaZ)

    # Call intersection function
    comp1.intersect(comp2)

    # Testing derivatives

    # Get intersection curve
    for curve in comp1.curves:
        if 'int' in curve:

            # Extract intersection curve
            intCurve = comp1.curves[curve]

            # Remesh the curve in a linear spacing
            newIntCurve = intCurve.remesh()

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

            newCoorIntb = np.random.rand(newIntCurve.coor.shape[0],newIntCurve.coor.shape[1])

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

            # Save the curve
            # newIntCurve.export_tecplot('curve_%03d'%ii)

            # Find the trailing edge point
            pt0 = newIntCurve.barsConn[0,0]###
            pt1 = newIntCurve.barsConn[-1,-1]###
            X0 = newIntCurve.coor[0,pt0-1]
            X1 = newIntCurve.coor[0,pt1-1]
            if X0 > X1:
                barID = 0 ###
                pointID = pt0
            else:
                barID = -1 ###
                pointID = pt1

            # Now compute the derivative for the Y coordinate of the first point of the intersection
            X = newIntCurve.coor[0,pointID-1]
            Y = newIntCurve.coor[1,pointID-1]
            Z = newIntCurve.coor[2,pointID-1]
            newCoorIntb[:,:] = 0.0
            newCoorIntb[1,pointID-1] = 1.0

            # Compute the remesh derivatives
            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)
            dYdZ = np.sum(coorAb[2,:])

    # Remove translation
    comp1.translate(0.0,0.0,-deltaZ)

    # Print results
    print 'results'
    print 'Y:',Y
    print 'dYdZ:',dYdZ
    print ''

    # Return results
    return Y, dYdZ

#========================================================

#========================================================

def cylinder_level(level):

    # Load plate as component 1
    comp1 = pysurf.TSurfGeometry('../inputs/plate.cgns',['geom','geom_1d'])

    comp1.scale(1.0/200.0)
    comp1.rotate(0.0,1)
    comp1.translate(1.5, 0.8, -0.3)

    # Load cylinder as component 2
    comp2 = pysurf.TSurfGeometry('../inputs/cylinderAlone_L%d.cgns'%level,['geom'])

    comp2.scale(2.0)

    nStates = 11
    Z = np.linspace(0.0, 0.6, nStates)

    Y = np.zeros(nStates)
    dYdZ = np.zeros(nStates)

    for ii in range(nStates):
        print ''
        print 'translation'
        print Z[ii]
        Y[ii], dYdZ[ii] = curve_intersection(comp1,comp2,Z[ii],ii)
        print 'results'
        print Y[ii]
        print dYdZ[ii]
        print ''

    results = np.vstack([Z,Y,dYdZ])
    with open('results_L%d.pickle'%level,'w') as fid:
        pickle.dump(results,fid)

#========================================================

# TESTING FUNCTION
class TestCylinderDerivs(unittest.TestCase):

    def test_cylinder_3(self):
        cylinder_level(3)

# MAIN PROGRAM
if __name__ == "__main__":
    unittest.main()

    # levels = range(4)
    #
    # for level in levels:
    #     cylinder_level(3-level)
