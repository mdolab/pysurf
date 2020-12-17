# IMPORTS
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

def curve_intersection_internal(deltaZ,ii):

    '''
    This version calls the intersections functions directly, bypassing all methods
    of the geometry objects.
    '''

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    comp1.translate(0,-100.0,deltaZ)

    # Call intersection function
    intNames = comp1.intersect(comp2)

    # Testing derivatives

    # Get intersection curve
    for curve in comp1.curves:
        if 'int' in curve:
            intCurve = comp1.curves[curve]
            intCurve2 = intCurve.remesh()
            intCurve.export_tecplot(outputName='init_curve')

            # Save initial coordinates of the curve
            coorInt0 = np.array(intCurve2.coor,order='F')

            np.random.seed(123)

            # Running forward mode
            coor1d = np.array(np.random.rand(comp1.coor.shape[0],comp1.coor.shape[1]),order='F')
            coor1d = coor1d/np.sqrt(np.sum(coor1d**2))
            coor2d = np.array(np.random.rand(comp2.coor.shape[0],comp2.coor.shape[1]),order='F')
            coor2d = coor2d/np.sqrt(np.sum(coor2d**2))

            coorIntd = pysurf.tsurf_tools._compute_pair_intersection_d(comp1,
                                                                       comp2,
                                                                       intCurve,
                                                                       coor1d,
                                                                       coor2d,
                                                                       distTol)

            #newCoorIntd = pysurf.tsurf_tools._remesh_d(intCurve, coorIntd)
            newCoorIntd = coorIntd # Use this to bypass the remesh (but comment the remesh step)

            # Running backward mode
            newCoorIntb = np.array(np.random.rand(intCurve.coor.shape[0],intCurve.coor.shape[1]),order='F')
            newCoorIntb_copy = np.array(newCoorIntb)

            #coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)
            coorIntb = newCoorIntb # Use this to bypass the remesh (but comment the remesh step)

            coor1b, coor2b = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)

            # Dot product test
            dotProd = 0.0
            dotProd = dotProd + np.sum(newCoorIntb_copy*newCoorIntd)
            dotProd = dotProd - np.sum(coor1b*coor1d)
            dotProd = dotProd - np.sum(coor2b*coor2d)

            # Run perturbed geometry to get finite differences
            stepSize = 1e-1
            comp1.update(comp1.coor + stepSize*coor1d)
            comp2.update(comp2.coor + stepSize*coor2d)
            Intersections = pysurf.tsurf_tools._compute_pair_intersection(comp1,
                                                                          comp2,
                                                                          distTol)
            Intersections[0].export_tecplot(outputName='pert_curve')
            Intersections[0] = Intersections[0].remesh(nNewNodes=intCurve.coor.shape[1])

            coorIntPert = np.array(Intersections[0].coor,order='F')
            comp1.update(comp1.coor - stepSize*coor1d)
            comp2.update(comp2.coor - stepSize*coor2d)

            coorIntd_FD = (coorIntPert - coorInt0)/stepSize
            newCoorIntd_FD = coorIntd_FD

            # Finite difference test
            FD_error = np.max(np.abs(newCoorIntd - newCoorIntd_FD))

            print 'FD test (this should be around',stepSize,' or less)'
            print FD_error

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
            Z0 = newIntCurve.coor[pt0, 2]
            Z1 = newIntCurve.coor[pt1, 2]
            if Z0 > Z1:
                pointID = pt0
            else:
                pointID = pt1

            # Now compute the derivative for the Y coordinate of the first point of the intersection
            X = newIntCurve.coor[pointID, 0]
            Y = newIntCurve.coor[pointID, 1]
            Z = newIntCurve.coor[pointID, 2]
            newCoorIntb[:,:] = 0.0
            newCoorIntb[pointID, 1] = 1.0

            # Get nodes of the parent triangles
            parentTriaA = newIntCurve

            # Compute the remesh derivatives
            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)
            dYdZ = np.sum(coorAb[:, 2])

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

def curve_intersection(deltaZ,ii):

    '''
    This version uses the geometry object methods to compute derivatives.
    '''

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    comp1.translate(0,-100.0,deltaZ)

    # Call intersection function
    intCurves = comp1.intersect(comp2)
    intCurve = intCurves[0]

    # Testing derivatives

    # Get intersection curve
    intCurve2 = intCurve.remesh()
    intCurve.export_tecplot(outputName='init_curve')

    # Save initial coordinates of the curve
    coorInt0 = np.array(intCurve2.coor,order='F')

    np.random.seed(123)

    # FORWARD MODE

    coor1d, curveCoor1d = comp1.set_randomADSeeds(mode='forward')
    coor2d, curveCoor2d = comp2.set_randomADSeeds(mode='forward')

    comp1.intersect_d(comp2, intCurve, distTol=distTol)

    # Get derivative seeds defined on both sides (which should be the same in this case)
    newCoorIntd = intCurve.get_forwardADSeeds()

    # REVERSE MODE

    # Store seeds of the intersection curve
    newCoorIntb = intCurve.set_randomADSeeds(mode='reverse')

    comp1.intersect_b(comp2, intCurve, distTol=distTol, accumulateSeeds=False)
    coor1b, curveCoor1b = comp1.get_reverseADSeeds()
    coor2b, curveCoor2b = comp2.get_reverseADSeeds()

    # Dot product test
    dotProd = 0.0
    dotProd = dotProd + np.sum(newCoorIntb*newCoorIntd)
    dotProd = dotProd - np.sum(coor1b*coor1d)
    dotProd = dotProd - np.sum(coor2b*coor2d)

    '''
    # Run perturbed geometry to get finite differences
    stepSize = 1e-1
    comp1.update(comp1.coor + stepSize*coor1d)
    comp2.update(comp2.coor + stepSize*coor2d)
    Intersections = pysurf.tsurf_tools._compute_pair_intersection(comp1,
                                                                  comp2,
                                                                  distTol)
    Intersections[0].export_tecplot(outputName='pert_curve')
    Intersections[0] = Intersections[0].remesh(nNewNodes=intCurve.coor.shape[1])

    coorIntPert = np.array(Intersections[0].coor,order='F')
    comp1.update(comp1.coor - stepSize*coor1d)
    comp2.update(comp2.coor - stepSize*coor2d)

    coorIntd_FD = (coorIntPert - coorInt0)/stepSize
    newCoorIntd_FD = coorIntd_FD

    # Finite difference test
    FD_error = np.max(np.abs(newCoorIntd - newCoorIntd_FD))

    print 'FD test (this should be around',stepSize,' or less)'
    print FD_error
    '''

    # Print results; should be 0
    print 'dotProd test'
    print dotProd
    np.testing.assert_almost_equal(dotProd, 0.)


#========================================================

# TESTING FUNCTION
class TestIntersectionDerivs(unittest.TestCase):

    def test_intersection_derivative(self):
        #curve_intersection_internal(50, 0)
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
