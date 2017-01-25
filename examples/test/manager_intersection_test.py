# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

os.system('rm *.plt')

# Load components
comp1 = pysurf.TSurfGeometry('../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
comp2 = pysurf.TSurfGeometry('../inputs/fuselage_crm4.cgns',['fuse'])
#comp1 = pysurf.TSurfGeometry('../inputs/crm.cgns',['w_upp','w_low'])
#comp2 = pysurf.TSurfGeometry('../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

name1 = 'wing'
name2 = 'body'

comp1.rename(name1)
comp2.rename(name2)

# Create manager object and add the geometry objects to it
manager = pysurf.Manager()
manager.add_geometry(comp1)
manager.add_geometry(comp2)

distTol = 1e-7

#========================================================

def curve_intersection(deltaZ,ii):

    '''
    This version uses the geometry object methods to compute derivatives.
    '''

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    #comp1.translate(0,-100.0,deltaZ)
    manager.geoms['wing'].translate(0,0.0,deltaZ)

    # FORWARD PASS
    
    # Call intersection function
    manager.intersect(distTol=distTol)

    curveName = manager.intCurves.keys()[0]
    manager.intCurves[curveName].export_tecplot(curveName)

    newCurveNames = manager.split_intCurve(curveName)

    for curve in newCurveNames:
        manager.intCurves[curve].export_tecplot(curve)

    mergedCurveName = 'merged_curve'
    manager.merge_intCurves(newCurveNames, mergedCurveName)
    curveName = mergedCurveName

    manager.intCurves[mergedCurveName].export_tecplot('curve')

    curveName = manager.remesh_intCurve(curveName)
    manager.intCurves[curveName].export_tecplot(curveName)

    # DERIVATIVE SEEDS

    # Generate random seeds
    manager.geoms[name1].set_randomADSeeds(mode='forward')
    manager.geoms[name2].set_randomADSeeds(mode='forward')
    manager.intCurves[curveName].set_randomADSeeds(mode='reverse')

    # Store relevant seeds
    coor1d, _ = manager.geoms[name1].get_forwardADSeeds()
    coor2d, _ = manager.geoms[name2].get_forwardADSeeds()
    intCoorb = manager.intCurves[curveName].get_reverseADSeeds(clean=False)

    # FORWARD AD

    # Call AD code
    manager.forwardAD()

    # Get relevant seeds
    intCoord = manager.intCurves[curveName].get_forwardADSeeds()

    # REVERSE AD

    # Call AD code
    manager.reverseAD()

    # Get relevant seeds
    coor1b,_ = manager.geoms[name1].get_reverseADSeeds()
    coor2b,_ = manager.geoms[name2].get_reverseADSeeds()

    # Dot product test
    dotProd = 0.0
    dotProd = dotProd + np.sum(intCoorb*intCoord)
    print dotProd
    dotProd = dotProd - np.sum(coor1b*coor1d)
    print dotProd
    dotProd = dotProd - np.sum(coor2b*coor2d)
    print dotProd

    print 'dotProd test'
    print dotProd
    np.testing.assert_almost_equal(dotProd, 0.)

    '''
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

            # FORWARD MODE

            comp1.set_randomADSeeds()
            comp2.set_randomADSeeds()

            coor1d, _ = comp1.get_forwardADSeeds()
            coor2d, _ = comp2.get_forwardADSeeds()

            comp1.intersect_d(comp2, distTol=distTol)

            # Get derivative seeds defined on both sides (which should be the same in this case)
            newCoorInt1d = comp1.curves[curve].get_forwardADSeeds()
            newCoorInt2d = comp2.curves[curve].get_forwardADSeeds()

            # REVERSE MODE

            # Store seeds of the intersection curve
            newCoorInt1b = comp1.curves[curve].get_reverseADSeeds(clean=False)
            newCoorInt2b = comp2.curves[curve].get_reverseADSeeds(clean=False)
            
            comp1.intersect_b(comp2, distTol=distTol, accumulateSeeds=False)
            coor1b, _ = comp1.get_reverseADSeeds()
            coor2b, _ = comp2.get_reverseADSeeds()

            # Dot product test
            dotProd = 0.0
            dotProd = dotProd + np.sum(newCoorInt1b*newCoorInt1d)
            dotProd = dotProd + np.sum(newCoorInt2b*newCoorInt2d)
            dotProd = dotProd - np.sum(coor1b*coor1d)
            dotProd = dotProd - np.sum(coor2b*coor2d)

            '/''
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
            '/''

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

            # Compute the remesh derivatives
            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb,
                                                                             distTol)
            dYdZ = np.sum(coorAb[2,:])

    '''

    # Remove translation
    manager.geoms['wing'].translate(0.0,+100.0,-deltaZ)
    #comp1.translate(0.0,+100.0,-deltaZ)

    '''
    # Print results
    print 'results'
    print 'Y:',Y
    print 'dYdZ:',dYdZ
    print ''

    # Return results
    return Y, dYdZ
    '''

#========================================================

# TESTING FUNCTION
class TestIntersectionDerivs(unittest.TestCase):

    def test_intersection_derivative(self):
        #curve_intersection_internal(50, 0)
        curve_intersection(0, 0)

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
