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

#========================================================

def test_curve_intersection(deltaZ,ii):

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

            # Running backward mode

            coorIntb = np.random.rand(intCurve.coor.shape[0],intCurve.coor.shape[1])

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb)

            # Running forward mode

            coorAd = np.array(np.random.rand(comp1.coor.shape[0],comp1.coor.shape[1]),order='F')
            coorBd = np.array(np.random.rand(comp2.coor.shape[0],comp2.coor.shape[1]),order='F')

            coorIntd = pysurf.tsurf_tools._compute_pair_intersection_d(comp1,
                                                                       comp2,
                                                                       intCurve,
                                                                       coorAd,
                                                                       coorBd)

            # Dot product test
            dotProd = 0.0
            dotProd = dotProd + np.sum(coorIntb*coorIntd)
            dotProd = dotProd - np.sum(coorAb*coorAd)
            dotProd = dotProd - np.sum(coorBb*coorBd)

            # Print results
            print 'dotProd test'
            print dotProd

            # Save the curve
            intCurve.export_plot3d('curve_%03d'%ii)

            # Find the upper skin trailing edge point
            pt0 = intCurve.barsConn[0,0]
            pt1 = intCurve.barsConn[-1,-1]
            Z0 = intCurve.coor[2,pt0-1]
            Z1 = intCurve.coor[2,pt1-1]
            if Z0 > Z1:
                pointID = pt0
            else:
                pointID = pt1

            Y = intCurve.coor[1,pointID-1]
            coorIntb[:,:] = 0.0
            coorIntb[1,pointID-1] = 1.0
            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb)
            dYdZ = np.sum(coorAb[2,:])

    # Remove translation
    comp1.translate(0,+100.0,-deltaZ)

    # Return results
    return Y, dYdZ

#========================================================

# MIN PROGRAM
nStates = 20
Z = np.linspace(0.0, 140.0, nStates)

Y = np.zeros(len(Z))
dYdZ = np.zeros(len(Z))

for ii in range(len(Z)):
    print ''
    print 'translation'
    print Z[ii]
    Y[ii], dYdZ[ii] = test_curve_intersection(Z[ii],ii)
    print 'results'
    print Y[ii]
    print dYdZ[ii]
    print ''

results = np.vstack([Z,Y,dYdZ])
with open('results.pickle','w') as fid:
    pickle.dump(results,fid)
