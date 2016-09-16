# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

# Load components
wing = pysurf.TSurfGeometry('../inputs/crm.cgns',['w_upp','w_low','te_low_curve','te_upp_curve'])
body = pysurf.TSurfGeometry('../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

comp1 = wing
comp2 = body

#========================================================

def test_curve_intersection(deltaZ):

    # Set communicator
    comm = MPI.COMM_WORLD

    # Apply translation
    wing.translate(0,-10.0,deltaZ)

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


            # Remesh here
            # wing._remesh_b(np.zeros(wing.coor.shape)+1.)


            # Now compute the derivative for the Y coordinate of the first point of the intersection
            Y = intCurve.coor[1,0]
            coorIntb[:,:] = 0.0
            coorIntb[1,0] = 1.0
            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb)
            dYdZ = np.sum(coorAb[2,:])

    # Remove translation
    comp1.translate(0,+10.0,-deltaZ)

    # Return results
    return Y, dYdZ

#========================================================

Y, dYdZ = test_curve_intersection(0.0)
print Y
print dYdZ
