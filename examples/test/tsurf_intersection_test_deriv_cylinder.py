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
plate = pysurf.TSurfGeometry('../inputs/plate.cgns',['geom','geom_1d'])
cylinder = pysurf.TSurfGeometry('../inputs/cylinder.cgns',['geom'])
#cylinder = pysurf.TSurfGeometry('../inputs/cubeAndCylinder.cgns',['cylinder','cylinder_1d'])

plate.scale(1.0/200.0)
plate.rotate(5.0,1)
plate.translate(1.5, 0.8, -0.3)
'''
cylinder.translate(-0.5, -0.25, -0.5)
cylinder.scale(2.0)
cylinder.rotate(-90,2)
'''
comp1 = plate
comp2 = cylinder
'''
comp1.curves['geom_1d'].export_plot3d('plate_outline')
comp2.curves['cylinder_1d'].export_plot3d('cyl_outline')
'''

#========================================================

def test_curve_intersection(deltaZ,ii):

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
            intCurve = comp1.curves[curve]

            # Running forward mode

            coorAd = np.array(np.random.rand(comp1.coor.shape[0],comp1.coor.shape[1]),order='F')
            coorBd = np.array(np.random.rand(comp2.coor.shape[0],comp2.coor.shape[1]),order='F')


            coorIntd = pysurf.tsurf_tools._compute_pair_intersection_d(comp1,
                                                                       comp2,
                                                                       intCurve,
                                                                       coorAd,
                                                                       coorBd)

            newCoorIntd = pysurf.tsurf_tools._remesh_d(intCurve, coorIntd)


            # Running backward mode

            newCoorIntb = np.random.rand(intCurve.coor.shape[0],intCurve.coor.shape[1])

            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb)

            # Dot product test
            dotProd = 0.0
            dotProd = dotProd + np.sum(newCoorIntb*newCoorIntd)
            dotProd = dotProd - np.sum(coorAb*coorAd)
            dotProd = dotProd - np.sum(coorBb*coorBd)

            # Print results; should be 0
            print 'dotProd test'
            print dotProd

            # Remesh the curve in a linear spacing
            newIntCurve = intCurve.remesh()

            # Save the curve
            newIntCurve.export_plot3d('curve_%03d'%ii)

            '''
            # Check which side of the curve we should pick
            if intCurve.coor[2,intCurve.barsConn[0,0]-1] < intCurve.coor[2,intCurve.barsConn[-1,-1]-1]:
                pointID = intCurve.barsConn[0,0]-1
            else:
                pointID = intCurve.barsConn[-1,-1]-1
            '''
            # Find the trailing edge
            pointID = np.argmax(newIntCurve.coor[0,:])

            # Now compute the derivative for the Y coordinate of the first point of the intersection
            Y = newIntCurve.coor[1,pointID]
            newCoorIntb[:,:] = 0.0
            newCoorIntb[1,pointID] = 1.0

            # Compute the remesh derivatives
            coorIntb = pysurf.tsurf_tools._remesh_b(intCurve, newCoorIntb)

            coorAb, coorBb = pysurf.tsurf_tools._compute_pair_intersection_b(comp1,
                                                                             comp2,
                                                                             intCurve,
                                                                             coorIntb)
            dYdZ = np.sum(coorAb[2,:])

    # Remove translation
    comp1.translate(0.0,0.0,-deltaZ)

    # Return results
    return Y, dYdZ

#========================================================

# MIN PROGRAM
nStates = 21
Z = np.linspace(0.0, 0.6, nStates)

Y = np.zeros(nStates)
dYdZ = np.zeros(nStates)

for ii in range(nStates):
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
