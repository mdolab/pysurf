# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

os.system('rm *.plt')

# Create one curve
'''
initCurveName = 'test_curve'
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.0],
                 [0.7, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.5, 0.5, 0.0]],order='F').T

barsConn = np.array([[4,5],
                     [6,7],
                     [2,3],
                     [3,4],
                     [1,2],
                     [5,6],
                     [7,8],
                     [8,1]],dtype='int32',order='F').T
'/''
barsConn = np.array([[1,2],
                     [2,3],
                     [3,4],
                     [4,5],
                     [5,6],
                     [6,7],
                     [7,8],
                     [8,1]],dtype='int32',order='F').T
'/''
initCurve = pysurf.TSurfCurve(initCurveName, coor, barsConn)
'''
initCurveDict = pysurf.tsurf_tools.read_tecplot_curves('int_body_wing_000.plt_')
initCurveName = initCurveDict.keys()[0]
initCurve = initCurveDict[initCurveName]

# Create a manager and add the test curve to it
manager = pysurf.Manager()
manager.add_curve(initCurve)

# Export initial curve
manager.intCurves[initCurveName].export_tecplot(initCurveName)

# Now let's remesh this curve
newCurveName = manager.remesh_intCurve(initCurveName)

# Save the remeshed curves
manager.intCurves[newCurveName].export_tecplot(newCurveName)

# DERIVATIVE SEEDS

# Generate random seeds
manager.intCurves[initCurveName].set_randomADSeeds(mode='forward')
manager.intCurves[newCurveName].set_randomADSeeds(mode='reverse')

# Store relevant seeds
initCurveCoord = manager.intCurves[initCurveName].get_forwardADSeeds()
newCurveCoorb = manager.intCurves[newCurveName].get_reverseADSeeds(clean=False)

# FORWARD AD

# Call AD code
manager.forwardAD()

# Get relevant seeds
newCurveCoord = manager.intCurves[newCurveName].get_forwardADSeeds()

# REVERSE AD

# Call AD code
manager.reverseAD()

# Get relevant seeds
initCurveCoorb = manager.intCurves[initCurveName].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(initCurveCoorb*initCurveCoord)
dotProd = dotProd - np.sum(newCurveCoorb*newCurveCoord)
    
print 'dotProd test'
print dotProd

# FINITE DIFFERENCE TEST

# Define step size
stepSize = 1e-9

# Perturb the initial curve
initCurve.coor = initCurve.coor + stepSize*initCurveCoord

# Create new manager
manager2 = pysurf.Manager()
manager2.add_curve(initCurve)

# Now let's split this curve
newCurveName = manager2.remesh_intCurve(initCurveName)

# Compute derivatives with finite differencing
coor0 = manager.intCurves[newCurveName].coor
coor = manager2.intCurves[newCurveName].coor

newCurveCoord_FD = (coor - coor0)/stepSize

# Finite difference test
FD_error = np.max(np.abs(newCurveCoord - newCurveCoord_FD))
print 'FD test'
print FD_error

print np.abs(newCurveCoord - newCurveCoord_FD)[:,5:20]
