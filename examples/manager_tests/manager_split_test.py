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
initCurveName = 'test_curve'
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.0],
                 [0.7, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.5, 0.5, 0.0]],order='F').T
barsConn = np.array([[1,2],
                     [2,3],
                     [3,4],
                     [4,5],
                     [5,6],
                     [6,7],
                     [7,8],
                     [8,1]],dtype='int32',order='F').T
initCurve = pysurf.TSurfCurve(initCurveName, coor, barsConn)

# Create a manager and add the test curve to it
manager = pysurf.Manager()
manager.add_curve(initCurve)

# Export initial curve
manager.intCurves[initCurveName].export_tecplot(initCurveName)

# Now let's split this curve
newCurveNames = manager.split_intCurve(initCurveName)

# Save the split curves
for curve in newCurveNames:
    manager.intCurves[curve].export_tecplot(curve)

# DERIVATIVE SEEDS

# Generate random seeds
manager.intCurves[initCurveName].set_randomADSeeds(mode='forward')
for curve in newCurveNames:
    manager.intCurves[curve].set_randomADSeeds(mode='reverse')

# Store relevant seeds
initCurveCoord = manager.intCurves[initCurveName].get_forwardADSeeds()
splitCurvesCoorb = []
for curve in newCurveNames:
    coorb = manager.intCurves[curve].get_reverseADSeeds(clean=False)
    splitCurvesCoorb.append(coorb)

# FORWARD AD

# Call AD code
manager.forwardAD()

# Get relevant seeds
splitCurvesCoord = []
for curve in newCurveNames:
    coord = manager.intCurves[curve].get_forwardADSeeds()
    splitCurvesCoord.append(coord)

# REVERSE AD

# Call AD code
manager.reverseAD()

# Get relevant seeds
initCurveCoorb = manager.intCurves[initCurveName].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(initCurveCoorb*initCurveCoord)
for ii in range(len(newCurveNames)):
    dotProd = dotProd - np.sum(splitCurvesCoorb[ii]*splitCurvesCoord[ii])
    
print 'dotProd test'
print dotProd

# FINITE DIFFERENCE TEST

# Define step size
stepSize = 1e-7

# Perturb the initial curve
initCurve.set_points(initCurve.get_points() + stepSize*initCurveCoord)

# Create new manager
manager2 = pysurf.Manager()
manager2.add_curve(initCurve)

# Now let's split this curve
newCurveNames = manager2.split_intCurve(initCurveName)

# Compute derivatives with finite differencing
splitCurvesCoord_FD = []
for curve in newCurveNames:
    coor0 = manager.intCurves[curve].get_points()
    coor = manager2.intCurves[curve].get_points()

    curr_FD = (coor - coor0)/stepSize
    splitCurvesCoord_FD.append(curr_FD)

# Finite difference test
FD_error = 0.0
for ii in range(len(newCurveNames)):
    curr_error = np.max(np.abs(splitCurvesCoord[ii] - splitCurvesCoord_FD[ii]))
    FD_error = max(curr_error, FD_error)
print 'FD test'
print FD_error
