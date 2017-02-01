# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

os.system('rm *.plt')

# Create a list of curves to be merged
initCurveList = []
initCurveNames = []

initCurveName = 'test_curve1'
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
                     [3,4]],dtype='int32',order='F').T
initCurve = pysurf.TSurfCurve(initCurveName, coor, barsConn)
initCurveList.append(initCurve)
initCurveNames.append(initCurveName)

initCurveName = 'test_curve2'
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.0],
                 [0.7, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.5, 0.5, 0.0]],order='F').T
barsConn = np.array([[4,5],
                     [5,6],
                     [6,7]],dtype='int32',order='F').T
initCurve = pysurf.TSurfCurve(initCurveName, coor, barsConn)
initCurveList.append(initCurve)
initCurveNames.append(initCurveName)

initCurveName = 'test_curve3'
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.0],
                 [0.7, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.5, 0.5, 0.0]],order='F').T
barsConn = np.array([[7,8],
                     [8,1]],dtype='int32',order='F').T
initCurve = pysurf.TSurfCurve(initCurveName, coor, barsConn)
initCurveList.append(initCurve)
initCurveNames.append(initCurveName)

# Create a manager and add the test curves to it
manager = pysurf.Manager()
for curve in initCurveList:
    manager.add_curve(curve)

# Export initial curves
for curve in manager.intCurves:
    manager.intCurves[curve].export_tecplot(curve)

# Now let's merge these curves
mergedCurveName = 'merged_curve'
manager.merge_intCurves(initCurveNames,mergedCurveName)

# Save the merged curve
manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName)

# DERIVATIVE SEEDS

# Generate random seeds
for curve in initCurveNames:
    manager.intCurves[curve].set_randomADSeeds(mode='forward')
manager.intCurves[mergedCurveName].set_randomADSeeds(mode='reverse')

# Store relevant seeds
initCurveCoord = []
for curve in initCurveNames:
    coord = manager.intCurves[curve].get_forwardADSeeds()
    initCurveCoord.append(coord)
mergedCurveCoorb = manager.intCurves[mergedCurveName].get_reverseADSeeds(clean=False)


# FORWARD AD

# Call AD code
manager.forwardAD()

# Get relevant seeds
mergedCurveCoord = manager.intCurves[mergedCurveName].get_forwardADSeeds()

# REVERSE AD

# Call AD code
manager.reverseAD()

# Get relevant seeds
initCurveCoorb = []
for curve in initCurveNames:
    coorb = manager.intCurves[curve].get_reverseADSeeds(clean=False)
    initCurveCoorb.append(coorb)

# Dot product test
dotProd = 0.0
for ii in range(len(initCurveNames)):
    dotProd = dotProd + np.sum(initCurveCoorb[ii]*initCurveCoord[ii])
dotProd = dotProd - np.sum(mergedCurveCoorb*mergedCurveCoord)
   
print 'dotProd test'
print dotProd

# FINITE DIFFERENCE TEST
# Define step size
stepSize = 1e-7

# Perturb the initial curves
for ii in range(len(initCurveList)):
    initCurveList[ii].set_points(initCurveList[ii].get_points() + stepSize*initCurveCoord[ii])

# Create new manager
manager2 = pysurf.Manager()
for curve in initCurveList:
    manager2.add_curve(curve)

# Now let's merge these curves
mergedCurveName = 'merged_curve'
manager2.merge_intCurves(initCurveNames,mergedCurveName)

# Get coordinates of the merged curve on the original and perturbed case
mergedCurveCoor0 = manager.intCurves[mergedCurveName].get_points()
mergedCurveCoor = manager2.intCurves[mergedCurveName].get_points()

# Compute derivatives with finite differencing
mergedCurveCoord_FD = (mergedCurveCoor - mergedCurveCoor0)/stepSize

# Finite difference test
FD_error = np.max(np.abs(mergedCurveCoord - mergedCurveCoord_FD))
print 'FD test'
print FD_error

'''
# Define step size
stepSize = 1e-7

# Perturb the initial curves
for ii in range(len(initCurveList)):
    initCurveList[ii].coor = initCurveList[ii].coor + stepSize*initCurveCoord[ii]

# Create new manager
manager2 = pysurf.Manager()
for curve in initCurveList:
    manager2.add_curve(curve)

# Now let's merge these curves
mergedCurveName = 'merged_curve'
manager2.merge_intCurves(initCurveNames,mergedCurveName)

# Get coordinates of the merged curve on the original and perturbed case
mergedCurveCoor0 = manager.intCurves[mergedCurveName].coor
mergedCurveCoor = manager2.intCurves[mergedCurveName].coor

# Compute derivatives with finite differencing
mergedCurveCoord_FD = (mergedCurveCoor - mergedCurveCoor0)/stepSize

# Finite difference test
FD_error = np.max(np.abs(mergedCurveCoord - mergedCurveCoord_FD))
print 'FD test'
print FD_error
'''
