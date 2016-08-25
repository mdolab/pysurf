# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import copy
import os

# TESTING FUNCTION

# Clean old curve files
os.system('rm curve*')

# Set communicator
comm = MPI.COMM_WORLD

# Wing example
wing = pysurf.TSurfComponent('../inputs/crm.cgns', ['w_upp',#'w_low',
                                                    'te_upp_curve','te_low_curve',
                                                    'w_root_upp_curve','w_tip_upp_curve',
                                                    'w_root_low_curve','w_tip_low_curve',
                                                    'w_root_te_curve','w_tip_te_curve',
                                                    'w_le_curve'])

# Split curves
#pysurf.tsurf_geometry.split_curves(wing.Curves)

'''
# Convert from list to dict
curveDict = wing.Curves

print ' Number of curves found:',len(curveDict)

# Export curves in plot3d format
curveID = 0
for curve in curveDict.itervalues():
    curveID = curveID + 1 
    curve.export_plot3d('curve_%03d'%curveID)
'''

# List relevant curves
tfiCurves = ['te_upp_curve','w_tip_upp_curve','w_le_curve','w_root_upp_curve']

# Initialize list that will contain remeshed curves
curveList = []

# Remesh curves
ii = 0 # Counter to use different refinement
for curveName in tfiCurves:

    print 'Remeshing curve:',curveName

    # Define refinement level
    if np.mod(ii,2) == 0:
        nNewNodes = 400
    else:
        nNewNodes = 200

    # Remesh curve
    wing.Curves[curveName].remesh(nNewNodes=200)

    # Add curve to the list
    curveList.append(wing.Curves[curveName])


# Organize curves
leftCurve, bottomCurve, rightCurve, topCurve = pysurf.tfi_mesh.link_curves(curveList)

# Set TFI problem
mesh = pysurf.tfi_mesh.TFIMesh(leftCurve, bottomCurve, rightCurve, topCurve, wing)

# Solve TFI problem
mesh.generate_mesh()

# Export results
mesh.export('mesh.xyz',addCurves=True)
