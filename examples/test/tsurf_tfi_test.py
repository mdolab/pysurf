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
wing = pysurf.TSurfGeometry('../inputs/initial_full_wing_crm4.cgns')

# Load TE curves and append them to the wing component
curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves('../smoothness_check/curve_te_upp.plt')
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('../smoothness_check/curve_te_low.plt')
curveName = curve_te_upp.keys()[0]
wing.add_curve(curveName, curve_te_upp[curveName])
curveName = curve_te_low.keys()[0]
wing.add_curve(curveName, curve_te_low[curveName])

print wing.curves.keys()

# Remove tip section from the trailing edge curves
pysurf.tsurf_tools.split_curve_with_curve(wing.curves, 'curve_te_upp', wing.curves['curve_tip_upp'])
wing.remove_curve('curve_te_upp_01')
wing.rename_curve('curve_te_upp_00','curve_te_upp')
pysurf.tsurf_tools.split_curve_with_curve(wing.curves, 'curve_te_low', wing.curves['curve_tip_low'])
wing.remove_curve('curve_te_low_00')
wing.rename_curve('curve_te_low_01','curve_te_low')


#wing.curves['w_le_curve'].condense_disconnect_curves()
#wing.curves['w_root_upp_curve'].condense_disconnect_curves()
#wing.curves['w_tip_upp_curve'].condense_disconnect_curves(guessNode=5)
#wing.curves['w_tip_upp_curve'].export_tecplot('tipcurve')

# Split curves
#pysurf.tsurf_geometry.split_curves(wing.curves)

'''
# Convert from list to dict
curveDict = wing.curves

print ' Number of curves found:',len(curveDict)

# Export curves in plot3d format
curveID = 0
for curve in curveDict.itervalues():
    curveID = curveID + 1
    curve.export_plot3d('curve_%03d'%curveID)
'''

# List relevant curves
tficurves = ['curve_te_upp','curve_tip_upp','curve_le','curve_root_upp']

# Initialize list that will contain remeshed curves
curveList = []

# Remesh curves
ii = 0 # Counter to use different refinement
for curveName in tficurves:

    print 'Remeshing curve:',curveName

    # Define refinement level
    if np.mod(ii,2) == 0:
        nNewNodes = 400
    else:
        nNewNodes = 200

    # Remesh curve
    wing.curves[curveName] = wing.curves[curveName].remesh(nNewNodes=nNewNodes)

    print wing.curves[curveName].coor.shape

    # Add curve to the list
    curveList.append(wing.curves[curveName])


# Organize curves
leftCurve, bottomCurve, rightCurve, topCurve = pysurf.tfi_mesh.link_curves(curveList)

# Set TFI problem
mesh = pysurf.tfi_mesh.TFIMesh(leftCurve, bottomCurve, rightCurve, topCurve, wing)

# Solve TFI problem
mesh.generate_mesh()

# Export results
mesh.export('mesh.xyz',addcurves=True)
