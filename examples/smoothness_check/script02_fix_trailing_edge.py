# IMPORTS
import pysurf

# Set anme of file that contains the trailing edge
fileName = 'extracted_curve_000.plt'

# Load trailing edge
curves = pysurf.tsurf_tools.read_tecplot_curves(fileName)

# Split upper and lower trailing edges
pysurf.tsurf_tools.split_curves(curves)

# Get curve objects
curveNames = curves.keys()
curve1 = curves[curveNames[0]]
curve2 = curves[curveNames[1]]

# Assing names based on their position
if curve1.coor[0,2] > curve2.coor[0,2]:
    curve1.name = 'curve_te_upp'
    curve2.name = 'curve_te_low'
else:
    curve2.name = 'curve_te_upp'
    curve1.name = 'curve_te_low'

# Save curves
curve1.export_tecplot(curve1.name)
curve2.export_tecplot(curve2.name)
