# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import os

# TESTING FUNCTION

os.system('rm curve*')

cube = pysurf.TSurfComponent('../inputs/wingBody.cgns', ['wing'], MPI.COMM_WORLD)
#cube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', MPI.COMM_WORLD)

featureCurves = pysurf.tsurf_geometry.extract_curves_from_surface(cube)

curveID = 0
for curve in featureCurves:
    curveID = curveID + 1 
    curve.export_plot3d('curve_%03d'%curveID)

#os.system('tec360 layout_cube.lay')
