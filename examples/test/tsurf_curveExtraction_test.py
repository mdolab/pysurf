# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import os

# TESTING FUNCTION

os.system('rm curve*')

cube = pysurf.TSurfComponent('../inputs/wingBody.cgns', ['wing'], MPI.COMM_WORLD)
#cube = pysurf.TSurfComponent('../inputs/simpleCube.cgns', MPI.COMM_WORLD)
#cube = pysurf.TSurfComponent('../inputs/crm.cgns', ['w_upp','w_low','w_ted'], MPI.COMM_WORLD)

cube.extract_curves()

curveID = 0
for curve in cube.Curves.itervalues():
    curveID = curveID + 1 
    curve.export_plot3d('curve_%03d'%curveID)

#os.system('tec360 layout_cube.lay')
