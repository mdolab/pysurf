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

# Define components


# Cube and cylinder example
comp1 = pysurf.TSurfComponent('../inputs/cube.cgns', comm)
comp2 = pysurf.TSurfComponent('../inputs/cylinder.cgns', comm)


'''
# Big cube and small cube example
comp1 = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)
comp2 = pysurf.TSurfComponent('../inputs/simpleCube.cgns', comm)
comp2.translate(.25, 0., 0.)
comp2.rotate(30, 0)
comp2.rotate(30, 1)
comp2.rotate(30, 2)
comp2.scale(1.5)
'''

'''
# Wing-body example
comp1 = pysurf.TSurfComponent('../inputs/crm.cgns', ['w_upp','w_low','w_ted'])
comp2 = pysurf.TSurfComponent('../inputs/crm.cgns', ['b_fwd','b_cnt','b_rrf'])
comp1.translate(0.0, -10.0, 0.0) # Move wing inboard so we have a well-defined intersection
'''

# Call intersection function
comp1.intersect(comp2)

# Split intersections
pysurf.tsurf_tools.split_curves(comp2.Curves)

print comp2.Curves.keys()

# Export curves in plot3d format
curveID = 0
for curve in comp2.Curves.itervalues():
    # Save only intersection curves
    if 'int' in curve.name:
        curveID = curveID + 1
        curve.export_plot3d('curve_%03d'%curveID)
