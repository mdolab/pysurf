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

'''
# Cube and cylinder example
comp1 = pysurf.TSurfComponent('../inputs/cube.cgns', comm)
comp2 = pysurf.TSurfComponent('../inputs/cylinder.cgns', comm)
'''

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

# Wing-body example
comp1 = pysurf.TSurfComponent('../inputs/crm.cgns', ['w_upp','w_low','w_ted'])
comp2 = pysurf.TSurfComponent('../inputs/crm.cgns', ['b_fwd','b_cnt','b_rrf'])
comp1.translate(0.0, -10.0, 0.0) # Move wing inboard so we have a well-defined intersection

# Call intersection function
Intersections = pysurf.compute_intersections([comp1, comp2])

print Intersections[0].barsConn[:,:6]
id1 = Intersections[0].barsConn[0,0]
id2 = Intersections[0].barsConn[1,0]
id3 = Intersections[0].barsConn[0,1]
id4 = Intersections[0].barsConn[1,1]
print Intersections[0].coor[:,id1]
print Intersections[0].coor[:,id2]
print Intersections[0].coor[:,id3]
print Intersections[0].coor[:,id4]

# Convert from list to dict
intersectionDict = {}
for curveID in range(len(Intersections)):
    curveName = 'curve_'+'%02d'%curveID
    intersectionDict[curveName] = Intersections[curveID]

# Split intersections
#pysurf.tsurf_geometry.split_curves(intersectionDict)

print ' Number of intersections found:\n          ', len(intersectionDict)

# Export curves in plot3d format
curveID = 0
for curve in intersectionDict.itervalues():
    curveID = curveID + 1 
    curve.export_plot3d('curve_%03d'%curveID)
