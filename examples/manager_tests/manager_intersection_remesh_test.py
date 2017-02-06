# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

os.system('rm *.plt')

# Load components
comp1 = pysurf.TSurfGeometry('../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
comp2 = pysurf.TSurfGeometry('../inputs/fuselage_crm4.cgns',['fuse'])
#comp1 = pysurf.TSurfGeometry('../inputs/crm.cgns',['w_upp','w_low'])
#comp2 = pysurf.TSurfGeometry('../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

#name1 = 'wing'
#name2 = 'body'
#comp1.rename(name1)
#comp2.rename(name2)

name1 = comp1.name
name2 = comp2.name

comp1.translate(0.0, 0.0, 0.0001)

# Create manager object and add the geometry objects to it
manager = pysurf.Manager()
manager.add_geometry(comp1)
manager.add_geometry(comp2)

distTol = 1e-7

# FORWARD PASS
    
# Call intersection function
manager.intersect(distTol=distTol)
curveName = manager.intCurves.keys()[0]

# Call remeshing function
#curveName = manager.remesh_intCurve(curveName)
curveName = manager.remesh_intCurve(curveName,{'nNewNodes':200})

manager.intCurves[curveName].export_tecplot(curveName)

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoor1d = manager.geoms[name1].set_randomADSeeds(mode='forward')
coor2d, curveCoor2d = manager.geoms[name2].set_randomADSeeds(mode='forward')
intCoorb = manager.intCurves[curveName].set_randomADSeeds(mode='reverse')

# FORWARD AD

# Call AD code
manager.forwardAD()

# Get relevant seeds
intCoord = manager.intCurves[curveName].get_forwardADSeeds()

# REVERSE AD

# Call AD code
manager.reverseAD()
    
# Get relevant seeds
coor1b, curveCoor1b = manager.geoms[name1].get_reverseADSeeds()
coor2b, curveCoor2b = manager.geoms[name2].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(intCoorb*intCoord)
dotProd = dotProd - np.sum(coor1b*coor1d)
dotProd = dotProd - np.sum(coor2b*coor2d)
for curveName in curveCoor1d:
    dotProd = dotProd - np.sum(curveCoor1b[curveName]*curveCoor1d[curveName])
for curveName in curveCoor2d:
    dotProd = dotProd - np.sum(curveCoor2b[curveName]*curveCoor2d[curveName])

print 'dotProd test'
print dotProd

# FINITE DIFFERENCE
stepSize = 1e-9

comp1.update(comp1.coor + stepSize*coor1d)
for curveName in curveCoor1d:
    comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize*curveCoor1d[curveName])
comp2.update(comp2.coor + stepSize*coor2d)
for curveName in curveCoor2d:
    comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize*curveCoor2d[curveName])

# Create new manager with perturbed components
manager2 = pysurf.Manager()
manager2.add_geometry(comp1)
manager2.add_geometry(comp2)

# Call intersection function
manager2.intersect(distTol=distTol)
curveName = manager2.intCurves.keys()[0]

# Call remeshing function
#curveName = manager2.remesh_intCurve(curveName,{'nNewNodes':intCoord.shape[1]+1})
curveName = manager2.remesh_intCurve(curveName,{'nNewNodes':200})

# Get coordinates of the intersection
coor0 = manager.intCurves[curveName].get_points()
coor = manager2.intCurves[curveName].get_points()

# Compute derivatives with FD
intCoord_FD = (coor - coor0)/stepSize

# Print results
print 'FD test'
print np.max(np.abs(intCoord-intCoord_FD))
