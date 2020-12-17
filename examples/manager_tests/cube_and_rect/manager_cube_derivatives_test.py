# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# MESH PARAMETERS
numSkinNodes = 129
LE_spacing = 0.01
TE_spacing = 0.001

# TESTING FUNCTION

os.system('rm *.plt')

# Load components
comp1 = pysurf.TSurfGeometry('../../inputs/cube_uns.cgns',['geom'])
comp2 = pysurf.TSurfGeometry('../../inputs/rect_uns.cgns',['geom'])
comp2.rotate(5,2)

name1 = comp1.name
name2 = comp2.name

#name1 = 'wing'
#name2 = 'body'

#comp1.rename(name1)
#comp2.rename(name2)

# Create manager object and add the geometry objects to it
manager0 = pysurf.Manager()
manager0.add_geometry(comp1)
manager0.add_geometry(comp2)

distTol = 1e-7

#======================================================
# FORWARD PASS

def forward_pass(manager):

    '''
    This function will apply all geometry operations to the given manager.
    '''

    # INTERSECT

    # Call intersection function
    intCurveNames = manager.intersect(distTol=distTol)
    intCurveName = intCurveNames[0]

    manager.intCurves[intCurveName].shift_end_nodes(criteria='maxX')

    # REMESH
    '''
    optionsDict = {
        'nNewNodes':200,
        'spacing':'hypTan',
        'initialSpacing':0.001,
        'finalSpacing':0.001,
    }
    '''
    optionsDict = {
        'nNewNodes':41,
        'spacing':'linear',
        'initialSpacing':0.005,
        'finalSpacing':0.005,
    }

    remeshedCurveName = manager.remesh_intCurve(intCurveName,optionsDict)

    manager.intCurves[remeshedCurveName].export_tecplot(remeshedCurveName)

    return remeshedCurveName

# END OF forward_pass
#======================================================

# Call the forward pass function to the original manager
remeshedCurveName = forward_pass(manager0)

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoor1d = manager0.geoms[name1].set_randomADSeeds(mode='forward')
coor2d, curveCoor2d = manager0.geoms[name2].set_randomADSeeds(mode='forward')
intCoorb = manager0.intCurves[remeshedCurveName].set_randomADSeeds(mode='reverse')

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
intCoord = manager0.intCurves[remeshedCurveName].get_forwardADSeeds()

# REVERSE AD

# Call AD code
manager0.reverseAD()

# Get relevant seeds
coor1b, curveCoor1b = manager0.geoms[name1].get_reverseADSeeds()
coor2b, curveCoor2b = manager0.geoms[name2].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(intCoorb*intCoord)
dotProd = dotProd - np.sum(coor1b*coor1d)
dotProd = dotProd - np.sum(coor2b*coor2d)
for curveName in curveCoor1b:
    dotProd = dotProd - np.sum(curveCoor1d[curveName]*curveCoor1b[curveName])
for curveName in curveCoor2b:
    dotProd = dotProd - np.sum(curveCoor2d[curveName]*curveCoor2b[curveName])

print('dotProd test (this will be repeated at the end as well)')
print(dotProd)

# FINITE DIFFERENCE
stepSize = 1e-7

# Apply perturbations to the geometries
comp1.update(comp1.coor + stepSize*coor1d)
for curveName in curveCoor1d:
    comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize*curveCoor1d[curveName])
comp2.update(comp2.coor + stepSize*coor2d)
for curveName in curveCoor2d:
    comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize*curveCoor2d[curveName])

# Create new manager with perturbed components
manager1 = pysurf.Manager()
manager1.add_geometry(comp1)
manager1.add_geometry(comp2)

# Do the forward pass to the perturbed case
remeshedCurveName = forward_pass(manager1)

# Get coordinates of the intersection in both cases
coor0 = manager0.intCurves[remeshedCurveName].get_points()
coor1 = manager1.intCurves[remeshedCurveName].get_points()

# Compute derivatives with FD
intCoord_FD = (coor1 - coor0)/stepSize

# Print results
print('dotProd test')
print(dotProd)
print('FD test')
print(np.max(np.abs(intCoord-intCoord_FD)))
#print np.abs(intCoord-intCoord_FD)
