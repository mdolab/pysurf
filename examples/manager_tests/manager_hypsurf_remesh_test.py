# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle
import copy

# TESTING FUNCTION

os.system('rm *.plt')

# Load components
geom = pysurf.TSurfGeometry('../inputs/cylinder.cgns')
name1 = geom.name

# Flip BC curve
geom.curves['source'].flip()

# Create manager object and add the geometry object to it
manager0 = pysurf.Manager()
manager0.add_geometry(geom)

distTol = 1e-7

def forward_pass(manager):

    # FORWARD PASS
    
    # Erase previous data
    manager.clear_all()

    # First we will create a new "intersection" curve based on the curves defined in the cylinder geometry
    ref_curve = geom.curves['source']
    #curve = copy.deepcopy(ref_curve)
    curve = pysurf.TSurfCurve(ref_curve.name, ref_curve.coor, ref_curve.barsConn)
    curve.extra_data['parentGeoms'] = [geom.name, geom.name]
    manager.add_curve(curve)
    curveName = curve.name
    
    # Remesh the curve
    remeshedCurveName = manager.remesh_intCurve(curveName,{'nNewNodes':15,
                                                           'spacing':'hypTan',
                                                           'initialSpacing':0.01,
                                                           'finalSpacing':0.01,})

    # Set marching options
    
    options1 = {
        
        'bc1' : 'splay',
        'bc2' : 'curve:bc1',
        'dStart' : 0.01,
        'numLayers' : 30,
        'extension' : 2.8,
        'epsE0' : 5.5,
        'theta' : 0.0,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 3,
        'nuArea' : 0.16,
        'numAreaPasses' : 3,
        'sigmaSplay' : 0.2,
        'cMax' : 1.0,
        'ratioGuess' : 10.0,
        
    }
    
    options2 = {
        
        'bc1' : 'splay',
        'bc2' : 'splay',
        'dStart' : 0.01,
        'numLayers' : 30,
        'extension' : 2.8,
        'epsE0' : 5.5,
        'theta' : 0.0,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 3,
        'nuArea' : 0.16,
        'numAreaPasses' : 3,
        'sigmaSplay' : 0.2,
        'cMax' : 1.0,
        'ratioGuess' : 10.0,
        
    }
    
    # Call meshing routine
    meshName = 'mesh'
    meshNames = manager.march_intCurve_surfaceMesh(remeshedCurveName, options1, options2, meshName)

    for meshName in meshNames:
        manager.meshes[meshName].exportPlot3d(meshName+'.xyz')

    return curveName, meshNames

#===============================================

# RUN FORWARD CODE
curveName, meshNames = forward_pass(manager0)

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoord = manager0.geoms[name1].set_randomADSeeds(mode='forward')

intCoord = manager0.intCurves[curveName].set_randomADSeeds(mode='forward')

meshb = []
for meshName in meshNames:
    meshb.append(manager0.meshes[meshName].set_randomADSeeds(mode='reverse'))

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
meshd = []
for meshName in meshNames:
    meshd.append(manager0.meshes[meshName].get_forwardADSeeds())

# REVERSE AD

# Call AD code
manager0.reverseAD()
    
# Get relevant seeds
coor1b, curveCoorb = manager0.geoms[name1].get_reverseADSeeds()
intCoorb = manager0.intCurves[curveName].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(intCoorb*intCoord)
dotProd = dotProd + np.sum(coor1b*coor1d)
for curveName in curveCoord:
    dotProd = dotProd + np.sum(curveCoord[curveName]*curveCoorb[curveName])
for ii in range(len(meshNames)):
    dotProd = dotProd - np.sum(meshb[ii]*meshd[ii])

print 'dotProd test'
print dotProd

# FINITE DIFFERENCE
stepSize = 1e-7

# Store initial coordinates of the seed curve
seedCurveCoor = geom.curves['source'].get_points()

# Perturb the geometries
geom.update(geom.coor + stepSize*coor1d)
for curveName in curveCoord:
    geom.curves[curveName].set_points(geom.curves[curveName].get_points() + stepSize*curveCoord[curveName])
geom.curves['source'].set_points(seedCurveCoor + stepSize*intCoord)

# Create new manager with perturbed components
manager2 = pysurf.Manager()
manager2.add_geometry(geom)

# Execute code with the perturbed geometry
curveName, meshNames = forward_pass(manager2)

# Get coordinates of the mesh nodes
meshd_FD = []
for meshName in meshNames:
    mesh0 = manager0.meshes[meshName].mesh
    mesh = manager2.meshes[meshName].mesh
    curr_meshd = (mesh - mesh0)/stepSize
    meshd_FD.append(curr_meshd)

for ii in range(len(meshNames)):

    # Print results
    print 'FD test'
    print np.max(np.abs(meshd[ii]-meshd_FD[ii]))

print 'dotProd test'
print dotProd
