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

    # First we will trick the manager by assigning a fake intersection curve
    curve = geom.curves['source']
    curve.extra_data['parentGeoms'] = [geom.name, geom.name]
    manager.add_curve(curve)
    curveName = curve.name
    
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
        'remesh':True
        
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
    meshNames = manager.march_intCurve_surfaceMesh(curveName, options1, options2, meshName)

    for meshName in meshNames:
        manager.meshes[meshName].exportPlot3d(meshName+'.xyz')

    return curveName, meshNames

#===============================================

# RUN FORWARD CODE
curveName, meshNames = forward_pass(manager0)

# DERIVATIVE SEEDS

# Generate random seeds
manager0.geoms[name1].set_randomADSeeds(mode='forward')
manager0.intCurves[curveName].set_randomADSeeds(mode='forward')
for meshName in meshNames:
    manager0.meshes[meshName].set_randomADSeeds(mode='reverse')

# Store relevant seeds
coor1d = manager0.geoms[name1].get_forwardADSeeds()
intCoord = manager0.intCurves[curveName].get_forwardADSeeds()
meshb = []
for meshName in meshNames:
    meshb.append(manager0.meshes[meshName].get_reverseADSeeds(clean=False))

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
coor1b = manager0.geoms[name1].get_reverseADSeeds()
intCoorb = manager0.intCurves[curveName].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(intCoorb*intCoord)
dotProd = dotProd + np.sum(coor1b*coor1d)
for ii in range(len(meshNames)):
    dotProd = dotProd - np.sum(meshb[ii]*meshd[ii])

print 'dotProd test'
print dotProd

# FINITE DIFFERENCE
stepSize = 1e-7

geom.update(geom.coor + stepSize*coor1d)
geom.curves['source'].set_points(geom.curves['source'].get_points() + stepSize*intCoord)

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
    print 'FD test for mesh',meshNames[ii]
    print np.max(np.abs(meshd[ii]-meshd_FD[ii]))

print 'dotProd test'
print dotProd
