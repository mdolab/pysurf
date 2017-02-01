# IMPORTS
from __future__ import division
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

name1 = comp1.name
name2 = comp2.name

#name1 = 'wing'
#name2 = 'body'

#comp1.rename(name1)
#comp2.rename(name2)

# ADDING GUIDE CURVES
# Create curve dictionary based on imported curves
# !!! Make sure to call `extract_curves.py` before running this script
curves = []
curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_000.plt_'))
curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_001.plt_'))
curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_002.plt_'))
curves.append(pysurf.tsurf_tools.read_tecplot_curves('extracted_curve_003.plt_'))

# Create an empty list in which we'll store the long (longitudinal) edges of
# the rect
long_curves = []
counter = 0

# Examine each of the imported curves
for ext_curve in curves:

    # Split these curves based on sharpness to get all edges of the rect
    split_curve = pysurf.tsurf_tools.split_curve_single(ext_curve[ext_curve.keys()[0]], 'int', criteria='sharpness')

    # Loop over these split curves
    for name in split_curve:

        # Give the curves new names so they do not conflict with each other
        split_curve[name].name = 'int_'+'{}'.format(counter).zfill(3)
        counter += 1

        # Extract the curve points and compute the length
        pts = split_curve[name].get_points()
        length = pts[2, 0] - pts[2, -1]

        # Hardcoded logic here based on the length of edges
        if np.abs(length) > 1:

            # Flip the long curve if it's facing the 'wrong' direction
            if length > 0:
                split_curve[name].flip()

            # Add the long curve to the list
            long_curves.append(split_curve[name])

# Create a list of guideCurves based on the extracted long curves
# Note here we need the strings, not the curve objects.
# We use the same loop to add the guide curves to the rectangle object
guideCurves = []
for ext_curve in long_curves:
    print ext_curve.name
    comp2.add_curve(ext_curve)
    if ext_curve.name != 'int_011':
        guideCurves.append(ext_curve.name)

# Rotate the rectangle in 5 degrees
comp2.rotate(5,2)

# Create manager object and add the geometry objects to it
manager0 = pysurf.Manager()
manager0.add_geometry(comp1)
manager0.add_geometry(comp2)

distTol = 1e-7

# Set up integer to export different meshes
mesh_pass = 0

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
    optionsDict = {
        'nNewNodes':41,
        'spacing':'linear',
        'initialSpacing':0.005,
        'finalSpacing':0.005,
    }
    remeshedCurveName = manager.remesh_intCurve(intCurveName,optionsDict)

    manager.intCurves[remeshedCurveName].export_tecplot(remeshedCurveName)

    # MARCH SURFACE MESHES
    meshName = 'mesh'
    
    options_rect = {
    
        'bc1' : 'curve:int_011',
        'bc2' : 'curve:int_011',
        'dStart' : 0.03,
        'numLayers' : 17,
        'extension' : 3.5,
        'epsE0' : 4.5,
        'theta' : -0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'nuArea' : 0.16,
        'numAreaPasses' : 20,
        'sigmaSplay' : 0.3,
        'cMax' : 10000.0,
        'ratioGuess' : 1.5,
        'guideCurves':guideCurves,
        
    }

    options_cube = {
    
        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.02,
        'numLayers' : 17,
        'extension' : 2.5,
        'epsE0' : 4.5,
        'theta' : -0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'nuArea' : 0.16,
        'numAreaPasses' : 20,
        'sigmaSplay' : 0.3,
        'cMax' : 10000.0,
        'ratioGuess' : 1.5,
        
    }

    meshName = 'mesh'
    meshNames = manager.march_intCurve_surfaceMesh(remeshedCurveName, options0=options_cube, options1=options_rect, meshName=meshName)

    # EXPORT
    for meshName in meshNames:
        manager.meshes[meshName].exportPlot3d(meshName+'_'+str(mesh_pass)+'.xyz')

    return meshNames

# END OF forward_pass
#======================================================

# Call the forward pass function to the original manager
meshNames = forward_pass(manager0)
mesh_pass = mesh_pass + 1

# DERIVATIVE SEEDS

# Generate random seeds
manager0.geoms[name1].set_randomADSeeds(mode='forward')
manager0.geoms[name2].set_randomADSeeds(mode='forward')
for mesh in manager0.meshes.itervalues():
    mesh.set_reverseAD_randomSeeds()

# Store relevant seeds
coor1d = manager0.geoms[name1].get_forwardADSeeds()
coor2d = manager0.geoms[name2].get_forwardADSeeds()
meshb = []
for mesh in manager0.meshes.itervalues():
    meshb.append(mesh.get_reverseAD_outputSeeds())

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
meshd = []
for mesh in manager0.meshes.itervalues():
    meshd.append(mesh.get_forwardAD_outputSeeds())

# REVERSE AD

# Call AD code
manager0.reverseAD()
    
# Get relevant seeds
coor1b = manager0.geoms[name1].get_reverseADSeeds()
coor2b = manager0.geoms[name2].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
for ii in range(len(meshd)):
    dotProd = dotProd + np.sum(meshd[ii]*meshb[ii])
dotProd = dotProd - np.sum(coor1b*coor1d)
dotProd = dotProd - np.sum(coor2b*coor2d)

print 'dotProd test (this will be repeated at the end as well)'
print dotProd

# FINITE DIFFERENCE
stepSize = 1e-7

# Apply perturbations to the geometries
comp1.update(comp1.coor + stepSize*coor1d)
comp2.update(comp2.coor + stepSize*coor2d)

# Create new manager with perturbed components
manager1 = pysurf.Manager()
manager1.add_geometry(comp1)
manager1.add_geometry(comp2)

# Do the forward pass to the perturbed case
meshNames = forward_pass(manager1)
mesh_pass = mesh_pass + 1

# Get coordinates of the meshes in both cases to compute FD derivatives
meshd_FD = []
for meshName in manager0.meshes:
    meshCoor0 = manager0.meshes[meshName].mesh[:,:,:]
    meshCoor1 = manager1.meshes[meshName].mesh[:,:,:]

    # Compute derivatives with FD
    meshCoord_FD = (meshCoor1 - meshCoor0)/stepSize

    # Append to the dictionary
    meshd_FD.append(meshCoord_FD)

def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

# Find the largest difference in derivatives
FD_error = 0.0
for ii in range(len(manager0.meshes)):
    print meshNames[ii]
    curr_error = np.max(np.abs(meshd[ii] - meshd_FD[ii]))
    '''
    view_mat(np.abs(meshd[ii][0,:,:] - meshd_FD[ii][0,:,:]))
    view_mat(np.abs(meshd[ii][1,:,:] - meshd_FD[ii][1,:,:]))
    view_mat(np.abs(meshd[ii][2,:,:] - meshd_FD[ii][2,:,:]))
    view_mat(meshd_FD[ii][0,:,:])
    view_mat(meshd_FD[ii][1,:,:])
    view_mat(meshd_FD[ii][2,:,:])
    '''
    FD_error = max(FD_error, curr_error)

# Print results
print 'dotProd test'
print dotProd
print 'FD test'
print FD_error
