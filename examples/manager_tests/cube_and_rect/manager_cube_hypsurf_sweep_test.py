# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

'''
This script tracks the "spanwise" position of a surface mesh node on the cube
as we move the rectangle vertically.

It generates results.pickle, which will be read by plotDeriv.py
'''

# WING POSITIONS
deltaZ = np.linspace(-0.1, 0.15, 11)

# TRACKING POINT
# Give i coordinate that we will use to create the tracking slice
i_node = 26
j_node = 15

# TESTING FUNCTION

os.system('rm *.plt')
os.system('rm *.xyz')

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
manager = pysurf.Manager()
manager.add_geometry(comp1)
manager.add_geometry(comp2)

distTol = 1e-7

# Set up integer to export different meshes
mesh_pass = 0

#======================================================
# FORWARD PASS

def compute_position(rect_deltaZ, i_track, j_track):

    '''
    This function will apply all geometry operations to compute the
    position af an arbitrary point in the wing surface mesh and its derivative.

    i_track and j_track are the indices of the node on the wing mesh that will be tracked
    to compute position and derivatives.
    '''

    # Clear the manager object first
    manager.clear_all()

    # Translate the wing
    manager.geoms[name2].translate(0.0, rect_deltaZ, 0.0)

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
        manager.meshes[meshName].exportPlot3d(meshName+'_%03d'%mesh_pass+'.xyz')

    # DERIVATIVE SEEDS

    # Get spanwise position of the tracked node
    Y = manager.meshes[meshNames[0]].mesh[2, i_node, j_node]

    # Set up derivative seeds so we can compute the derivative of the spanwise position
    # of the tracked node of the intersection
    meshb = np.zeros(manager.meshes[meshNames[0]].mesh.shape)
    meshb[2, i_node, j_node] = 1.0

    # Set derivatives to the intersection curve
    manager.meshes[meshNames[0]].set_reverseADSeeds(meshb)

    # REVERSE AD
    
    # Call AD code
    manager.reverseAD()
    
    # Get relevant seeds
    coor1b, curveCoor1b = manager.geoms[name1].get_reverseADSeeds()
    coor2b, curveCoor2b = manager.geoms[name2].get_reverseADSeeds()
    
    # Condense derivatives to take into acount the translation
    dYdZ = np.sum(coor2b[1,:])
    for curveName in curveCoor2b:
        dYdZ = dYdZ + np.sum(curveCoor2b[curveName][1,:])

    # Translate the rectangle back
    manager.geoms[name2].translate(0.0, -rect_deltaZ, 0.0)

    # Return the spanwise position of the highest intersection node and its derivative
    return Y, dYdZ

# END OF compute_position
#======================================================


# Now we will execute the function for every wing position

# Initialize arrays to hold results
numPos = len(deltaZ)
Y = np.zeros(numPos)
dYdZ = np.zeros(numPos)

# Compute intersection for every wing position
for ii in range(numPos):
    print ''
    print 'translation'
    print deltaZ[ii]
    Y[ii], dYdZ[ii] = compute_position(deltaZ[ii], i_node, j_node)
    mesh_pass = mesh_pass + 1
    print 'results'
    print Y[ii]
    print dYdZ[ii]
    print ''

results = np.vstack([deltaZ,Y,dYdZ])
with open('results.pickle','w') as fid:
    pickle.dump(results,fid)

