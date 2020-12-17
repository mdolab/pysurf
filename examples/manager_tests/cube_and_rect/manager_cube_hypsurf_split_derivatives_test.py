'''
In this example, we split the 4 edges of the intersection curve to
perform individual remesh steps. This avoids the undefined derivatives
at the rectangle edges.
'''

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle

# TESTING FUNCTION

os.system('rm *.plt')

# Load components
comp1 = pysurf.TSurfGeometry('../../inputs/cube_uns.cgns',['geom'])
comp2 = pysurf.TSurfGeometry('../../inputs/rect_uns.cgns',['geom'])

name1 = comp1.name
name2 = comp2.name

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
        length = pts[0, 2] - pts[-1, 2]

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

    # SPLIT

    # Split the intersection curve
    splitCurveNames = manager.split_intCurve(intCurveName,
                                             criteria='sharpness')

    # REMESH
    remeshedCurveNames = []
    for splitCurveName in splitCurveNames:

        # Remesh each splitted curve individually
        optionsDict = {
            'nNewNodes':21,
            'spacing':'linear',
            'initialSpacing':0.005,
            'finalSpacing':0.005,
        }

        # The remesh function returns a name that we will append to the list
        # of curve names so we can merge them later
        remeshedCurveNames.append(manager.remesh_intCurve(splitCurveName,optionsDict))

    # MERGE
    mergedCurveName = 'intersection'
    manager.merge_intCurves(remeshedCurveNames, mergedCurveName)

    # REORDER
    manager.intCurves[mergedCurveName].shift_end_nodes(criteria='maxX')

    manager.intCurves[mergedCurveName].export_tecplot('intersection_'+str(mesh_pass))

    for curve in manager.geoms[name2].curves.itervalues():
        curve.export_tecplot(curve.name)

    # MARCH SURFACE MESHES
    meshName = 'mesh'
    
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

    meshNames = manager.march_intCurve_surfaceMesh(mergedCurveName, options0=options_cube, options1=options_rect, meshName=meshName)

    return meshNames

# END OF forward_pass
#======================================================

# Call the forward pass function to the original manager
meshNames = forward_pass(manager0)
mesh_pass = mesh_pass + 1

manager0.export_meshes('.')

# DERIVATIVE SEEDS

# Generate random seeds
coor1d, curveCoor1d = manager0.geoms[name1].set_randomADSeeds(mode='forward')
coor2d, curveCoor2d = manager0.geoms[name2].set_randomADSeeds(mode='forward')

meshb = []
for meshGen in manager0.meshGenerators.itervalues():
    meshb.append(meshGen.meshObj.set_randomADSeeds(mode='reverse'))

# FORWARD AD

# Call AD code
manager0.forwardAD()

# Get relevant seeds
meshd = []
for meshGen in manager0.meshGenerators.itervalues():
    meshd.append(meshGen.meshObj.get_forwardADSeeds())

# REVERSE AD

# Call AD code
manager0.reverseAD()
    
# Get relevant seeds
coor1b, curveCoor1b = manager0.geoms[name1].get_reverseADSeeds()
coor2b, curveCoor2b = manager0.geoms[name2].get_reverseADSeeds()

# Dot product test
dotProd = 0.0
for ii in range(len(meshd)):
    dotProd = dotProd + np.sum(meshd[ii]*meshb[ii])
dotProd = dotProd - np.sum(coor1b*coor1d)
print dotProd
dotProd = dotProd - np.sum(coor2b*coor2d)
print dotProd
for curveName in curveCoor1b:
    dotProd = dotProd - np.sum(curveCoor1d[curveName]*curveCoor1b[curveName])
    print dotProd
for curveName in curveCoor2b:
    dotProd = dotProd - np.sum(curveCoor2d[curveName]*curveCoor2b[curveName])
    print dotProd

print 'dotProd test (this will be repeated at the end as well)'
print dotProd

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
meshNames = forward_pass(manager1)
mesh_pass = mesh_pass + 1

# Get coordinates of the meshes in both cases to compute FD derivatives
meshd_FD = []
for meshName in manager0.meshGenerators:
    meshCoor0 = manager0.meshGenerators[meshName].meshObj.get_points()
    meshCoor1 = manager1.meshGenerators[meshName].meshObj.get_points()

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
for ii in range(len(manager0.meshGenerators)):
    numNodes = manager0.meshGenerators[manager0.meshGenerators.keys()[ii]].numNodes
    curr_error = np.abs(meshd[ii] - meshd_FD[ii]).reshape((numNodes,-1,3),order='F')

    view_mat(np.abs(curr_error[:,:,0]))
    view_mat(np.abs(curr_error[:,:,1]))
    view_mat(np.abs(curr_error[:,:,2]))

    FD_error = max(FD_error, np.max(curr_error))

# Print results
print 'dotProd test'
print dotProd
print 'FD test'
print FD_error
