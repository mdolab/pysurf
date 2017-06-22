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
LE_spacing = 0.001
TE_spacing = 0.01

# WING POSITIONS
deltaZ = np.linspace(0.0, 3.0, 2)

# TRACKING POINT
# Give i coordinate that we will use to create the tracking slice
i_node = 190
j_node = 30

# TESTING FUNCTION

os.system('rm *.plt')
os.system('rm *.xyz')

# Load components
comp1 = pysurf.TSurfGeometry('../../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
comp2 = pysurf.TSurfGeometry('../../inputs/fuselage_crm4.cgns',['fuse'])

#name1 = 'wing'
#name2 = 'body'

#comp1.rename(name1)
#comp2.rename(name2)

name1 = comp1.name
name2 = comp2.name

# Load TE curves and append them to the wing component
curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves('curve_te_upp.plt_')
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('curve_te_low.plt_')
curveName = curve_te_upp.keys()[0]
comp1.add_curve(curve_te_upp[curveName])
curveName = curve_te_low.keys()[0]
comp1.add_curve(curve_te_low[curveName])
comp1.curves[curveName].flip()

# Create manager object and add the geometry objects to it
manager = pysurf.Manager()
manager.add_geometry(comp1)
manager.add_geometry(comp2)

distTol = 1e-7

#======================================================

def compute_position(wing_deltaZ):

    '''
    This function will apply all geometry operations to compute the
    intersection and its derivative at the new location.
    '''

    # Clear the manager object first
    manager.clean_all()

    # Translate the wing
    manager.geoms[name1].translate(0.0, 0.0, wing_deltaZ)

    # INTERSECT

    # Call intersection function
    intCurveNames = manager.intersect(distTol=distTol)
    intCurveName = intCurveNames[0]

    # SPLIT

    # Split curves based on TE and LE curves
    optionsDict = {'splittingCurves' : [comp1.curves['curve_le'],
                                        comp1.curves['curve_te_upp'],
                                        comp1.curves['curve_te_low']]}

    splitCurveNames = manager.split_intCurve(intCurveName,
                                             optionsDict,
                                             criteria='curve')

    # REMESH

    # Find the highest z-coordinate of the entire intersection (vertical position)
    maxZ = -99999
    for curve in splitCurveNames:
        curr_maxZ = np.max(manager.intCurves[curve].coor[:,2])
        maxZ = max(maxZ, curr_maxZ)

    # Now we can identify and remesh each curve properly
    for curveName in splitCurveNames:

        # Get pointer to the curve object
        curve = manager.intCurves[curveName]

        # The trailing edge curve will have less nodes than the other ones
        if curve.numNodes < 20:

            # This is the trailing edge curve.
            # Just apply an uniform spacing
            optionsDict = {'nNewNodes':9}
            TE_curveName = manager.remesh_intCurve(curveName,optionsDict)

        else:

            # We have an upper or lower skin curve.
            # First let's identify if the curve is defined from
            # LE to TE or vice-versa

            curveCoor = curve.get_points()
            deltaX = curveCoor[-1,0] - curveCoor[0,0]

            if deltaX > 0:
                LE_to_TE = True
            else:
                LE_to_TE = False

            # Compute the highest vertical coordinate of the curve
            curr_maxZ = np.max(curve.coor[:,2])

            # Now we can determine if we have upper or lower skin
            if curr_maxZ < maxZ:

                # We are at the lower skin

                if LE_to_TE:

                    optionsDict = {'nNewNodes':numSkinNodes,
                                   'spacing':'hypTan',
                                   'initialSpacing':LE_spacing,
                                   'finalSpacing':TE_spacing}

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {'nNewNodes':numSkinNodes,
                                   'spacing':'hypTan',
                                   'initialSpacing':TE_spacing,
                                   'finalSpacing':LE_spacing}

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

            else:

                # We are at the upper skin

                if LE_to_TE:

                    optionsDict = {'nNewNodes':numSkinNodes,
                                   'spacing':'hypTan',
                                   'initialSpacing':LE_spacing,
                                   'finalSpacing':TE_spacing}

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {'nNewNodes':numSkinNodes,
                                   'spacing':'hypTan',
                                   'initialSpacing':TE_spacing,
                                   'finalSpacing':LE_spacing}

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

    # Now we can merge the new curves
    curveNames = [TE_curveName, LS_curveName, US_curveName]
    mergedCurveName = 'intersection'
    manager.merge_intCurves(curveNames, mergedCurveName)

    # REORDER
    manager.intCurves[mergedCurveName].shift_end_nodes(criteria='maxX')

    # Flip the curve for marching if necessary
    mergedCurveCoor = manager.intCurves[mergedCurveName].get_points()
    deltaZ = mergedCurveCoor[1,2] - mergedCurveCoor[0,2]

    if deltaZ > 0:
        manager.intCurves[mergedCurveName].flip()

    # Export final curve
    manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName+'_%03d'%numPasses)

    # MARCH SURFACE MESHES
    meshName = 'mesh'

    options_body = {

        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.01,
        'numLayers' : 49,
        'extension' : 1.3,
        'epsE0' : 4.5,
        'theta' : -0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'nuArea' : 0.16,
        'numAreaPasses' : 20,
        'sigmaSplay' : 0.3,
        'cMax' : 5.0,
        'ratioGuess' : 10.0,

    }

    options_wing = {

        'bc1' : 'curve:curve_te_upp',
        'bc2' : 'curve:curve_te_upp',
        'dStart' : 0.01,
        'numLayers' : 49,
        'extension' : 1.4,
        'epsE0' : 12.5,
        'theta' : -0.8,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 4,
        'nuArea' : 0.16,
        'numAreaPasses' : 0,
        'sigmaSplay' : 0.3,
        'cMax' : 10.0,
        'ratioGuess' : 10.0,
        'guideCurves' : ['curve_te_low'],
        'remesh' : True

    }

    meshNames = manager.march_intCurve_surfaceMesh(mergedCurveName, options0=options_wing, options1=options_body, meshName=meshName)

    for meshName in meshNames:
        manager.meshGenerators[meshName].export_plot3d(meshName+'_%03d'%numPasses+'.xyz')

    # DERIVATIVE SEEDS
    '''
    # Get spanwise position of the tracked node
    Y = manager.meshGenerators[meshNames[0]].meshObj.get_points()[1, i_node]

    # Set up derivative seeds so we can compute the derivative of the spanwise position
    # of the tracked node of the intersection
    meshb = np.zeros(manager.meshes[meshNames[0]].mesh.shape)
    meshb[1, i_node, j_node] = 1.0

    # Set derivatives to the mesh
    manager.meshes[meshNames[0]].set_reverseADSeeds(meshb)

    # REVERSE AD

    # Call AD code
    manager.reverseAD()

    # Get relevant seeds
    coor1b, curveCoor1b = manager.geoms[name1].get_reverseADSeeds()
    coor2b, curveCoor2b = manager.geoms[name2].get_reverseADSeeds()

    # Condense derivatives to take into acount the translation
    dYdZ = 0.0
    dYdZ = np.sum(coor1b[:,2])
    for curveName in curveCoor1b:
        dYdZ = dYdZ + np.sum(curveCoor1b[curveName][:,2])

    # Move the wing back to its original position
    manager.geoms[name1].translate(0.0, 0.0, -wing_deltaZ)
    '''

    Y = 0
    dYdZ = 0

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
numPasses = 0
for ii in range(numPos):
    print ''
    print 'translation'
    print deltaZ[ii]
    Y[ii], dYdZ[ii] = compute_position(deltaZ[ii])
    numPasses = numPasses + 1
    print 'results'
    print Y[ii]
    print dYdZ[ii]
    print ''

results = np.vstack([deltaZ,Y,dYdZ])
with open('results.pickle','w') as fid:
    pickle.dump(results,fid)
