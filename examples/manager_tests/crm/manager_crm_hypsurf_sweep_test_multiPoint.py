# IMPORTS
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
deltaZ = np.linspace(0.000001, 3.0, 11)

# TRACKING POINT
# Give i and j coordinates that we will use to create the tracking slice
ij_list = [[190,  0],
           [190, 10],
           [190, 20],
           [190, 30],
           [190, 40],
           [190, 48]]

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
    manager.clear_all()

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
    manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName)

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
        manager.meshGenerators[meshName].export_plot3d(meshName+'.xyz')

    # DERIVATIVE SEEDS

    # Initialize arrays to hold results
    Y = np.zeros(len(ij_list))
    dYdZ = np.zeros(len(ij_list))

    # We will run one node position at a time
    for nodeID,ij in enumerate(ij_list):

        # Get coordinates
        i_node = ij[0]
        j_node = ij[1]

        # Compute node position in the flattened array
        nodePos = (2*numSkinNodes + 9 - 2)*j_node + i_node

        # Get spanwise position of the tracked node
        Y[nodeID] = manager.meshGenerators[meshNames[0]].meshObj.get_points()[nodePos,1]

        # Set up derivative seeds so we can compute the derivative of the spanwise position
        # of the tracked node of the intersection
        meshb = np.zeros(manager.meshGenerators[meshNames[0]].meshObj.get_points().shape)
        meshb[nodePos, 1] = 1.0

        # Set derivatives to the mesh
        manager.meshGenerators[meshNames[0]].meshObj.set_reverseADSeeds(meshb)

        # REVERSE AD
        
        # Call AD code
        manager.reverseAD()
        
        # Get relevant seeds 
        coor1b, curveCoor1b = manager.geoms[name1].get_reverseADSeeds()
        coor2b, curveCoor2b = manager.geoms[name2].get_reverseADSeeds()
        
        # Condense derivatives to take into acount the translation
        dYdZ[nodeID] = 0.0
        dYdZ[nodeID] = dYdZ[nodeID] + np.sum(coor1b[:,2])
        for curveName in curveCoor1b:
            dYdZ[nodeID] = dYdZ[nodeID] + np.sum(curveCoor1b[curveName][:,2])

        # Clean derivatives for a new pass
        manager.clean_reverseADSeeds()

    # Move the wing back to its original position
    manager.geoms[name1].translate(0.0, 0.0, -wing_deltaZ)

    # Return the spanwise position of the highest intersection node and its derivative
    return Y, dYdZ

# END OF compute_position
#======================================================

# Now we will execute the function for every wing position
for ii in range(len(deltaZ)):
    print('')
    print('translation')
    print(deltaZ[ii])
    print('')
    Y, dYdZ = compute_position(deltaZ[ii])
    print('')
    print('results')
    print(Y)
    print(dYdZ)
    print('')

    results = [deltaZ[ii], np.vstack([Y,dYdZ])]
    with open('results_%03d.pickle'%(ii),'w') as fid:
        pickle.dump(results,fid)
