# IMPORTS

from __future__ import division
import numpy as np
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
import pysurf
import os
from mpi4py import MPI
import copy

# USER INPUTS
guideCurves = ['curve_te_low']

# Define translation cases for the wing
nStates = 2

wingTranslation = [[0.0, 0.0, 0.0]]

wingRotation = [0.0 for s in range(len(wingTranslation))]

# Define video parameters
fps = 2

# Load components
wing = pysurf.TSurfGeometry('../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
body = pysurf.TSurfGeometry('../inputs/fuselage_crm4.cgns',['fuse'])

# Load TE curves and append them to the wing component
curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves('curve_te_upp.plt')
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('curve_te_low.plt')
curveName = curve_te_upp.keys()[0]
wing.add_curve(curveName, curve_te_upp[curveName])
curveName = curve_te_low.keys()[0]
wing.add_curve(curveName, curve_te_low[curveName])

# Flip some curves
wing.curves['curve_te_low'].flip()

# Define function to generate a wing body mesh for a given wing position and angle

def generateWingBodyMesh(wingTranslation, wingRotation, meshIndex):

    # OPTIONS
    generate_wing = True
    generate_body = True

    # DEFINE FUNCTION TO GENERATE MESH

    def generateMesh(curve, geom, output_name):

        # Set options
        options = {

            'bc1' : bc1,
            'bc2' : bc2,
            'dStart' : sBaseline,
            'numLayers' : numLayers,
            'extension' : extension,
            'epsE0' : epsE0,
            'theta' : theta,
            'alphaP0' : alphaP0,
            'numSmoothingPasses' : numSmoothingPasses,
            'nuArea' : nuArea,
            'numAreaPasses' : numAreaPasses,
            'sigmaSplay' : sigmaSplay,
            'cMax' : cMax,
            'ratioGuess' : ratioGuess,

        }

        # Set guideCurves for the wing case
        if output_name == 'wing.xyz':
            options['guideCurves'] = guideCurves
            options['remesh'] = True

        mesh = pysurf.hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

        mesh.createMesh()

        mesh.exportPlot3d(output_name)


    #================================================

    # Apply transformations to the wing
    print wingTranslation
    wing.translate(wingTranslation[0],wingTranslation[1],wingTranslation[2])

    # Compute intersection
    Intersections = pysurf.tsurf_tools.compute_intersections([wing,body])

    # Reorder curve so that it starts at the trailing edge
    for intName in Intersections:
        Intersections[intName].shift_end_nodes(criteria='maxX')

    # Split curves
    optionsDict = {'splittingCurves' : [wing.curves['curve_le'], wing.curves['curve_te_upp'], wing.curves['curve_te_low']]}

    optionsDict['curvesToBeSplit'] = [Intersections.keys()[0]]

    pysurf.tsurf_tools.split_curves(Intersections, optionsDict, criteria='curve')

    # Remesh curve to get better spacing
    # Dirty hard-coded logic that works for this specific case.
    # In the future we'll have a more sophisticated and robust method to sort
    # the curve results.
    curveNames = Intersections.keys()
    maxCoor = -99999.
    for cur in curveNames:
        curveMax = np.max(Intersections[cur].coor, axis=1)[2]
        if maxCoor < curveMax:
            maxCoor = curveMax

    for cur in curveNames:
        if Intersections[cur].coor.shape[1] < 20:
            Intersections[cur] = Intersections[cur].remesh(nNewNodes=9)
        else:
            end_to_end = Intersections[cur].extract_points()[0, :] - Intersections[cur].extract_points()[-1, :]
            if np.max(Intersections[cur].coor, axis=1)[2] < maxCoor:
                if end_to_end[0] > 0:
                    Intersections[cur] = Intersections[cur].remesh(nNewNodes=129, spacing='hypTan', initialSpacing=0.01, finalSpacing=.001)
                else:
                    Intersections[cur] = Intersections[cur].remesh(nNewNodes=129, spacing='hypTan', initialSpacing=0.001, finalSpacing=.01)
                    Intersections[cur].flip()
            else:
                if end_to_end[0] < 0:
                    Intersections[cur] = Intersections[cur].remesh(nNewNodes=129, spacing='hypTan', initialSpacing=0.001, finalSpacing=.01)
                else:
                    Intersections[cur] = Intersections[cur].remesh(nNewNodes=129, spacing='hypTan', initialSpacing=0.01, finalSpacing=.001)
                    Intersections[cur].flip()

    # Merge curves
    mergedCurve = pysurf.tsurf_tools.merge_curves(Intersections,'intersection')

    mergedCurve.shift_end_nodes(criteria='maxX')

    # Export intersection curve
    mergedCurve.export_tecplot('curve')

    # Add intersection curve to each component
    wing.add_curve('intersection',mergedCurve)
    body.add_curve('intersection',mergedCurve)

    # GENERATE WING MESH

    if generate_wing:

        print 'Generating wing mesh'

        # Flip BC curve
        #wing.curves['intersection'].flip()
        #wing.curves['te'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'curve:curve_te_upp'
        bc2 = 'curve:curve_te_upp'

        # Set parameters
        epsE0 = 12.5
        theta = -0.8
        alphaP0 = 0.25
        numSmoothingPasses = 4
        nuArea = 0.16
        numAreaPasses = 0
        sigmaSplay = 0.
        cMax = 10.0
        ratioGuess = 20

        # Options
        sBaseline = 0.01
        numLayers = 49
        extension = 1.4

        # Call meshing function
        generateMesh(curve, wing, 'wing.xyz')

    # GENERATE BODY MESH

    if generate_body:

        print 'Generating body mesh'

        # Flip BC curve
        body.curves['intersection'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'continuous'
        bc2 = 'continuous'

        # Set parameters
        epsE0 = 4.5
        theta = -0.5
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 20
        sigmaSplay = 0.3
        cMax = 10000.0
        ratioGuess = 20

        # Options
        sBaseline = 0.01
        numLayers = 49
        extension = 1.3

        # Call meshing function
        generateMesh(curve, body, 'body.xyz')


    # Run tecplot in bash mode to save picture
    # os.system('tec360 -b plotMesh.mcr')
    # os.system('tec360 layout_mesh.lay')


    # Rename image file
    # os.system('mv images/image.png images/image%03d.png'%(meshIndex))


    # Revert transformations to the wing
    wing.translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])
    #wing.curves['te_low_curve'].translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])
    #wing.curves['te_upp_curve'].translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])

# end of genereteWingBodyMesh
    ###########################################

# Generate mesh for each case

for index in range(len(wingTranslation)):

    # Generate mesh for each case
    #try:
    generateWingBodyMesh(wingTranslation[index],wingRotation[index],index+1)
    #except:
        #pass

#Generate a video with all the images using the ffmpeg command
#os.system('ffmpeg -f image2 -r ' + str(fps) + ' -i images/image%03d.png -vcodec mpeg4 -y movie.mp4')
