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
guideCurves = ['curve_te_upp','curve_te_low']

# Define translation cases for the wing
nStates = 11

wingTranslation = np.zeros((nStates,3))
wingTranslation[:,1] = np.linspace(-10.0, -100.0, nStates)
wingTranslation[:,2] = np.linspace(0.0, 140.0, nStates)

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
    generate_body = False

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
    curveNames = Intersections.keys()

    Intersections[curveNames[1]] = Intersections[curveNames[1]].remesh(nNewNodes=5)
    Intersections[curveNames[2]] = Intersections[curveNames[2]].remesh(nNewNodes=150, spacing='hypTan', initialSpacing=0.005, finalSpacing=.05)
    Intersections[curveNames[0]] = Intersections[curveNames[0]].remesh(nNewNodes=150, spacing='hypTan', initialSpacing=0.05, finalSpacing=.005)

    # Check if we need to flip the curve
    # if Intersections[curveNames[1]].coor[2,0] > Intersections[curveNames[1]].coor[2,-1]:
    #     Intersections[curveNames[0]].flip()
    #     Intersections[curveNames[1]].flip()
    #     Intersections[curveNames[2]].flip()

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
        epsE0 = 8.5
        theta = 0.0
        alphaP0 = 0.25
        numSmoothingPasses = 4
        nuArea = 0.16
        numAreaPasses = 0
        sigmaSplay = 0.
        cMax = 10.0
        ratioGuess = 20

        # Options
        sBaseline = 0.01
        numLayers = 60
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
        theta = 0.0
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 5
        sigmaSplay = 0.3
        cMax = 10.0
        ratioGuess = 20

        # Options
        sBaseline = 0.01
        numLayers = 60
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
