# IMPORTS

import numpy as np
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
import pysurf
import os
from mpi4py import MPI
import copy

# USER INPUTS
guideCurves = ['te_low_curve','te_upp_curve',]

# Define translation cases for the wing
nStates = 11

wingTranslation = np.zeros((nStates,3))
wingTranslation[:,1] = np.linspace(-10.0, -100.0, nStates)
wingTranslation[:,2] = np.linspace(0.0, 140.0, nStates)

wingTranslation = [wingTranslation[0,:]]

wingRotation = [0.0 for s in range(len(wingTranslation))]

# Define video parameters
fps = 2

# Load components
wing = pysurf.TSurfGeometry('../../inputs/crm.cgns',['w_upp','w_low','w_ted','te_low_curve','te_upp_curve','w_le_curve'])
body = pysurf.TSurfGeometry('../../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

# Flip some curves
wing.curves['te_upp_curve'].flip()

# Define function to generate a wing body mesh for a given wing position and angle

def generateWingBodyMesh(wingTranslation, wingRotation, meshIndex):

    # OPTIONS
    generate_wing = False
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

        mesh = pysurf.hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

        mesh.createMesh()

        mesh.exportPlot3d(output_name)


    #================================================

    # Apply transformations to the wing
    print(wingTranslation)
    wing.translate(wingTranslation[0],wingTranslation[1],wingTranslation[2])

    # Compute intersection
    Intersections = pysurf.tsurf_tools.compute_intersections([wing,body])

    # Reorder curve so that it starts at the trailing edge
    for intName in Intersections:
        Intersections[intName].shift_end_nodes(criteria='maxX')

    # Split curves
    pysurf.tsurf_tools.split_curves(Intersections)
    pysurf.tsurf_tools.split_curve_with_curve(Intersections, Intersections.keys()[1], wing.curves['w_le_curve'])

    # Remesh curve to get better spacing
    curveNames = Intersections.keys()
    Intersections[curveNames[0]] = Intersections[curveNames[0]].remesh(nNewNodes=5)
    Intersections[curveNames[1]] = Intersections[curveNames[1]].remesh(nNewNodes=150, spacing='hypTan', initialSpacing=0.1, finalSpacing=.5)
    Intersections[curveNames[2]] = Intersections[curveNames[2]].remesh(nNewNodes=150, spacing='hypTan', initialSpacing=0.5, finalSpacing=.1)

    # Currently 1 is upper, 2 is lower surface

    # Check if we need to flip the curve
    # if Intersections[curveNames[1]].coor[2,0] > Intersections[curveNames[1]].coor[2,-1]:
    #     Intersections[curveNames[0]].flip()
    #     Intersections[curveNames[1]].flip()
    #     Intersections[curveNames[2]].flip()

    # Merge curves
    mergedCurve = pysurf.tsurf_tools.merge_curves(Intersections,'intersection')

    mergedCurve.shift_end_nodes(criteria='maxX')

    # Export intersection curve
    mergedCurve.export_plot3d('curve')

    # Add intersection curve to each component
    wing.add_curve('intersection',mergedCurve)
    body.add_curve('intersection',mergedCurve)

    # GENERATE WING MESH

    if generate_wing:

        print('Generating wing mesh')

        # Flip BC curve
        #wing.curves['intersection'].flip()
        #wing.curves['te'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'curve:te_upp_curve'
        bc2 = 'curve:te_upp_curve'

        # Set parameters
        epsE0 = 9.5
        theta = 1.0
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 0
        sigmaSplay = 0.
        cMax = 1000.0
        ratioGuess = 20

        # Options
        sBaseline = .2
        numLayers = 60
        extension = 1.4

        # Call meshing function
        generateMesh(curve, wing, 'wing.xyz')

    # GENERATE BODY MESH

    if generate_body:

        print('Generating body mesh')

        # Flip BC curve
        body.curves['intersection'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'continuous'
        bc2 = 'continuous'

        # Set parameters
        epsE0 = 6.5
        theta = 0.0
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 5
        sigmaSplay = 0.3
        cMax = 10.0
        ratioGuess = 20

        # Options
        sBaseline = 0.2
        numLayers = 60
        extension = 1.3

        # Call meshing function
        generateMesh(curve, body, 'body.xyz')


    # Run tecplot in bash mode to save picture
    # os.system('tec360 -b plotMesh.mcr')
    os.system('tec360 layout_mesh.lay')


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
os.system('ffmpeg -f image2 -r ' + str(fps) + ' -i images/image%03d.png -vcodec mpeg4 -y movie.mp4')
