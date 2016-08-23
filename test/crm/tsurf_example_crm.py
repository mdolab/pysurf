# IMPORTS

from __future__ import division
import numpy as np
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
from pysurf import TSurfComponent, compute_intersections, plot3d_interface, hypsurf
import os
from mpi4py import MPI
import copy

# USER INPUTS

# Define translation cases for the wing
nStates = 20

#wingTranslation = [[0.0, -3.0, s] for s in np.hstack([np.linspace(0,0.25,nStates),np.linspace(0.25,0,nStates)])]
#wingTranslation = [[0.0, -100.0, 145.0]]

wingTranslation = np.zeros((nStates,3))
wingTranslation[:,1] = np.linspace(-10.0, -100.0, nStates)
wingTranslation[:,2] = np.linspace(0.0, 140.0, nStates)

wingRotation = [0.0 for s in range(len(wingTranslation))]

# Define video parameters
fps = 2

# Load components
wing = TSurfComponent('../../inputs/crm.cgns',['w_upp','w_low','te_low_curve','te_upp_curve'])
body = TSurfComponent('../../inputs/crm.cgns',['b_fwd','b_cnt','b_rrf'])

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

        mesh = hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

        mesh.createMesh()

        mesh.exportPlot3d(output_name)


    #================================================

    # Apply transformations to the wing
    wing.translate(wingTranslation[0],wingTranslation[1],wingTranslation[2])
    wing.Curves['te_low_curve'].translate(wingTranslation[0],wingTranslation[1],wingTranslation[2])
    wing.Curves['te_upp_curve'].translate(wingTranslation[0],wingTranslation[1],wingTranslation[2])

    # Compute intersection
    Intersections = compute_intersections([wing, body])

    quit()

    # Reorder curve so that it starts at the trailing edge
    Intersections[0].shift_end_nodes(criteria='maxX')

    # Remesh curve to get better spacing
    Intersections[0].remesh(nNewNodes=300)

    # Check if we need to flip the curve
    if Intersections[0].coor[2,0] > Intersections[0].coor[2,-1]:
        Intersections[0].flip()

    curve = plot3d_interface.Curve(Intersections[0].coor, Intersections[0].barsConn)
    curve.export_plot3d('curve')

    # Add intersection curve to each component
    wing.add_curve('intersection',Intersections[0])
    body.add_curve('intersection',Intersections[0])

    # GENERATE WING MESH

    if generate_wing:

        print 'Generating wing mesh'

        # Flip BC curve
        #wing.Curves['intersection'].flip()
        #wing.Curves['te'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'curve:te_low_curve'
        bc2 = 'curve:te_upp_curve'

        # Set parameters
        epsE0 = 9.5
        theta = 1.0
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 0
        sigmaSplay = 0.
        cMax = 10.0
        ratioGuess = 20

        # Options
        sBaseline = 3.5
        numLayers = 30
        extension = 1.5

        # Call meshing function
        generateMesh(curve, wing, 'wing.xyz')

    # GENERATE BODY MESH

    if generate_body:

        print 'Generating body mesh'

        # Flip BC curve
        body.Curves['intersection'].flip()

        # Set problem
        curve = 'intersection'
        bc1 = 'splay'
        bc2 = 'splay'

        # Set parameters
        epsE0 = 4.5
        theta = 0.0
        alphaP0 = 0.25
        numSmoothingPasses = 0
        nuArea = 0.16
        numAreaPasses = 0
        sigmaSplay = 0.3
        cMax = 1.0
        ratioGuess = 20

        # Options
        sBaseline = 3.0
        numLayers = 20
        extension = 1.3

        # Call meshing function
        generateMesh(curve, body, 'body.xyz')


    # Run tecplot in bash mode to save picture
    os.system('tec360 -b plotMesh.mcr')

    # Rename image file
    os.system('mv images/image.png images/image%03d.png'%(meshIndex))


    # Revert transformations to the wing
    wing.translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])
    wing.Curves['te_low_curve'].translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])
    wing.Curves['te_upp_curve'].translate(-wingTranslation[0],-wingTranslation[1],-wingTranslation[2])

# end of genereteWingBodyMesh
    ###########################################

# Generate mesh for each case

for index in range(len(wingTranslation)):

    # Generate mesh for each case
    try:
        generateWingBodyMesh(wingTranslation[index],wingRotation[index],index+1)
    except:
        pass

#Generate a video with all the images using the ffmpeg command
os.system('ffmpeg -f image2 -r ' + str(fps) + ' -i images/image%03d.png -vcodec mpeg4 -y movie.mp4')
