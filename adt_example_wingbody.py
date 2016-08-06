# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
from pysurf import ADTComponent
import classes
import os
from mpi4py import MPI
import copy

# OPTIONS
generate_wing = True
generate_body = True
layout_file = 'layout_wingBody_adt.lay'

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


# GENERATE WING MESH

if generate_wing:

    # Read inputs from CGNS file
    wing = ADTComponent('inputs/wingBody.cgns',['wing'])

    # Flip BC curve
    wing.Curves['intersection'].flip()
    wing.Curves['lower_te'].flip()
    wing.Curves['upper_te'].flip()

    # Set problem
    curve = 'intersection'
    bc1 = 'curve:upper_te'
    bc2 = 'curve:lower_te'

    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = .005
    numLayers = 20
    extension = 2.5

    # Call meshing function
    generateMesh(curve, wing, 'wing.xyz')

# GENERATE BODY MESH

if generate_body:

    # Read inputs from CGNS file
    body = ADTComponent('inputs/wingBody.cgns',['geom'])

    # Set problem
    curve = 'intersection'
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = .005
    numLayers = 10
    extension = 1.5

    # Call meshing function
    generateMesh(curve, body, 'body.xyz')


###########################################

# Open tecplot
os.system('tec360 ' + layout_file)
