# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
from DiscreteSurf.python import adt_geometry
import classes
import os
from mpi4py import MPI

# OPTIONS
generate_wing = True
generate_body = False
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
    componentsDict = adt_geometry.getCGNScomponents('inputs/wingBody.cgns')
    componentsDict.pop("geom")

    # Flip BC curve
    componentsDict['intersection'].flip()
    componentsDict['lower_te'].flip()
    componentsDict['upper_te'].flip()
    
    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['intersection']
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
    extension = 3.5

    # Call meshing function
    generateMesh(curve, geom, 'wing.xyz')

# GENERATE BODY MESH

if generate_body:

    # Read inputs from CGNS file
    componentsDict = adt_geometry.getCGNScomponents('inputs/wingBody.cgns')
    componentsDict.pop("wing")

    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['intersection']
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
    generateMesh(curve, geom, 'body.xyz')


###########################################

# Open tecplot
os.system('tec360 ' + layout_file)