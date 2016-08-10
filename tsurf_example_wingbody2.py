# IMPORTS

from __future__ import division
import hypsurf
import numpy as np
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
# from surface import Surface  # GeoMACH
from pysurf import TSurfComponent, compute_intersections, plot3d_interface
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

# Load components
wing = TSurfComponent('inputs/wingBody.cgns',['wing','upper_te','lower_te', 'te'])
body = TSurfComponent('inputs/wingBody.cgns',['geom'])

#wing = TSurfComponent('inputs/cubeAndCylinder.cgns',['cylinder'])
#body = TSurfComponent('inputs/cubeAndCylinder.cgns',['geom'])

# Compute intersection
Intersections = compute_intersections([wing, body])

curve = plot3d_interface.Curve(Intersections[0].coor, Intersections[0].barsConn)
curve.export_plot3d('curve_before')

# Reorder curve so that it starts at the trailing edge
Intersections[0].shift_end_nodes(criteria='maxX')

# Remesh curve to get better spacing
Intersections[0].remesh(nNewNodes=300)

curve = plot3d_interface.Curve(Intersections[0].coor, Intersections[0].barsConn)
curve.export_plot3d('curve_after')

#quit()

# Add intersection curve to each component
wing.add_curve('intersection',Intersections[0])
body.add_curve('intersection',Intersections[0])

# GENERATE WING MESH

if generate_wing:

    # Flip BC curve
    wing.Curves['intersection'].flip()
    wing.Curves['te'].flip()
    #wing.Curves['lower_te'].flip()
    #wing.Curves['upper_te'].flip()

    # Set problem
    curve = 'intersection'
    bc1 = 'curve:te'
    bc2 = 'curve:te'

    # Set parameters
    epsE0 = 5.5#5.5
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
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 10.0
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
