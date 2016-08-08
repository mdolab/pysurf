# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
from pysurf import TSurfComponent
import os
from mpi4py import MPI

# OPTIONS

# Select example
# example = 'kink_on_plate'
# example = 'line_on_cylinder'
# example = 'cylinder_cap'
# example = 'line_on_cubeAndCylinder'
example = 'line_on_cube'

# EXAMPLE SELECTION

if example == 'kink_on_plate':

    # Read inputs from CGNS file
    geom = TSurfComponent('inputs/plate.cgns')

    # Set source curve
    curve = array([[0,0,0],
                   [1,1,0],
                   [2,2,0],
                   [3,2,0],
                   [4,1,0],
                   [5,0,0]])

    # Flip curve to change marching direction
    curve = curve[::-1,:]

    # Define boundary conditions
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 1.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.5
    cMax = 20.0
    ratioGuess = 20

    # Options
    sBaseline = 0.15
    numLayers = 20
    extension = 3.5

    # Give layout file
    layout_file = 'layout_plate.lay'

elif example == 'line_on_cylinder':

    # Read inputs from CGNS file
    geom = TSurfComponent('inputs/cylinder.cgns')

    # Flip BC curve
    geom.Curves['bc1'].flip()

    # Set problem
    curve = 'source'
    bc1 = 'splay'
    bc2 = 'curve:bc1'

    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.2
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
    numLayers = 50
    extension = 4.8

    layout_file = 'layout_cylinder.lay'

elif example == 'cylinder_cap':

    # Read inputs from CGNS file
    geom = TSurfComponent('inputs/cylinder.cgns')

    # Set reference curve
    n1 = 7
    n2 = 21

    x = hstack([linspace(0.5,0,n1), zeros(n2-1), linspace(0,0.5,n1)[1:]])
    y = hstack([-0.5*ones(n1), linspace(-0.5, 0.5, n2)[1:], 0.5*ones(n1-1)])
    z = zeros(n1+n2+n1-2)
    curve = vstack([x, y, z]).T

    # Set boundary conditions
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 2.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.2
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
    numLayers = 40
    extension = 2.5

    layout_file = 'layout_cylinder.lay'

elif example == 'line_on_cubeAndCylinder':

    # Read inputs from CGNS file
    #geom = TSurfComponent('inputs/cubeAndCylinder.cgns',['cylinder']) # mesh cylinder only
    geom = TSurfComponent('inputs/cubeAndCylinder.cgns',['geom']) # mesh cube only

    # Set problem
    curve = 'diag'
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 1.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.1
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
    numLayers = 50
    extension = 4.8

    layout_file = 'layout_cubeAndCylinder.lay'

elif example == 'line_on_cube':

    # Read inputs from CGNS file
    geom = TSurfComponent('inputs/cube.cgns') # mesh cube only

    # Set problem
    n = 40
    nums = linspace(0.75, 0.25, n)
    x = nums
    y = nums #0.2*ones(n)
    z = zeros(n)
    curve = vstack([x, y, z]).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 2.4
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.2
    cMax = 10.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
    numLayers = 60
    extension = 5.0

    layout_file = 'layout_cube.lay'

###########################################

# MAIN PROGRAM

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

mesh.exportPlot3d('output.xyz')

# EXPORT REFERENCE SURFACE

# Open tecplot
os.system('tec360 ' + layout_file)