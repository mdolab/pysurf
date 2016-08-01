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

# Select example
#example = 'kink_on_plate'
# example = 'line_on_cylinder'
#example = 'line_on_cubeAndCylinder'
#example = 'line_on_cube'
# example = 'wing'
example = 'body'


# EXAMPLE SELECTION

if example == 'kink_on_plate':

    # Set problem
    rBaseline = array([[0,0,0],
                       [1,1,0],
                       [2,2,0],
                       [3,2,0],
                       [4,1,0],
                       [5,0,0]])
    bc1 = 'splay'
    bc2 = 'splay'
    n = len(rBaseline)/3

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / (nu-1))
            pts[i, j, 1] = s1 * (-1. + 2 * j / (nv-1))

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 4
    curve1_pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            curve1_pts[i, j, 0] = 0.0
            curve1_pts[i, j, 1] = s1 * (-1. + 2 * i / (nv-1))

    # Define surface points
    s0 = 50
    s1 = 50
    nu, nv = 5, 4
    curve2_pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            curve2_pts[i, j, 0] = 5.0
            curve2_pts[i, j, 1] = s1 * (-1. + 2 * i / (nv-1))

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
    extension = 2.5

    surf = Surface('inputs/plate.cgns')  # ADT

elif example == 'line_on_cylinder':

    # Read inputs from CGNS file
    componentsDict = adt_geometry.getCGNScomponents('inputs/cylinder.cgns')

    # Flip BC curve
    componentsDict['bc1'].flip()

    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['source']
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

elif example == 'line_on_cubeAndCylinder':

    # Read inputs from CGNS file
    componentsDict = adt_geometry.getCGNScomponents('inputs/cubeAndCylinder.cgns')

    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['diag']
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

    # Set problem
    n = 40
    nums = linspace(0.75, 0.25, n)
    x = nums
    y = nums #0.2*ones(n)
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Set parameters
    epsE0 = 4.0
    theta = 1.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.0
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
    numLayers = 50
    extension = 5.0

    surf = Surface('inputs/cube.cgns')  # ADT

    layout_file = 'layout_cube.lay'

elif example == 'line_on_paraboloid':

    # Set problem
    n = 20
    nums = 2.5*linspace(-1, 1, n)
    x = -5.0*ones(n)
    y = nums
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T
    bc1 = 'curve'
    bc2 = 'curve'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)
            pts[i, j, 2] = 10*s0 * (-1. + 2 * i / nv)**2 + 10*s1 * (-1. + 2 * j / nv)**2

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 10, 4
    curve1_pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            curve1_pts[i, j, 0] = -10. + (nu - i) - 1.
            curve1_pts[i, j, 1] = 2.5

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 10, 4
    curve2_pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            curve2_pts[i, j, 0] = -10. + (nu - i) - 1.
            curve2_pts[i, j, 1] = -2.5

    # Set parameters
    epsE0 = 20.0
    theta = 5.0
    alphaP0 = 0.25
    numSmoothingPasses = 1
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.2
    cMax = 1.0
    ratioGuess = 2

    # Options
    sBaseline = 0.05
    numLayers = 15
    extension = 10

elif example == 'line_on_eggcrate':

    parAmplitude = 0.0
    sinAmplitude = 0.5

    # Set problem
    n = 40
    nums = linspace(3, -3, n)
    x = -2.5*ones(n)
    y = nums
    z = parAmplitude*(x**2 + y**2) + \
        sinAmplitude*(sin(x)**2 + sin(y)**2)# Z coordinates
    rBaseline = vstack([x, y, z]).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 10
    s1 = 10
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0*(-1. + 2 * i / (nu-1))    # X coordinates
            pts[i, j, 1] = s1*(-1. + 2 * j / (nv-1))    # Y coordinates
            pts[i, j, 2] = parAmplitude*(pts[i, j, 0]**2 + pts[i, j, 1]**2) + \
                           sinAmplitude*(sin(pts[i, j, 0])**2 + sin(pts[i, j, 1])**2)# Z coordinates

    # Set parameters
    epsE0 = 5.0
    theta = 1.0
    alphaP0 = 0.25
    numSmoothingPasses = 1
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.1
    cMax = 100.0
    ratioGuess = 4

    # Options
    sBaseline = 0.05
    numLayers = 40
    extension = 4.0

elif example == 'circle_on_curve':

    # Set problem
    n = 41
    nums = linspace(1, 0, n) * 2*pi
    x = cos(nums)
    y = sin(nums)
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T
    bc1 = 'continuous'
    bc2 = 'continuous'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)
            pts[i, j, 2] = 10*s1 * (-1. + 2 * j / nv)**2

    # Set parameters
    epsE0 = 1.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.02
    numLayers = 40
    extension = 10

elif example == 'circle_on_paraboloid':

    # Set problem
    n = 41
    nums = linspace(1, 0, n) * 2*pi
    x = cos(nums)+5.0
    y = sin(nums)
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T
    bc1 = 'continuous'
    bc2 = 'continuous'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)
            pts[i, j, 2] = 10*s0 * (-1. + 2 * i / nv)**2 + 10*s1 * (-1. + 2 * j / nv)**2

    # Set parameters
    epsE0 = 1.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 100.0
    ratioGuess = 3

    # Options
    sBaseline = 0.02
    numLayers = 100
    extension = 30

elif example == 'wing':

    # Read inputs from CGNS file
    componentsDict = adt_geometry.getCGNScomponents('inputs/wing.cgns')
    componentsDict.pop("geom")

    # Flip BC curve
    componentsDict['intersection'].flip()
    componentsDict['lower_te'].flip()
    componentsDict['upper_te'].flip()

    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['intersection']
    bc1 = 'curve:lower_te'
    bc2 = 'curve:upper_te'

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

    layout_file = 'wingbody.lay'

elif example == 'body':

    # Read inputs from CGNS file
    componentsDict = adt_geometry.getCGNScomponents('inputs/wing.cgns')
    componentsDict.pop("wing")

    # Flip BC curve
    # componentsDict['bc1'].flip()

    # Assign components to a Geometry object
    geom = classes.Geometry(componentsDict)

    # Set problem
    curve = componentsDict['intersection']
    bc1 = 'splay'  #'curve:UPPER_TE'
    bc2 = 'splay'  #'curve:LOWER_TE'

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

    layout_file = 'wingbody.lay'


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
