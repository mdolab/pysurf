# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack
from numpy.random import rand
import surface
import os

# OPTIONS

# Select example
# example = 'kink_on_plate'
# example = 'line_on_curve'
example = 'line_on_paraboloid'
# example = 'circle_on_curve'
#example = 'circle_on_paraboloid'
example = 'airfoil_on_cylinder'

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
    numLayers = 3
    extension = 3.25
    n = len(rBaseline)/3

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)

    # Set parameters
    epsE0 = 1.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.0
    cMax = 20.0
    ratioGuess = 20

    # Options
    sBaseline = 0.15
    numLayers = 20
    extension = 3.25

elif example == 'line_on_curve':

    # Set problem
    n = 20
    nums = linspace(1, -1, n)
    x = zeros(n)
    y = nums
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T.ravel()
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)    # X coordinates
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)    # Y coordinates
            pts[i, j, 2] = 100*s1 * (-1. + 2 * j / nv)**2  # Z coordinates

    # Set parameters
    epsE0 = 5.0
    theta = 3.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 10
    sigmaSplay = 0.3
    cMax = 0.5
    ratioGuess = 20

    # Options
    sBaseline = 0.002
    numLayers = 50
    extension = 5

elif example == 'line_on_paraboloid':

    # Set problem
    n = 20
    nums = linspace(1, -1, n)
    x = -0.5*ones(n)
    y = nums
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)    # X coordinates
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)    # Y coordinates
            pts[i, j, 2] = 100*s0 * (-1. + 2 * i / nv)**2 + 100*s1 * (-1. + 2 * j / nv)**2  # Z coordinates

    # Set parameters
    epsE0 = 20.0
    theta = 5.0
    alphaP0 = 0.25
    numSmoothingPasses = 1
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.8
    cMax = 100.0
    ratioGuess = 2

    # Options
    sBaseline = 0.05
    numLayers = 20
    extension = 10

elif example == 'circle_on_curve':

    # Set problem
    n = 41
    nums = linspace(1, 0, n) * 2*pi
    x = cos(nums)
    y = sin(nums)
    z = zeros(n)
    rBaseline = vstack([x, y, z])
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
    cMax = 1.0
    ratioGuess = 20

    # Options
    sBaseline = 0.02
    numLayers = 40
    extension = 10

elif example == 'airfoil_on_cylinder':

    # Set problem
    n = 81
    tc = 0.12
    x = (1-cos(linspace(0, 1, n)*pi))/2
    y = 5*tc*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
    x = hstack([x[::-1], x[1:]])
    y = hstack([-y[::-1], y[1:]])
    z = ones(len(x))
    rBaseline = vstack([x, y, z]).T

    bc1 = 'continuous'
    bc2 = 'continuous'

    # Define surface points
    s0 = 100
    s1 = 100
    radius = 1.5
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            theta = pi*j/(nv-1)
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)
            pts[i, j, 1] = -radius*cos(theta)
            pts[i, j, 2] = radius*sin(theta)

    # Set parameters
    epsE0 = 4.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 10.0
    ratioGuess = 20

    # Options
    sBaseline = 0.02
    numLayers = 10
    extension = 2

###########################################

# MAIN PROGRAM

# Create surface object
surf = surface.Surface(pts)

# Set reference geometry
ref_geom = {'surf':surf}

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

mesh = hypsurf.HypSurfMesh(curve=rBaseline, ref_geom=ref_geom, options=options)

mesh.createMesh()

mesh.exportPlot3d('output.xyz')

os.system('tec360 layout_mesh.lay')
