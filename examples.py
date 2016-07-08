# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones
from numpy.random import rand
import surface
import os

# OPTIONS

# Select example
#example = 'kink_on_plate'
#example = 'line_on_curve'
example = 'line_on_paraboloid'
#example = 'circle_on_curve'
#example = 'circle_on_paraboloid'

# EXAMPLE SELECTION

if example == 'kink_on_plate':

    # Set problem
    rBaseline = array([0,0,0, 1,1,0, 2,2,0, 3,2,0, 4,1,0, 5,0,0])
    NBaseline = array([[0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1]]).T
    bc1 = 'splay'
    bc2 = 'splay'
    sBaseline = 1e-2
    numLayers = 20
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
    cMax = 2.0
    ratioGuess = 20

    # Options
    sBaseline = 0.01
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
    NBaseline = array([array([0, 0, 1])] * n).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)    # X coordiantes
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)    # Y coordiantes
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
    rBaseline = vstack([x, y, z]).T.ravel()
    NBaseline = array([array([0, 0, 1])] * n).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / nu)    # X coordiantes
            pts[i, j, 1] = s1 * (-1. + 2 * j / nv)    # Y coordiantes
            pts[i, j, 2] = 100*s0 * (-1. + 2 * i / nv)**2 + 100*s1 * (-1. + 2 * j / nv)**2  # Z coordinates

    # Set parameters
    epsE0 = 15.0
    theta = 3.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 10
    sigmaSplay = 0.3
    cMax = 0.1
    ratioGuess = 20

    # Options
    sBaseline = 0.002
    numLayers = 50
    extension = 5

elif example == 'circle_on_curve':

    # Set problem
    n = 41
    nums = linspace(1, 0, n) * 2*pi
    x = cos(nums)
    y = sin(nums)
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T.ravel()
    NBaseline = array([array([0, 0, 1])] * n).T
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
    x = cos(nums)
    y = sin(nums)
    z = zeros(n)
    rBaseline = vstack([x, y, z]).T.ravel()
    NBaseline = array([array([0, 0, 1])] * n).T
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

###########################################

# MAIN PROGRAM

# Create surface object
surf = surface.Surface(pts)

# Project onto surface
surf.inverse_evaluate(rBaseline.reshape((n, 3)))
rBaseline = surf.get_points(None).flatten()

# Generate surface mesh
X,Y,Z = hypsurf.createMesh(rBaseline, bc1, bc2, sBaseline, numLayers, extension, surf,
                           epsE0 = epsE0, theta = theta,
                           alphaP0 = alphaP0, numSmoothingPasses = numSmoothingPasses,
                           nuArea = nuArea, numAreaPasses = numAreaPasses,
                           sigmaSplay = sigmaSplay,
                           cMax = cMax,
                           ratioGuess = ratioGuess)

# Plot results
fig = hypsurf.plotGrid(X,Y,Z,show=False)
fig.savefig('output.png')
hypsurf.exportPlot3d(X,Y,Z,'output.xyz')
os.system('tec360 layout_mesh.lay')
