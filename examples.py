# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, max
from numpy.random import rand
import surface
import os

# OPTIONS

# Select example
# example = 'kink_on_plate'
#example = 'line_on_curve'
example = 'line_on_paraboloid'
#example = 'line_on_eggcrate'
#example = 'circle_on_curve'
#example = 'circle_on_paraboloid'
#example = 'airfoil_on_cylinder'

# EXAMPLE SELECTION

if example == 'kink_on_plate':

    # Set problem
    rBaseline = array([[0,0,0],
                       [1,0,0],
                       [2,0,0],
                       [3,0,0],
                       [4,0,0],
                       [5,0,0]])
    bc1 = 'curve'
    bc2 = 'curve'
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
    sigmaSplay = 0.0
    cMax = 20.0
    ratioGuess = 20

    # Options
    sBaseline = 0.15
    numLayers = 20
    extension = 3.25

elif example == 'line_on_curve':

    # Set problem
    n = 80
    nums = 0.1*linspace(1, -1, n)
    x = -1.0*ones(n)
    y = nums
    z = 10*y**2
    rBaseline = vstack([x, y, z]).T
    bc1 = 'splay'
    bc2 = 'splay'

    # Define surface points
    s0 = 3
    s1 = 3
    nu, nv = 100, 100
    pts = zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            pts[i, j, 0] = s0 * (-1. + 2 * i / (nu-1))    # X coordinates
            pts[i, j, 1] = s1 * (-1. + 2 * j / (nv-1))    # Y coordinates
            pts[i, j, 2] = 10*(s1 * (-1. + 2 * j / (nv-1)))**2  # Z coordinates

    # Set parameters
    epsE0 = 2.5
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
    extension = 10.0

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
    s0 = 3
    s1 = 3
    radius = 1.0
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

# print curve1_pts[:, 0, :]

# Project curve points to the Surface
surf.inverse_evaluate(curve1_pts.reshape((-1, 3)))
curve1_pts = surf.get_points(None).reshape((-1, 4, 3))

# print curve1_pts[:, 0, :]
# exit()

surf.inverse_evaluate(curve2_pts.reshape((-1, 3)))
curve2_pts = surf.get_points(None).reshape((-1, 4, 3))

curve1 = surface.Surface(curve1_pts)
curve2 = surface.Surface(curve2_pts)


# Set reference geometry
ref_geom = {'surf':surf,
            'curve1':curve2,
            'curve2':curve1}

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

# EXPORT REFERENCE SURFACE

# Replace coordinates
mesh.mesh = zeros((3,pts.shape[0],pts.shape[1]))
mesh.mesh[0,:,:] = pts[:,:,0]
mesh.mesh[1,:,:] = pts[:,:,1]
mesh.mesh[2,:,:] = pts[:,:,2]

# Export wing
mesh.exportPlot3d('surf.xyz')

# Open tecplot
if example == 'line_on_eggcrate':
    os.system('tec360 layout_eggcrate.lay')
else:
    os.system('tec360 layout_mesh.lay')
