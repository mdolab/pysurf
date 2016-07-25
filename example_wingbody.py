# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, ceil
from numpy.random import rand
import surface
import os

trailingEdge = 'sharp'
generateBodyMesh = True
generateWingMesh = True

# BODY GEOMETRY

# Basic dimensions
x0 = -0.5 # Initial position
length = 2.0 # cylinder length
radius = 0.5 # cylinder ratio

# Set cylinder points
nu, nv = 100, 100
pts = zeros((nu, nv, 3))
for i in xrange(nu):
    for j in xrange(nv):
        theta = pi*j/(nv-1)
        pts[i, j, 0] = x0 + length*i/(nu-1)
        pts[i, j, 1] = -radius*cos(theta)
        pts[i, j, 2] = radius*sin(theta)

# Save body points
bodyPts = pts[:,:,:]

# Create surface object
bodySurf = surface.Surface(bodyPts)

# AIRFOIL COORDINATES
# Parameters
n = 81 # Number of points in the airfoil
tc = 0.12 # Airfoil thickness
chord = 1.0 # Airfoil chord
x0 = 0.0 # Leading edge position
y0 = -.3 # Leading edge position
z0 = radius # Leading edge position

# Generate non-dimensional coordinates
if trailingEdge == 'sharp':
    x = (1-cos(linspace(0, 1, n)*pi))/2
    y = 5*tc*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)

elif trailingEdge == 'blunt':
    x = (1-cos(linspace(0, 1, n)*pi))/2
    y = 5*tc*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

    # Estimate spacing between last two nodes
    dy = sqrt((x[-1] - x[-2])**2 + (y[-1] - y[-2])**2)

    # Estimate how many nodes we can fit on the blunt TE (symmetric part only)
    nTEnodes = ceil(y[-1]/dy) + 1

    # Generate coordinates for the vertical segment
    yTE = linspace(y[-1],0,nTEnodes)
    xTE = ones(nTEnodes)

    # Append the coordinates
    x = hstack([x, xTE[1:]])
    y = hstack([y, yTE[1:]])

# Scale and translate
x = x0 + hstack([x[::-1], x[1:]])*chord
y = y0 + hstack([-y[::-1], y[1:]])*chord
z = z0*ones(len(x))

# Assemble coordinates
rBaseline = vstack([x, y, z]).T

# Project to the cylinder surface
bodySurf.inverse_evaluate(rBaseline)
rBaseline = bodySurf.get_points(None)

# WING COORDINATES
# We will use the projected airfoil to extreude the wing

# Parameters
span = 1.0
nv = 4 # Number of span-wise control points

# Determine number of chord-wise points
nu = rBaseline.shape[0]

# Determine points
pts = zeros((nu, nv, 3))
for i in xrange(nu):
    for j in xrange(nv):
        pts[i, j, 0] = rBaseline[i,0]
        pts[i, j, 1] = rBaseline[i,1]
        pts[i, j, 2] = rBaseline[i,2] + span*(j/(nv-1) - 0.05)

# Save wing points
wingPts = pts[:,:,:]

# Create surface object
wingSurf = surface.Surface(wingPts)

#==============================================================#

# BODY PROBLEM

# Set reference geometry
ref_geom = {'surf':bodySurf}

# Set options
if trailingEdge == 'blunt':
    options = {
        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.02,
        'numLayers' : 17,
        'extension' : 1.7,
        'epsE0' : 10.0,
        'theta' : 0.0,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'numAreaPasses' : 0,
        'cMax' : 100.0,
    }

if trailingEdge == 'sharp':
    options = {
        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.02,
        'numLayers' : 17,
        'extension' : 1.7,
        'epsE0' : 1.5,
        'theta' : 0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 5,
        'numAreaPasses' : 5,
        'cMax' : 1.0,
    }

# Set meshing problem
bodyMesh = hypsurf.HypSurfMesh(curve=rBaseline, ref_geom=ref_geom, options=options)

if generateBodyMesh:

    # Print log
    print 'Generating body mesh'

    # Generate mesh
    bodyMesh.createMesh()

    # Export mesh
    bodyMesh.exportPlot3d('bodyMesh.xyz')

# WING PROBLEM

# Determine points
nu = 6
nv = 4
curve1_pts = zeros((nu, nv, 3))
for i in xrange(nu):
    for j in xrange(nv):
        curve1_pts[i, j, 0] = rBaseline[-1,0]
        curve1_pts[i, j, 1] = rBaseline[-1,1]
        curve1_pts[i, j, 2] = rBaseline[-1,2] + span*(i/(nv-1) - 0.05)

curve1 = surface.Surface(curve1_pts)

# Set reference geometry
ref_geom = {'surf':wingSurf,
            'curve1':curve1,
            'curve2':curve1}

# Set options
if trailingEdge == 'blunt':
    options = {
        'bc1' : 'continuous',
        'bc2' : 'continuous',
        'dStart' : 0.02,
        'numLayers' : 17,
        'extension' : 2,
        'epsE0' : 2.5,
        'theta' : 0.0,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'numAreaPasses' : 0,
        'cMax' : 1.0,
    }

if trailingEdge == 'sharp':
    options = {
        'bc1' : 'curve',
        'bc2' : 'curve',
        'dStart' : 0.002,
        'numLayers' : 4,
        'extension' : 1.05,
        'epsE0' : 1.5,
        'theta' : 1.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'numAreaPasses' : 0,
        'cMax' : 1.0,
    }

# Set meshing problem
wingMesh = hypsurf.HypSurfMesh(curve=rBaseline, ref_geom=ref_geom, options=options)

if generateWingMesh:

    # Print log
    print 'Generating wing mesh'

    # Generate mesh
    wingMesh.createMesh()

    # Export mesh
    wingMesh.exportPlot3d('wingMesh.xyz')

# EXPORTING REFERENCE SURFACES

# Replace coordinates of the body mesh
bodyMesh.mesh = zeros((3,bodyPts.shape[0],bodyPts.shape[1]))
bodyMesh.mesh[0,:,:] = bodyPts[:,:,0]
bodyMesh.mesh[1,:,:] = bodyPts[:,:,1]
bodyMesh.mesh[2,:,:] = bodyPts[:,:,2]

# Export body
bodyMesh.exportPlot3d('bodySurf.xyz')

# Replace coordinates of the wing mesh
wingMesh.mesh = zeros((3,wingPts.shape[0],wingPts.shape[1]))
wingMesh.mesh[0,:,:] = wingPts[:,:,0]
wingMesh.mesh[1,:,:] = wingPts[:,:,1]
wingMesh.mesh[2,:,:] = wingPts[:,:,2]

# Export wing
wingMesh.exportPlot3d('wingSurf.xyz')

# Open tecplot
os.system('tec360 layout_wingBody.lay')
