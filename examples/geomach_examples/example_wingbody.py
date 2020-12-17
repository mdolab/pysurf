# IMPORTS

import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones, sqrt, hstack, ceil
from numpy.random import rand
from GeoMACHSurf.geomach_geometry import Surface, Curve
import os

trailingEdge = "blunt"
generateBodyMesh = True
generateWingMesh = True

# BODY GEOMETRY

# Basic dimensions
x0 = -0.5  # Initial position
length = 2.0  # cylinder length
radius = 0.5  # cylinder ratio

# Set cylinder points
nu, nv = 100, 100
pts = zeros((nu, nv, 3))
for i in range(nu):
    for j in range(nv):
        theta = pi * j / (nv - 1)
        pts[i, j, 0] = x0 + length * i / (nu - 1)
        pts[i, j, 1] = -radius * cos(theta)
        pts[i, j, 2] = radius * sin(theta)

# Save body points
bodyPts = pts[:, :, :]

# Create surface object
bodySurf = Surface("body", bodyPts)

bodyDict = {"body": bodySurf}
bodyGeom = Geometry(bodyDict)

# AIRFOIL COORDINATES
# Parameters
n = 81  # Number of points in the airfoil
tc = 0.12  # Airfoil thickness
chord = 1.0  # Airfoil chord
x0 = 0.0  # Leading edge position
y0 = -0.0  # Leading edge position
z0 = radius  # Leading edge position

# Generate non-dimensional coordinates
if trailingEdge == "sharp":
    x = (1 - cos(linspace(0, 1, n) * pi)) / 2
    y = 5 * tc * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)

elif trailingEdge == "blunt":
    x = (1 - cos(linspace(0, 1, n) * pi)) / 2
    y = 5 * tc * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)

    # Estimate spacing between last two nodes
    dy = sqrt((x[-1] - x[-2]) ** 2 + (y[-1] - y[-2]) ** 2)

    # Estimate how many nodes we can fit on the blunt TE (symmetric part only)
    nTEnodes = int(ceil(y[-1] / dy) + 1)

    # Generate coordinates for the vertical segment
    yTE = linspace(y[-1], 0, nTEnodes)
    xTE = ones(nTEnodes)

    # Append the coordinates
    x = hstack([x, xTE[1:]])
    y = hstack([y, yTE[1:]])

# Scale and translate
x = x0 + hstack([x[::-1], x[1:]]) * chord
y = y0 + hstack([-y[::-1], y[1:]]) * chord
z = z0 * ones(len(x))

# Assemble coordinates
rBaseline = vstack([x, y, z]).T

# Project to the cylinder surface
rBaseline, _ = bodyGeom.project_on_surface(rBaseline)

# WING COORDINATES
# We will use the projected airfoil to extreude the wing

# Parameters
span = 1.0
nv = 4  # Number of span-wise control points

# Determine number of chord-wise points
nu = rBaseline.shape[0]

# Determine points
pts = zeros((nu, nv, 3))
for i in range(nu):
    for j in range(nv):
        pts[i, j, 0] = rBaseline[i, 0]
        pts[i, j, 1] = rBaseline[i, 1]
        pts[i, j, 2] = rBaseline[i, 2] + span * (j / (nv - 1) - 0.05)

if trailingEdge == "blunt":

    rBaselineTE = vstack((rBaseline[-nTEnodes:, :], rBaseline[:nTEnodes, :]))

    # Determine number of chord-wise points
    nu = rBaselineTE.shape[0]

    # Determine points
    ptsTE = zeros((nu, nv, 3))
    for i in range(nu):
        for j in range(nv):
            ptsTE[i, j, 0] = rBaselineTE[i, 0]
            ptsTE[i, j, 1] = rBaselineTE[i, 1]
            ptsTE[i, j, 2] = rBaselineTE[i, 2] + span * (j / (nv - 1) - 0.05)

    # Save wing points
    wingPtsTE = ptsTE[:, :, :]

    # Create surface object
    wingSurfTE = Surface("wingTE", wingPtsTE)

    # Save wing points
    wingPts = pts[nTEnodes - 1 : -nTEnodes + 1, :, :]

else:

    wingPts = pts.copy()

# Create surface object
wingSurf = Surface("wing", wingPts)

# Determine points
nu = 6
nv = 4
curve1_pts = zeros((nu, nv, 3))
for i in range(nu):
    for j in range(nv):
        curve1_pts[i, j, 0] = rBaseline[-1, 0]
        curve1_pts[i, j, 1] = rBaseline[-1, 1]
        curve1_pts[i, j, 2] = rBaseline[-1, 2] + span * (i / (nv - 1) - 0.05)

curve1 = Curve("curve1", curve1_pts)

wingDict = {"wing": wingSurf, "curve1": curve1}
if trailingEdge == "blunt":
    wingDict["wingTE"] = wingSurfTE
wingGeom = Geometry(wingDict)

# ==============================================================#

# BODY PROBLEM

# Set options
if trailingEdge == "blunt":
    options = {
        "bc1": "continuous",
        "bc2": "continuous",
        "dStart": 0.02,
        "numLayers": 17,
        "extension": 1.7,
        "epsE0": 10.0,
        "theta": 0.0,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "numAreaPasses": 0,
        "cMax": 100.0,
    }

if trailingEdge == "sharp":
    options = {
        "bc1": "continuous",
        "bc2": "continuous",
        "dStart": 0.02,
        "numLayers": 17,
        "extension": 1.7,
        "epsE0": 1.5,
        "theta": 0.5,
        "alphaP0": 0.25,
        "numSmoothingPasses": 5,
        "numAreaPasses": 5,
        "cMax": 1.0,
    }

# Set meshing problem
bodyMesh = hypsurf.HypSurfMesh(curve=rBaseline, ref_geom=bodyGeom, options=options)

if generateBodyMesh:

    # Print log
    print("Generating body mesh")

    # Generate mesh
    bodyMesh.createMesh()

    # Export mesh
    bodyMesh.exportPlot3d("bodyMesh.xyz")

# WING PROBLEM

# Set options
if trailingEdge == "blunt":
    options = {
        "bc1": "curve.curve1",
        "bc2": "curve.curve1",
        "dStart": 0.002,
        "numLayers": 17,
        "extension": 1.2,
        "epsE0": 2.0,
        "theta": 0.0,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "numAreaPasses": 0,
        "cMax": 1.0,
    }

if trailingEdge == "sharp":
    options = {
        "bc1": "curve.curve1",
        "bc2": "curve.curve1",
        "dStart": 0.002,
        "numLayers": 17,
        "extension": 1.2,
        "epsE0": 1.5,
        "theta": 1.5,
        "alphaP0": 0.25,
        "numSmoothingPasses": 0,
        "numAreaPasses": 0,
        "cMax": 1.0,
    }

# Set meshing problem
wingMesh = hypsurf.HypSurfMesh(curve=rBaseline, ref_geom=wingGeom, options=options)

if generateWingMesh:

    # Print log
    print("Generating wing mesh")

    # Generate mesh
    wingMesh.createMesh()

    # Export mesh
    wingMesh.exportPlot3d("wingMesh.xyz")

# EXPORTING REFERENCE SURFACES

# Replace coordinates of the body mesh
bodyMesh.mesh = zeros((3, bodyPts.shape[0], bodyPts.shape[1]))
bodyMesh.mesh[0, :, :] = bodyPts[:, :, 0]
bodyMesh.mesh[1, :, :] = bodyPts[:, :, 1]
bodyMesh.mesh[2, :, :] = bodyPts[:, :, 2]

# Export body
bodyMesh.exportPlot3d("bodySurf.xyz")

# Replace coordinates of the wing mesh
wingMesh.mesh = zeros((3, wingPts.shape[0], wingPts.shape[1]))
wingMesh.mesh[0, :, :] = wingPts[:, :, 0]
wingMesh.mesh[1, :, :] = wingPts[:, :, 1]
wingMesh.mesh[2, :, :] = wingPts[:, :, 2]

# Export wing
wingMesh.exportPlot3d("wingSurf.xyz")

if trailingEdge == "blunt":
    # Replace coordinates of the wing mesh
    wingMesh.mesh = zeros((3, wingPtsTE.shape[0], wingPtsTE.shape[1]))
    wingMesh.mesh[0, :, :] = wingPtsTE[:, :, 0]
    wingMesh.mesh[1, :, :] = wingPtsTE[:, :, 1]
    wingMesh.mesh[2, :, :] = wingPtsTE[:, :, 2]

    # Export wing
    wingMesh.exportPlot3d("wingSurfTE.xyz")

# Open tecplot
os.system("tec360 layout_wingBody.lay")
