# IMPORTS

from pysurf import hypsurf
import numpy as np
from numpy.random import rand
from pysurf import TSurfGeometry
from pysurf import tsurf_tools
import os
from mpi4py import MPI
import unittest
import pickle

geom = TSurfGeometry('../../../inputs/rect_uns.cgns',['geom'])

# ADDING GUIDE CURVES
# Create curve dictionary based on imported curves
# !!! Make sure to call `extract_curves.py` before running this script
curves = []
curves.append(tsurf_tools.read_tecplot_curves('extracted_curve_000.plt_'))
curves.append(tsurf_tools.read_tecplot_curves('extracted_curve_001.plt_'))
curves.append(tsurf_tools.read_tecplot_curves('extracted_curve_002.plt_'))
curves.append(tsurf_tools.read_tecplot_curves('extracted_curve_003.plt_'))

# Create an empty list in which we'll store the long (longitudinal) edges of
# the rect
long_curves = []
counter = 0

# Examine each of the imported curves
for ext_curve in curves:

    # Split these curves based on sharpness to get all edges of the rect
    split_curve = tsurf_tools.split_curve_single(ext_curve[ext_curve.keys()[0]], 'int', criteria='sharpness')

    # Loop over these split curves
    for name in split_curve:

        # Give the curves new names so they do not conflict with each other
        split_curve[name].name = 'int_'+'{}'.format(counter).zfill(3)
        counter += 1

        # Extract the curve points and compute the length
        pts = split_curve[name].get_points()
        length = pts[2, 0] - pts[2, -1]

        # Hardcoded logic here based on the length of edges
        if np.abs(length) > 1:

            # Flip the long curve if it's facing the 'wrong' direction
            if length > 0:
                split_curve[name].flip()

            # Add the long curve to the list
            long_curves.append(split_curve[name])

# Create a list of guideCurves based on the extracted long curves
# Note here we need the strings, not the curve objects.
# We use the same loop to add the guide curves to the rectangle object
guideCurves = []
for ext_curve in long_curves:
    print ext_curve.name
    geom.add_curve(ext_curve)
    if ext_curve.name != 'int_011':
        guideCurves.append(ext_curve.name)

# Rotate the main component
geom.rotate(5,2)

# Load intersection curve
intCurveDict = tsurf_tools.read_tecplot_curves('intersection.plt_')
curveName = intCurveDict.keys()[0]
curve = intCurveDict[curveName]
curve.shift_end_nodes(criteria='maxX')
curve.flip()

# Set problem
bc1 = 'curve:int_011'
bc2 = 'curve:int_011'

# Set parameters
epsE0 = 4.5
theta = -0.5
alphaP0 = 0.25
numSmoothingPasses = 0
nuArea = 0.16
numAreaPasses = 20
sigmaSplay = 0.3
cMax = 10000.0
ratioGuess = 1.5

# Options
sBaseline = 0.03#0.01
numLayers = 17
extension = 3.5
guideCurves = guideCurves

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
    'remesh':False,
    'guideCurves':guideCurves
        
}

mesh = hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

mesh.test_all()

mesh.exportPlot3d('output.xyz')


    # def test_line_on_cylinder(self):
    #     example = 'line_on_cylinder'
    #
    #     # Read inputs from CGNS file
    #     self.geom = TSurfGeometry('examples/inputs/cylinder.cgns')
    #
    #     # Flip BC self.curve
    #     self.geom.curves['bc1'].flip()
    #
    #     # Set problem
    #     self.curve = 'source'
    #     self.bc1 = 'splay'
    #     self.bc2 = 'self.curve:bc1'
    #
    #     # Set parameters
    #     self.epsE0 = 5.5
    #     self.theta = 0.0
    #     self.alphaP0 = 0.25
    #     self.numSmoothingPasses = 0
    #     self.nuArea = 0.16
    #     self.numAreaPasses = 0
    #     self.sigmaSplay = 0.2
    #     self.cMax = 1.0
    #     self.ratioGuess = 20
    #
    #     # Options
    #     self.sBaseline = 0.01
    #     self.numLayers = 50
    #     self.extension = 4.8
    #
    #     self.layout_file = 'layout_cylinder.lay'
    #
    #     self.create_mesh()
    #
    #     np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    # def test_cylinder_cap(self):
    #     example = 'cylinder_cap'
    #
    #     # Read inputs from CGNS file
    #     self.geom = TSurfGeometry('examples/inputs/cylinder.cgns')
    #
    #     # Set reference self.curve
    #     n1 = 7
    #     n2 = 21
    #
    #     x = np.hstack([np.linspace(0.5,0,n1), np.zeros(n2-1), np.linspace(0,0.5,n1)[1:]])
    #     y = np.hstack([-0.5*np.ones(n1), np.linspace(-0.5, 0.5, n2)[1:], 0.5*np.ones(n1-1)])
    #     z = np.zeros(n1+n2+n1-2)
    #     self.curve = np.vstack([x, y, z]).T
    #
    #     # Set boundary conditions
    #     self.bc1 = 'splay'
    #     self.bc2 = 'splay'
    #
    #     # Set parameters
    #     self.epsE0 = 2.5
    #     self.theta = 0.0
    #     self.alphaP0 = 0.25
    #     self.numSmoothingPasses = 0
    #     self.nuArea = 0.16
    #     self.numAreaPasses = 0
    #     self.sigmaSplay = 0.2
    #     self.cMax = 1.0
    #     self.ratioGuess = 20
    #
    #     # Options
    #     self.sBaseline = 0.01
    #     self.numLayers = 40
    #     self.extension = 2.5
    #
    #     self.layout_file = 'layout_cylinder.lay'
    #
    #     self.create_mesh()
    #
    #     np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)
    #
    # # def test_line_on_cubeAndCylinder(self):
    # #     example = 'line_on_cubeAndCylinder'
    # #
    # #     # Read inputs from CGNS file
    # #     #self.geom = TSurfGeometry('inputs/cubeAndCylinder.cgns',['cylinder']) # mesh cylinder only
    # #     self.geom = TSurfGeometry('examples/inputs/cubeAndCylinder.cgns',['self.geom']) # mesh cube only
    # #
    # #     # Set problem
    # #     self.curve = 'diag'
    # #     self.bc1 = 'splay'
    # #     self.bc2 = 'splay'
    # #
    # #     # Set parameters
    # #     self.epsE0 = 1.5
    # #     self.theta = 0.0
    # #     self.alphaP0 = 0.25
    # #     self.numSmoothingPasses = 0
    # #     self.nuArea = 0.16
    # #     self.numAreaPasses = 0
    # #     self.sigmaSplay = 0.1
    # #     self.cMax = 1.0
    # #     self.ratioGuess = 20
    # #
    # #     # Options
    # #     self.sBaseline = 0.01
    # #     self.numLayers = 50
    # #     self.extension = 4.8
    # #
    # #     self.layout_file = 'layout_cubeAndCylinder.lay'
    # #
    # #     self.create_mesh()
    # #
    # #     np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)
    #
    # def test_line_on_cube(self):
    #     example = 'line_on_cube'
    #
    #     # Read inputs from CGNS file
    #     self.geom = TSurfGeometry('examples/inputs/cube.cgns', MPI.COMM_WORLD) # mesh cube only
    #
    #     # Set problem
    #     n = 40
    #     nums = np.linspace(0.75, 0.25, n)
    #     x = nums
    #     y = nums #0.2*np.ones(n)
    #     z = np.zeros(n)
    #     self.curve = np.vstack([x, y, z]).T
    #     self.bc1 = 'splay'
    #     self.bc2 = 'splay'
    #
    #     # Set parameters
    #     self.epsE0 = 2.4
    #     self.theta = 0.0
    #     self.alphaP0 = 0.25
    #     self.numSmoothingPasses = 0
    #     self.nuArea = 0.16
    #     self.numAreaPasses = 0
    #     self.sigmaSplay = 0.2
    #     self.cMax = 10.0
    #     self.ratioGuess = 20
    #
    #     # Options
    #     self.sBaseline = 0.01
    #     self.numLayers = 60
    #     self.extension = 5.1
    #
    #     self.layout_file = 'layout_cube.lay'
    #
    #     self.create_mesh()
    #
    #     np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

