# IMPORTS

from __future__ import division
from pysurf import hypsurf
import numpy as np
from numpy.random import rand
from pysurf import TSurfGeometry
from pysurf import tsurf_tools
import os
from mpi4py import MPI
import unittest
import pickle

#example = 'kink_on_plate'
#example = 'line_on_cylinder'
#example = 'line_on_cylinder_with_guides'
example = 'crm_wing'

if example == 'kink_on_plate':

    # Read inputs from CGNS file
    geom = TSurfGeometry('../inputs/plate.cgns')

    # Set source curve
    '''
    curve = np.array([[0.,0,0],
                      [1,1,0],
                      [2,2,0],
                      [3,2,0],
                      [4,1,0],
                      [5,0,0]])
    '''
    numNodes = 70
    coor = np.zeros((numNodes,3))
    coor[:,0] = 1.0
    coor[:,1] = np.linspace(0.6, 0.1, numNodes)*5
    coor[:,2] = np.linspace(0.1, 0.6, numNodes)*5
    
    # Create curve object
    curve = tsurf_tools.create_curve_from_points(coor, 'source')

    # Flip curve to change marching direction
    curve.flip()
    
    # Define boundary conditions
    bc1 = 'splay'
    bc2 = 'splay'
    
    # Set parameters
    epsE0 = 3.0
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 5
    nuArea = 0.16
    numAreaPasses = 5
    sigmaSplay = 0.05
    cMax = 1.0
    ratioGuess = 5
    
    # Options
    sBaseline = 0.005
    numLayers = 75
    extension = 6.0
    guideCurves = []

elif example == 'line_on_cylinder':

    # Read inputs from CGNS file
    geom = TSurfGeometry('../inputs/cylinder.cgns')
    geom.curves['bc1'].translate(-0.00001,0,0)    

    # Set initial curve
    curve = geom.curves['source']
    curve.flip()
    curve = curve.remesh(nNewNodes=15)

    # Set problem
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
    numLayers = 30
    extension = 2.8
    guideCurves = []

elif example == 'line_on_cylinder_with_guides':

    # Read inputs from CGNS file
    geom = TSurfGeometry('../inputs/cylinder.cgns')
    
    # Set initial curve
    numNodes = 50
    coor = np.zeros((numNodes,3))
    coor[:,0] = 1.03
    coor[:,1] = np.linspace(0.6, 0.0, numNodes)
    coor[:,2] = np.linspace(0.0, 0.6, numNodes)
    curve = tsurf_tools.create_curve_from_points(coor, 'seed')

    # Set problem
    bc1 = 'splay'
    bc2 = 'splay'
        
    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 3
    nuArea = 0.16
    numAreaPasses = 5
    sigmaSplay = 0.2
    cMax = 1.0
    ratioGuess = 20
    
    # Options
    sBaseline = 0.01#0.01
    numLayers = 30
    extension = 2.8
    guideCurves = ['bc1']

elif example == 'crm_wing':

    geom = TSurfGeometry('../inputs/initial_full_wing_crm4.cgns',['wing'])

    # Load TE curves and append them to the wing component
    curve_te_upp = tsurf_tools.read_tecplot_curves('../manager_tests/crm/curve_te_upp.plt_')
    curve_te_low = tsurf_tools.read_tecplot_curves('../manager_tests/crm/curve_te_low.plt_')
    curveName = curve_te_upp.keys()[0]
    geom.add_curve(curve_te_upp[curveName])
    curveName = curve_te_low.keys()[0]
    geom.add_curve(curve_te_low[curveName])
    geom.curves[curveName].flip()

    # Load intersection curve
    intCurveDict = tsurf_tools.read_tecplot_curves('../manager_tests/crm/intersection.plt_')
    curveName = intCurveDict.keys()[0]
    curve = intCurveDict[curveName]
    curve.shift_end_nodes(criteria='maxX')

    # Set problem
    bc1 = 'curve:curve_te_upp'
    bc2 = 'curve:curve_te_upp'
        
    # Set parameters
    epsE0 = 12.5
    theta = -0.8
    alphaP0 = 0.25
    numSmoothingPasses = 4
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.3
    cMax = 10.0
    ratioGuess = 10.0
    
    # Options
    sBaseline = 0.01#0.01
    numLayers = 49
    extension = 1.4
    guideCurves = ['curve_te_low']

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
    'remesh':True,
    'guideCurves':guideCurves
        
}

mesh = hypsurf.HypSurfMesh(curve=curve, ref_geom=geom, options=options)

mesh.test_all()

mesh.export_plot3d('output.xyz')


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

