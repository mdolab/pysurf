# IMPORTS

from __future__ import division
from pysurf import hypsurf
import numpy as np
from numpy.random import rand
from pysurf import TSurfGeometry
import os
from mpi4py import MPI
import unittest
import pickle

class TestCreateMesh(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestCreateMesh, self).__init__(*args, **kwargs)
        with open('examples_dict.pickle', 'r') as f:
            self.master_dict = pickle.load(f)

    def test_kink_on_plate(self):
        example = 'kink_on_plate'

        # Read inputs from CGNS file
        self.geom = TSurfGeometry('examples/inputs/plate.cgns')

        # Set source self.curve
        self.curve = np.array([[0,0,0],
                       [1,1,0],
                       [2,2,0],
                       [3,2,0],
                       [4,1,0],
                       [5,0,0]])

        # Flip self.curve to change marching direction
        self.curve = self.curve[::-1,:]

        # Define boundary conditions
        self.bc1 = 'constX'
        self.bc2 = 'constX'

        # Set parameters
        self.epsE0 = 1.0
        self.theta = 0.0
        self.alphaP0 = 0.25
        self.numSmoothingPasses = 0
        self.nuArea = 0.16
        self.numAreaPasses = 0
        self.sigmaSplay = 0.5
        self.cMax = 20.0
        self.ratioGuess = 20

        # Options
        self.sBaseline = 0.15
        self.numLayers = 5
        self.extension = 1.5

        # Give layout file
        self.layout_file = 'layout_plate.lay'

        self.create_mesh()

        np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    def test_line_on_cylinder(self):
        example = 'line_on_cylinder'

        # Read inputs from CGNS file
        self.geom = TSurfGeometry('examples/inputs/cylinder.cgns')

        # Flip BC self.curve
        self.geom.curves['bc1'].flip()

        # Set problem
        self.curve = 'source'
        self.bc1 = 'splay'
        self.bc2 = 'self.curve:bc1'

        # Set parameters
        self.epsE0 = 5.5
        self.theta = 0.0
        self.alphaP0 = 0.25
        self.numSmoothingPasses = 0
        self.nuArea = 0.16
        self.numAreaPasses = 0
        self.sigmaSplay = 0.2
        self.cMax = 1.0
        self.ratioGuess = 20

        # Options
        self.sBaseline = 0.01
        self.numLayers = 50
        self.extension = 4.8

        self.layout_file = 'layout_cylinder.lay'

        self.create_mesh()

        np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    def test_cylinder_cap(self):
        example = 'cylinder_cap'

        # Read inputs from CGNS file
        self.geom = TSurfGeometry('examples/inputs/cylinder.cgns')

        # Set reference self.curve
        n1 = 7
        n2 = 21

        x = np.hstack([np.linspace(0.5,0,n1), np.zeros(n2-1), np.linspace(0,0.5,n1)[1:]])
        y = np.hstack([-0.5*np.ones(n1), np.linspace(-0.5, 0.5, n2)[1:], 0.5*np.ones(n1-1)])
        z = np.zeros(n1+n2+n1-2)
        self.curve = np.vstack([x, y, z]).T

        # Set boundary conditions
        self.bc1 = 'splay'
        self.bc2 = 'splay'

        # Set parameters
        self.epsE0 = 2.5
        self.theta = 0.0
        self.alphaP0 = 0.25
        self.numSmoothingPasses = 0
        self.nuArea = 0.16
        self.numAreaPasses = 0
        self.sigmaSplay = 0.2
        self.cMax = 1.0
        self.ratioGuess = 20

        # Options
        self.sBaseline = 0.01
        self.numLayers = 40
        self.extension = 2.5

        self.layout_file = 'layout_cylinder.lay'

        self.create_mesh()

        np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    # def test_line_on_cubeAndCylinder(self):
    #     example = 'line_on_cubeAndCylinder'
    #
    #     # Read inputs from CGNS file
    #     #self.geom = TSurfGeometry('inputs/cubeAndCylinder.cgns',['cylinder']) # mesh cylinder only
    #     self.geom = TSurfGeometry('examples/inputs/cubeAndCylinder.cgns',['self.geom']) # mesh cube only
    #
    #     # Set problem
    #     self.curve = 'diag'
    #     self.bc1 = 'splay'
    #     self.bc2 = 'splay'
    #
    #     # Set parameters
    #     self.epsE0 = 1.5
    #     self.theta = 0.0
    #     self.alphaP0 = 0.25
    #     self.numSmoothingPasses = 0
    #     self.nuArea = 0.16
    #     self.numAreaPasses = 0
    #     self.sigmaSplay = 0.1
    #     self.cMax = 1.0
    #     self.ratioGuess = 20
    #
    #     # Options
    #     self.sBaseline = 0.01
    #     self.numLayers = 50
    #     self.extension = 4.8
    #
    #     self.layout_file = 'layout_cubeAndCylinder.lay'
    #
    #     self.create_mesh()
        #
        # np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    def test_line_on_cube(self):
        example = 'line_on_cube'

        # Read inputs from CGNS file
        self.geom = TSurfGeometry('examples/inputs/cube.cgns', MPI.COMM_WORLD) # mesh cube only

        # Set problem
        n = 40
        nums = np.linspace(0.75, 0.25, n)
        x = nums
        y = nums #0.2*np.ones(n)
        z = np.zeros(n)
        self.curve = np.vstack([x, y, z]).T
        self.bc1 = 'splay'
        self.bc2 = 'splay'

        # Set parameters
        self.epsE0 = 2.4
        self.theta = 0.0
        self.alphaP0 = 0.25
        self.numSmoothingPasses = 0
        self.nuArea = 0.16
        self.numAreaPasses = 0
        self.sigmaSplay = 0.2
        self.cMax = 10.0
        self.ratioGuess = 20

        # Options
        self.sBaseline = 0.01
        self.numLayers = 60
        self.extension = 5.1

        self.layout_file = 'layout_cube.lay'

        self.create_mesh()

        np.testing.assert_almost_equal(self.master_dict[example], self.mesh.mesh)

    def create_mesh(self):

        # Set options
        options = {

            'bc1' : self.bc1,
            'bc2' : self.bc2,
            'dStart' : self.sBaseline,
            'numLayers' : self.numLayers,
            'extension' : self.extension,
            'epsE0' : self.epsE0,
            'theta' : self.theta,
            'alphaP0' : self.alphaP0,
            'numSmoothingPasses' : self.numSmoothingPasses,
            'nuArea' : self.nuArea,
            'numAreaPasses' : self.numAreaPasses,
            'sigmaSplay' : self.sigmaSplay,
            'cMax' : self.cMax,
            'ratioGuess' : self.ratioGuess,

            }

        mesh = hypsurf.HypSurfMesh(curve=self.curve, ref_geom=self.geom, options=options)

        mesh.createMesh()

        # mesh.exportPlot3d('output.xyz')

        # EXPORT REFERENCE SURFACE

        # Open tecplot
        # os.system('tec360 ' + self.layout_file)

        self.mesh = mesh



if __name__ == "__main__":
    unittest.main()
