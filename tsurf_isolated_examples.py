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

#example = 'kink_on_plate'
#example = 'line_on_cylinder'
example = 'line_on_cylinder_with_guides'

if example == 'kink_on_plate':

    # Read inputs from CGNS file
    geom = TSurfGeometry('examples/inputs/plate.cgns')

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
    curve = np.zeros((numNodes,3))
    curve[:,0] = 1.0
    curve[:,1] = np.linspace(0.6, 0.1, numNodes)*5
    curve[:,2] = np.linspace(0.1, 0.6, numNodes)*5
    
    # Flip curve to change marching direction
    curve = curve[::-1,:]
    
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
    geom = TSurfGeometry('examples/inputs/cylinder.cgns')
    
    # Flip BC curve
    geom.curves['source'].flip()

    # Set problem
    curve = 'source'
    bc1 = 'splay'
    bc2 = 'curve:bc1'
        
    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 3
    nuArea = 0.16
    numAreaPasses = 3
    sigmaSplay = 0.2
    cMax = 1.0
    ratioGuess = 20
    
    # Options
    sBaseline = 0.01#0.01
    numLayers = 30
    extension = 2.8
    guideCurves = []

elif example == 'line_on_cylinder_with_guides':

    # Read inputs from CGNS file
    geom = TSurfGeometry('examples/inputs/cylinder.cgns')
    
    # Flip BC curve
    geom.curves['source'].flip()

    # Set problem
    numNodes = 50
    curve = np.zeros((numNodes,3))
    curve[:,0] = 1.03
    curve[:,1] = np.linspace(0.6, 0.0, numNodes)
    curve[:,2] = np.linspace(0.0, 0.6, numNodes)
    bc1 = 'splay'
    bc2 = 'splay'
        
    # Set parameters
    epsE0 = 5.5
    theta = 0.0
    alphaP0 = 0.25
    numSmoothingPasses = 0
    nuArea = 0.16
    numAreaPasses = 0
    sigmaSplay = 0.2
    cMax = 1000.0
    ratioGuess = 20
    
    # Options
    sBaseline = 0.01#0.01
    numLayers = 30
    extension = 2.8
    guideCurves = ['bc1']

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

mesh.createMesh()

mesh.exportPlot3d('output.xyz')

hs = hypsurf.hypsurfAPI.hypsurfapi

### DOT PRODUCT TEST FOR SMOOTHING
coor = mesh.mesh[:,:,3].T.flatten()#geom.curves[curve].coor.reshape(-1, order='F')
hs.smoothing(coor, 5., alphaP0, 1, numLayers)

# FORWARD MODE
coord = np.random.random_sample(coor.shape)
coord_copy = coord.copy()
rOut, rOutd = hs.smoothing_d(coor, coord, 5., alphaP0, 1, numLayers)

# REVERSE MODE
rOutb = np.random.random_sample(coor.shape)**2
rOutb_copy = rOutb.copy()
coorb, rOut = hs.smoothing_b(coor, 5., alphaP0, 1, numLayers, rOutb)

print
print 'Dot product test for Smoothing. This should be zero:'
dotprod = np.sum(coord_copy*coorb) - np.sum(rOutd*rOutb_copy)
print dotprod
print
np.testing.assert_almost_equal(dotprod, 0.)

### DOT PRODUCT TEST FOR PROJECTION

r0 = mesh.mesh[:,:,3].T.flatten()

# Clean all projection dictionaries
mesh.projDict = []
mesh.curveProjDict1 = []
mesh.curveProjDict2 = []
for index in range(len(mesh.curveProjDictGuide)):
    mesh.curveProjDictGuide[index] = []

# FORWARD PASS
rNext, NNext = mesh.projection(r0)

# Store projection dictionaries
projDict = mesh.projDict[:]
curveProjDict1 = mesh.curveProjDict1[:]
curveProjDict2 = mesh.curveProjDict2[:]

# FORWARD DERIVATIVES
r0d = np.array(np.random.rand(r0.shape[0]),order='F')
r0d = r0d/np.sqrt(np.sum(r0d**2))
r0d_copy = r0d.copy()

mesh.coord[:,:] = np.array(np.random.rand(mesh.coord.shape[0],mesh.coord.shape[1]),order='F')
mesh.coord[:,:] = mesh.coord/np.sqrt(np.sum(mesh.coord**2))
coord_copy = mesh.coord.copy()

curveCoord_copy = {}
for curveName in mesh.curveCoord:
    mesh.curveCoord[curveName][:,:] = np.array(np.random.rand(mesh.curveCoord[curveName].shape[0],mesh.curveCoord[curveName].shape[1]),order='F')
    mesh.curveCoord[curveName][:,:] = mesh.curveCoord[curveName][:,:]/np.sqrt(np.sum(mesh.curveCoord[curveName]**2))
    curveCoord_copy[curveName] = mesh.curveCoord[curveName].copy()

rNextd, NNextd = mesh.projection_d(r0, r0d, rNext, NNext, 0)
r0d = r0d_copy
mesh.coord = coord_copy
for curveName in mesh.curveCoord:
    mesh.curveCoord[curveName] = curveCoord_copy[curveName]

# REVERSE DERIVATIVES
rNextb = np.random.rand(rNext.shape[0])
NNextb = np.random.rand(NNext.shape[0],NNext.shape[1])
rNextb_copy = rNextb.copy()
NNextb_copy = NNextb.copy()

mesh.coorb[:,:] = 0.0

# Restore projection dictionaries
mesh.projDict = projDict
mesh.curveProjDict1 = curveProjDict1
mesh.curveProjDict2 = curveProjDict2

for curveName in mesh.curveCoord:
    mesh.curveCoorb[curveName][:,:] = 0.0

r0b = mesh.projection_b(r0, rNext, rNextb, NNext, NNextb, 0)
rNextb = rNextb_copy
NNextb = NNextb_copy

dotprod = 0.0
dotprod = dotprod + np.sum(r0d*r0b)
dotprod = dotprod + np.sum(mesh.coord*mesh.coorb)
for curveName in mesh.curveCoord:
    dotprod = dotprod + np.sum(mesh.curveCoord[curveName]*mesh.curveCoorb[curveName])
dotprod = dotprod - np.sum(rNextd*rNextb)
dotprod = dotprod - np.sum(NNextd*NNextb)

print
print 'Dot product test for Projection. This should be zero:'
print dotprod
print

# FORWARD DERIVATIVES WITH FINITE DIFFERENCING
stepSize = 1e-7

r0_FD = r0 + r0d*stepSize
mesh.ref_geom.update(mesh.ref_geom.coor + np.array(mesh.coord*stepSize,order='F'))
for curveName in mesh.curveCoord:
    mesh.ref_geom.curves[curveName].coor = mesh.ref_geom.curves[curveName].coor + mesh.curveCoord[curveName]*stepSize

rNext_FD, NNext_FD = mesh.projection(r0_FD)

rNextd_FD = (rNext_FD-rNext)/stepSize
NNextd_FD = (NNext_FD-NNext)/stepSize

FD_check = np.max([np.max(rNextd - rNextd_FD), np.max(NNextd - NNextd_FD)])

print
print 'Finite difference test for Projection. This should be around 1e-8:'
print FD_check
print


### DOT PRODUCT TEST FOR COMPUTEMATRICES
n = 4
r0 = np.random.rand(3*n) * 10
rm1 = r0 / 2
n0 = np.random.rand(3, n)
s0 = np.random.rand(n) * 4
sm1 = s0 / 2
layerindex = 1

# FORWARD MODE
r0_d = np.random.random_sample(r0.shape)
n0_d = np.random.random_sample(n0.shape)
s0_d = np.random.random_sample(s0.shape)
rm1_d = np.random.random_sample(rm1.shape)
sm1_d = np.random.random_sample(sm1.shape)

r0_d_copy = r0_d.copy()
n0_d_copy = n0_d.copy()
s0_d_copy = s0_d.copy()
rm1_d_copy = rm1_d.copy()
sm1_d_copy = sm1_d.copy()

retainSpacing = True

rnext, rnext_d = hs.computematrices_d(r0, r0_d, n0, n0_d, s0, s0_d, rm1, rm1_d, sm1, sm1_d, layerindex, theta, sigmaSplay, bc1, bc2, mesh.guideIndices, retainSpacing, numLayers, epsE0)

# REVERSE MODE
rnext_b = np.random.random_sample(r0.shape)

rnext_b_copy = rnext_b.copy()

r0_b, n0_b, s0_b, rm1_b, sm1_b, rnext = hs.computematrices_b(r0, n0, s0, rm1, sm1, layerindex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, mesh.guideIndices, retainSpacing, rnext_b)


# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(rnext_d*rnext_b_copy)
dotProd = dotProd - np.sum(r0_b*r0_d_copy)
dotProd = dotProd - np.sum(n0_b*n0_d_copy)
dotProd = dotProd - np.sum(s0_b*s0_d_copy)
dotProd = dotProd - np.sum(rm1_b*rm1_d_copy)
dotProd = dotProd - np.sum(sm1_b*sm1_d_copy)

print ''
print 'Dot product test for ComputeMatrices. This should be zero:'
print dotProd
print ''

### DOT PRODUCT TEST FOR AREAFACTOR

r0 = mesh.mesh[:,:,3].T.flatten()
d = sBaseline
layerindex = 1

# FORWARD MODE
r0_d = np.random.random_sample(r0.shape)
d_d = np.random.random_sample(1)[0]

r0_d_copy = r0_d.copy()
d_d_copy = d_d

s,s_d,maxstretch =  hs.areafactor_test_d(r0, r0_d, d, d_d, nuArea, numAreaPasses, bc1, bc2, mesh.guideIndices)

# REVERSE MODE
s_b = np.random.random_sample(s.shape)

s_b_copy = s_b.copy()

r0_b,d_b = hs.areafactor_test_b(r0, d, nuArea, numAreaPasses, bc1, bc2, mesh.guideIndices, s, s_b, maxstretch)


# Dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(r0_d_copy*r0_b)
dotProd = dotProd + d_d_copy*d_b
dotProd = dotProd - np.sum(s_d*s_b_copy)

print ''
print 'Dot product test for AreaFactor. This should be zero:'
print dotProd
print ''

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

