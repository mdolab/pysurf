# 3D Hyperbolic Surface Grid generator
# Written by:
# Ney Rafael Secco
# John Jasa
# neysecco@umich.edu
# July 2016

from __future__ import division
from time import time
import numpy as np
import pdb
import hypsurfAPI
import pysurf
import copy

'''
TO DO

- set up scaling factor based on distance
- blend the angle-based dissipation coefficient
'''

class HypSurfMesh(object):
    """
    Hyperbolic surface mesh created by projecting points extruded from a
    curve onto a surface.

    """

    def __init__(self, curve, ref_geom, options={}, meshName='mesh'):

        '''
        This initializes a hyperbolic surface mesh

        INPUTS

        curve: Curve object used to march the mesh from.

        ref_geom: This is a pysurf geometry object. It's surface will be used
        as reference for the mesh growth. The curves contained in this object could
        be used as source curves, guide curve, or boundary conditions
        '''

        # Store the mesh name
        self.name = meshName

        # Flag that we still do not know the user provided marching parameter. It could be
        # either extension or growth ratio
        self.extension_given = None

        # Apply options
        self._getDefaultOptions()
        self._applyUserOptions(options)

        # Store the curve object given by the user
        self.curve = curve

        # Extract curve coordinates
        rStart = curve.get_points().T

        # Detect nodes that should follow guide curves
        guideIndices = []

        if self.optionsDict['guideCurves']:
            for curve in self.optionsDict['guideCurves']:
                guideCurve = ref_geom.curves[curve]
                guideIndices.append(closest_node(guideCurve, self.curve))

        newGuideIndices = np.array(sorted(guideIndices))
        origGuideCurves = copy.copy(self.optionsDict['guideCurves'])
        for i, index in enumerate(newGuideIndices):
            j = np.where(index==guideIndices)[0][0]
            self.optionsDict['guideCurves'][i] = origGuideCurves[j]
        self.guideIndices = newGuideIndices

        if 0 in self.guideIndices:
            raise NameError('the guide curve is near the first node of the seed curve. Use it as a curve boundary condition instead.')

        # Store the reference geometry object
        self.ref_geom = ref_geom

        # Store curve names if we have curve BCs
        # we assume that curve BCs will be defined as: 'curve:<curve_name>'
        if self.optionsDict['bc1'].lower().startswith('curve'):
            self.ref_curve1 = self.optionsDict['bc1'][6:]
        else:
            self.ref_curve1 = []

        if self.optionsDict['bc2'].lower().startswith('curve'):
            self.ref_curve2 = self.optionsDict['bc2'][6:]
        else:
            self.ref_curve2 = []

        # Get the number of nodes and layers
        self.numNodes = rStart.shape[0]
        self.numLayers = self.optionsDict['numLayers']
        self.mesh = np.zeros((3, self.numNodes, self.numLayers))
        self.meshd = np.zeros((3, self.numNodes, self.numLayers))
        self.meshb = np.zeros((3, self.numNodes, self.numLayers))

        # Initialize list to gather dictionary with projection information
        # that will be used to compute derivatives

        self.projDict = []
        self.curveProjDict1 = []
        self.curveProjDict2 = []

        # Create list of dictionaries to store projection information for each guide curve
        self.curveProjDictGuide = [[] for s in range(len(guideIndices))]

    def createMesh(self, fortran_flag=True):
        '''
        This is the main function that generates the surface mesh from an initial curve

        INPUTS
        rStart: 2D array (numNodes x 3) -> coordinates of the initial curve
                rStart = [x1 y1 z1,
                          x2 y2 z2,
                          x3 y3 z3, ...]

        bc1: string -> boundary condition applied at the first node
             options: 'free', 'splay', 'constX', 'constY', 'constZ'

        bc2: string -> boundary condition applied at the last node
             options: 'free', 'splay', 'constX', 'constY', 'constZ'

        dStart: float -> initial layer spacing to be used when generating the first layer.
                This value will grow geometrically until the last layer

        numLayers: integer -> total number of layers (including the starting one) that
                   should be in the final mesh

        extension: float -> total marching distance, normalized by the largest size of
                   the bounding box surrounding the initial curve

        surf: surface object -> surface object that should contain inverse_evaluate,
              get_points, and get_normals methods
        '''

        from time import time
        st = time()


        # Make sure variables have the correct type
        dStart = float(self.optionsDict['dStart'])
        cMax = float(self.optionsDict['cMax'])
        ratioGuess = float(self.optionsDict['ratioGuess'])
        theta = float(self.optionsDict['theta'])
        bc1 = self.optionsDict['bc1']
        bc2 = self.optionsDict['bc2']
        sigmaSplay = float(self.optionsDict['sigmaSplay'])
        numLayers = int(self.optionsDict['numLayers'])
        epsE0 = float(self.optionsDict['epsE0'])
        theta = float(self.optionsDict['theta'])
        alphaP0 = float(self.optionsDict['alphaP0'])
        numSmoothingPasses = int(self.optionsDict['numSmoothingPasses'])
        numAreaPasses = int(self.optionsDict['numAreaPasses'])
        nuArea = float(self.optionsDict['nuArea'])

        # Get coordinates of the initial curve (this is a flattened array)
        rStart = self.curve.get_points().T
        rStart = rStart.flatten().astype(float)
        rStart = np.array(rStart,order='F')

        # Clean the projection dictionaries
        self.projDict = []
        self.curveProjDict1 = []
        self.curveProjDict2 = []

        if fortran_flag:

            # Perform the marching algorithm and output the results into R,
            # which contains the mesh.
            # fail is a flag set to true if the marching algo failed
            # ratios is the ratios of quality for the mesh

            R, fail, ratios, majorIndices = hypsurfAPI.hypsurfapi.march(self.projection, rStart, dStart, theta, sigmaSplay, bc1.lower(), bc2.lower(), epsE0, alphaP0, self.marchParameter, nuArea, ratioGuess, cMax, self.extension_given, self.guideIndices+1, self.optionsDict['remesh'], numSmoothingPasses, numAreaPasses, numLayers)

            # Obtain the pseudomesh, or subiterations mesh from the three stages of marching.
            # These are used in the adjoint formulation.
            self.R_initial_march = np.array(hypsurfAPI.hypsurfapi.r_initial_march,order='F')
            self.R_smoothed = np.array(hypsurfAPI.hypsurfapi.r_smoothed,order='F')
            self.R_projected = np.array(hypsurfAPI.hypsurfapi.r_projected,order='F')
            self.R_remeshed = np.array(hypsurfAPI.hypsurfapi.r_remeshed,order='F')
            self.R_final = np.array(hypsurfAPI.hypsurfapi.r_final,order='F')
            self.S0_hist = np.array(hypsurfAPI.hypsurfapi.s0_hist,order='F')
            self.Sm1_hist = np.array(hypsurfAPI.hypsurfapi.sm1_hist,order='F')
            self.N_projected = np.array(hypsurfAPI.hypsurfapi.n_projected,order='F')
            self.N_final = np.array(hypsurfAPI.hypsurfapi.n_final,order='F')

            # Release memory for next runs
            hypsurfAPI.hypsurfapi.releasememory()

        else: # We will use the Python version of hypsurf

            # Run the Python marching code
            import hypsurf_python
            R, fail, ratios = hypsurf_python.march(self.projection, rStart, dStart, theta, sigmaSplay, bc1, bc2, epsE0, alphaP0, self.marchParameter, nuArea, ratioGuess, cMax, self.extension_given, self.guideIndices, self.optionsDict['remesh'], numSmoothingPasses, numAreaPasses, numLayers)

        #=====================================
        print time() - st, 'secs'

        if self.optionsDict['plotQuality']:
            fail, ratios = self.qualityCheck(R)
            view_mat(ratios)

        # Convert to X, Y and Z
        X = R[:,::3]
        Y = R[:,1::3]
        Z = R[:,2::3]

        self.mesh = np.zeros((3, int(self.numNodes), int(X.T.shape[1])))
        self.mesh[0, :, :] = X.T
        self.mesh[1, :, :] = Y.T
        self.mesh[2, :, :] = Z.T

    def clean(self):

        '''
        This method removes the residual data generated by the original
        and differentiated code.
        '''

        # Clean projection dictionaries
        self.projDict = []
        self.curveProjDict1 = []
        self.curveProjDict2 = []
        for index in range(len(self.curveProjDictGuide)):
            self.curveProjDictGuide[index] = []

        # Clean reverse derivatives
        self.ref_geom.clean_reverseADSeeds()

    #================================================================
    #================================================================
    # FORWARD AD METHODS

    # The methods are defined in the order that they should be used.

    def set_forwardADSeeds(self, meshd):

        '''
        This function will just overwrite the derivative seeds
        used by the forward mode AD.

        ATTENTION: The user should set derivatives of the underlying triangulated surface
        by using its own methods, such as self.ref_geom.set_forwardADSeeds(), or
        self.ref_geom.set_randomADSeeds().

        INPUTS:

        meshd: float[3,numNodes,numLayers] -> Forward AD seeds of the surface mesh.
        '''

        # Verify shape
        if self.meshd.shape != meshd.shape:
            raise ValueError('The shape of derivatives array in not consistent.')

        self.meshd = np.array(meshd, order='F')

    def compute_forwardAD(self,updateSeeds=True):

        '''
        This method uses forward mode AD to propagate derivatives from inputs to outputs.

        This method updates self.meshd. Then user can use self.get_forwardAD_outputSeeds
        to retrieve this result.

        updateSeeds: boolean -> If True, the code will get derivative seeds from the curve
        object. Otherwise, it will use whatever is stored in self.rStartd.

        ATTENTION:
        The user should call the original method (self.createMesh) first, as this will populate
        the projection dictionaries (self.projDict, self.curveProjDict1, self.curveProjDict2,
        and self.curveProjDictGuides) with the necessary information for the differentiation.
        The user should may call self.set_forwardAD_seeds to define the derivative seeds of
        the input variables, but then it should set updateSeeds to False. The best way is to
        set the derivative seeds directly into the curve object using curve.set_forwardADSeeds.
        '''

        # Check if we have projection dictionaries stored
        if not self.projDict:
            print ''
            print 'ERROR: hypsurf.py - compute_forwardAD'
            print ' Cannot compute derivatives without running the original code first.'
            print ' Call mesh.createMesh, Then use mesh.compute_forwardAD'
            print ''
            exit()

        # Get coordinates of the initial curve (this is a flattened array)
        rStart = self.curve.get_points().T
        rStart = rStart.flatten().astype(float)
        rStart = np.array(rStart,order='F')

        # Get derivative seeds from the starting curve
        rStartd = self.curve.get_forwardADSeeds().T
        rStartd = rStartd.flatten().astype(float)
        rStartd = np.array(rStartd,order='F')

        # Call the Fortran function that computes derivatives
        R, Rd, fail, ratios = hypsurfAPI.hypsurfapi.march_d(self.projection_d,
                                                            rStart,
                                                            rStartd,
                                                            self.R_projected,
                                                            self.R_final,
                                                            self.N_projected,
                                                            self.N_final,
                                                            self.optionsDict['dStart'],
                                                            self.optionsDict['theta'],
                                                            self.optionsDict['sigmaSplay'],
                                                            self.optionsDict['bc1'],
                                                            self.optionsDict['bc2'],
                                                            self.optionsDict['epsE0'],
                                                            self.optionsDict['alphaP0'],
                                                            self.marchParameter,
                                                            self.optionsDict['nuArea'],
                                                            self.optionsDict['ratioGuess'],
                                                            self.optionsDict['cMax'],
                                                            self.guideIndices+1,
                                                            self.optionsDict['remesh'],
                                                            self.extension_given,
                                                            self.optionsDict['numSmoothingPasses'],
                                                            self.optionsDict['numAreaPasses'],
                                                            self.optionsDict['numLayers'])

        # Convert to Xd, Yd and Zd
        Xd = Rd[:,::3]
        Yd = Rd[:,1::3]
        Zd = Rd[:,2::3]

        # Set the forward AD output seeds
        meshd = np.zeros((3, self.numNodes, Xd.T.shape[1]))
        meshd[0, :, :] = Xd.T
        meshd[1, :, :] = Yd.T
        meshd[2, :, :] = Zd.T

        self.set_forwardADSeeds(meshd)

    def get_forwardADSeeds(self):

        # The only output are the derivatives of the surface mesh coordinates
        return self.meshd[:,:,:]

    #================================================================
    #================================================================
    # REVERSE AD METHODS

    # The methods are defined in the order that they should be used.

    def set_reverseADSeeds(self, meshb):

        '''
        This function will just overwrite the surface mesh seeds used by
        the reverse mode AD to propagate derivatives.

        This method updates self.meshb. Then user can use self.get_reversedADSeeds
        to retrieve this result.

        INPUTS:

        meshb -> float[3, numNodes, numLayers] : Reverse derivative seeds
        of the surface mesh coordinates.
        '''

        # Verify shape
        if self.meshb.shape != meshb.shape:
            raise ValueError('The shape of derivatives array in not consistent.')

        # Set values
        self.meshb = np.array(meshb,order='F')

    def compute_reverseAD(self):

        '''
        This method will use reverse mode AD to propagate derivatives from outputs to inputs.

        ATTENTION:
        The user should call the original method (self.createMesh) first, as this will populate
        the projection dictionaries (self.projDict, self.curveProjDict1, self.curveProjDict2,
        and self.curveProjDictGuides) with the necessary information for the differentiation.
        The user should also call self.set_reverseAD_outputSeeds to define the derivative seeds of
        the output variables.

        '''

        # Check if we have projection dictionaries stored
        if not self.projDict:
            print ''
            print 'ERROR: hypsurf.py - compute_reverseAD'
            print ' Cannot compute derivatives without running the original code first.'
            print ' 1- Call mesh.createMesh'
            print ' 2- Set seeds with mesh.set_reverseAD_outputSeeds'
            print ' 3- Then use mesh.compute_reverseAD'
            print ''
            exit()

        # Get coordinates of the initial curve (this is a flattened array)
        rStart = self.curve.get_points().T
        rStart = rStart.flatten().astype(float)
        rStart = np.array(rStart,order='F')

        # Get total number of projection steps
        numProjs = len(self.projDict)

        # Transform the set of 3D coordinates into a 2D array
        X = self.mesh[0, :, :].T
        Y = self.mesh[1, :, :].T
        Z = self.mesh[2, :, :].T

        R = np.zeros((self.optionsDict['numLayers'], 3*self.numNodes), order='F')
        R[:,::3] = X
        R[:,1::3] = Y
        R[:,2::3] = Z

        del X, Y, Z

        # Transform the set of 3D derivative seeds into a 2D array
        Xb = self.meshb[0, :, :].T
        Yb = self.meshb[1, :, :].T
        Zb = self.meshb[2, :, :].T

        Rb = np.zeros((self.optionsDict['numLayers'], 3*self.numNodes), order='F')
        Rb[:,::3] = Xb
        Rb[:,1::3] = Yb
        Rb[:,2::3] = Zb

        del Xb, Yb, Zb

        # Call the Fortran function that computes derivatives
        rStartb = hypsurfAPI.hypsurfapi.march_b(self.projection_b,
                                                rStart,
                                                self.R_initial_march,
                                                self.R_smoothed,
                                                self.R_projected,
                                                self.R_remeshed,
                                                self.R_final,
                                                self.N_projected,
                                                self.N_final,
                                                self.Sm1_hist,
                                                self.S0_hist,
                                                self.optionsDict['dStart'],
                                                self.optionsDict['theta'],
                                                self.optionsDict['sigmaSplay'],
                                                self.optionsDict['bc1'],
                                                self.optionsDict['bc2'],
                                                self.optionsDict['epsE0'],
                                                self.optionsDict['alphaP0'],
                                                self.marchParameter,
                                                self.optionsDict['nuArea'],
                                                self.optionsDict['ratioGuess'],
                                                self.optionsDict['cMax'],
                                                self.guideIndices+1,
                                                self.optionsDict['remesh'],
                                                self.extension_given,
                                                self.optionsDict['numSmoothingPasses'],
                                                self.optionsDict['numAreaPasses'],
                                                numProjs,
                                                R,
                                                Rb)

        # Reshape rStartb to a 3xnNodes shape
        rStartb = rStartb.reshape((-1,3)).T

        # Accumulate seeds back to the curve object
        self.curve.accumulate_reverseADSeeds(rStartb)
 
    def get_reverseADSeeds(self, clean=True):

        '''
        This will return the reverse AD seeds associated with the surface mesh
        '''

        meshb = self.meshb[:,:,:]

        if clean:
            self.meshb[:,:,:] = 0.0

        return meshb
       
    def set_randomADSeeds(self, mode='both', fixedSeed=True):

        '''
        This will set random normalized seeds to all variables.
        This can be used for testing purposes.

        mode: ['both','forward','reverse'] -> Which mode should have
        its derivatives replaced.
        '''

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # Set forward AD seeds
        if mode=='forward' or mode=='both':

            meshd = np.array(np.random.rand(self.mesh.shape[0],self.mesh.shape[1],self.mesh.shape[2]),order='F')
            meshd = meshd/np.sqrt(np.sum(meshd**2))
            self.meshd = meshd

        # Set reverse AD seeds
        if mode=='reverse' or mode=='both':

            # Set reverse AD seeds
            meshb = np.array(np.random.rand(self.mesh.shape[0],self.mesh.shape[1],self.mesh.shape[2]),order='F')
            meshb = meshb/np.sqrt(np.sum(meshb**2))
            self.meshb = meshb

        # Return generated seeds

        if mode == 'forward':
            return meshd

        elif mode == 'reverse':
            return meshb

        elif mode == 'both':
            return meshd, meshb

    #================================================================
    #================================================================
    # PROJECTION METHODS

    def projection(self, r):

        '''
        This function will project the nodes defined in r onto the reference geometry
        used by the current mesh object. This function may also store the projection
        dictionaries that should be used in derivative calculation. To avoid
        storing any dictionaries, use storeDict = 0.
        '''

        numNodes = int(r.shape[0] / 3)

        # Initialize normals array
        NNext = np.zeros((3, numNodes))

        # Save endpoints
        node1 = r[:3]
        node2 = r[-3:]

        # Project onto surface and compute surface normals
        rNext, NNext, projDict = self.ref_geom.project_on_surface(r.reshape((numNodes, 3)))
        rNext = rNext.flatten()
        NNext = NNext.T

        # Store projection dictionaries in their corresponding lists
        self.projDict = self.projDict + [projDict]

        # Replace end points if we use curve BCs
        if self.optionsDict['bc1'].lower().startswith('curve'):
            rNext[:3], NNextAux, curveProjDict1 = self.ref_geom.project_on_curve(node1.reshape((1, 3)), curveCandidates=[self.ref_curve1])
            NNext[:, 0] = NNextAux.T[:, 0]
            self.curveProjDict1 = self.curveProjDict1 + [curveProjDict1]

        if self.optionsDict['bc2'].lower().startswith('curve'):
            rNext[-3:], NNextAux, curveProjDict2 = self.ref_geom.project_on_curve(node2.reshape((1, 3)), curveCandidates=[self.ref_curve2])
            NNext[:, -1] = NNextAux.T[:, 0]
            self.curveProjDict2 = self.curveProjDict2 + [curveProjDict2]

        if np.any(self.guideIndices):
            for i, index in enumerate(self.guideIndices):
                curve = self.optionsDict['guideCurves'][i]
                node = r[3*index:3*index+3].reshape((1, 3))
                rNext[3*index:3*index+3], NNextAux, curveProjDict = self.ref_geom.project_on_curve(node, curveCandidates=[curve])
                NNext[:, index] = NNextAux.T[:, 0]
                self.curveProjDictGuide[i] = self.curveProjDictGuide[i] + [curveProjDict]

        return rNext, NNext

    def projection_d(self, r, rd, rNext, NNext, layerID):

        # Save endpoints
        node1 = r[:3].reshape((1, 3))
        node2 = r[-3:].reshape((1, 3))
        node1d = rd[:3].reshape((1, 3))
        node2d = rd[-3:].reshape((1, 3))

        # Pop the first projection dictionary from the list (since we are propagating derivatives forward)
        projDict = self.projDict[layerID]

        # Project onto surface and compute surface normals
        rNextd, NNextd = self.ref_geom.project_on_surface_d(r.reshape((self.numNodes, 3)),
                                                            rd.reshape((self.numNodes, 3)),
                                                            rNext.reshape((self.numNodes, 3)),
                                                            NNext.T,
                                                            projDict)
        rNextd = rNextd.flatten()
        NNextd = NNextd.T

        # Replace end points if we use curve BCs
        if self.optionsDict['bc1'].lower().startswith('curve'):
            curveProjDict1 = self.curveProjDict1[layerID]
            rNextd[:3], NNextAuxd = self.ref_geom.project_on_curve_d(node1, node1d,
                                                                     rNext[:3], NNext[:,0],
                                                                     curveProjDict1)
            NNextd[:, 0] = NNextAuxd.T[:, 0]

        if self.optionsDict['bc2'].lower().startswith('curve'):
            curveProjDict2 = self.curveProjDict2[layerID]
            rNextd[-3:], NNextAuxd = self.ref_geom.project_on_curve_d(node2, node2d,
                                                                      rNext[-3:], NNext[:,-1],
                                                                      curveProjDict2)
            NNextd[:, -1] = NNextAuxd.T[:, 0]

        if self.guideIndices.any():
            for i, index in enumerate(self.guideIndices):

                node = r[3*index:3*index+3].reshape((1, 3))
                noded = rd[3*index:3*index+3].reshape((1, 3))
                curveProjDict = self.curveProjDictGuide[i][layerID]

                rNextd[3*index:3*index+3], NNextAuxd = self.ref_geom.project_on_curve_d(node, noded,
                                                                                        rNext[3*index:3*index+3], NNext[:,index],
                                                                                        curveProjDict)
                NNextd[:, index] = NNextAuxd.T[:, 0]

        return rNextd, NNextd

    def projection_b(self, r, rNext, rNextb, NNext, NNextb, layerID):

        # Save endpoints
        node1 = r[:3].reshape((1, 3))
        node2 = r[-3:].reshape((1, 3))

        # Pop the last projection dictionary from the list (since we are propagating derivatives forward)
        layerID = int(layerID) # For some reason, f2py does not understand that this should be int
        projDict = self.projDict[layerID]

        # We need to turn off the derivative seeds of the nodes that follow guide curves or curve bcs
        # since they will be replaced later on.
        rNextb_filtered = np.array(rNextb, order='F')
        NNextb_filtered = np.array(NNextb, order='F')

        if self.optionsDict['bc1'].lower().startswith('curve'):
            rNextb_filtered[:3] = 0.0
            NNextb_filtered[:,0] = 0.0

        if self.optionsDict['bc2'].lower().startswith('curve'):
            rNextb_filtered[-3:] = 0.0
            NNextb_filtered[:,-1] = 0.0

        if self.guideIndices.any():
            for index in self.guideIndices:
                rNextb_filtered[3*index:3*index+3] = 0.0
                NNextb_filtered[:,index] = 0.0

        # Project onto surface and compute surface normals
        rb = self.ref_geom.project_on_surface_b(r.reshape((self.numNodes, 3)),
                                                rNext.reshape((self.numNodes, 3)),
                                                rNextb_filtered.reshape((self.numNodes, 3)),
                                                NNext.T,
                                                NNextb_filtered.T,
                                                projDict)

        # Reshape nodal derivatives
        rb = np.array(rb.flatten(),order='F')

        # Replace end points if we use curve BCs
        if self.optionsDict['bc1'].lower().startswith('curve'):
            curveProjDict1 = self.curveProjDict1[layerID]
            rb[:3] = self.ref_geom.project_on_curve_b(node1, rNext[:3], rNextb[:3],
                                                      NNext[:,0], NNextb[:,0],
                                                      curveProjDict1)

        if self.optionsDict['bc2'].lower().startswith('curve'):
            curveProjDict2 = self.curveProjDict2[layerID]
            rb[-3:] = self.ref_geom.project_on_curve_b(node2, rNext[-3:], rNextb[-3:],
                                                       NNext[:,-1], NNextb[:,-1],
                                                       curveProjDict2)

        if self.guideIndices.any():
            for i, index in enumerate(self.guideIndices):
                curveProjDict = self.curveProjDictGuide[i][layerID]
                node = r[3*index:3*index+3].reshape((1, 3))
                curve = self.optionsDict['guideCurves'][i]


                rbAux = self.ref_geom.project_on_curve_b(node,
                                                         rNext[3*index:3*index+3], rNextb[3*index:3*index+3],
                                                         NNext[:,index], NNextb[:,index],
                                                         curveProjDict)
                '''
                elemIDs = curveProjDict['elemIDs']
                curveMask = [1]
                nodeb = np.zeros(node.shape)
                self.ref_geom.curves[curve].project_b(node, nodeb, rNext[3*index:3*index+3], rNextb[3*index:3*index+3],
                                                      NNext[:,index], NNextb[:,index], elemIDs, curveMask)
                rbAux = nodeb
                '''
                rb[3*index:3*index+3] = rbAux

        return rb

    #================================================================
    #================================================================
    # EXPORTING METHODS

    def exportPlot3d(self, filename):

        '''
        This function exports a 3D mesh in plot3d format.
        The user specifies the span and number of nodes in the z direction.
        '''

        # IMPORTS
        from pysurf import plot3d_interface

        # Get grid size
        gridShape = self.mesh.shape[1:]

        # Get coordinates
        X = self.mesh[0,:,:].T
        Y = self.mesh[1,:,:].T
        Z = self.mesh[2,:,:].T

        # Initialize grid object
        myGrid = plot3d_interface.Grid()

        # Expand the coordinate matrices
        X3d = np.zeros((gridShape[1], gridShape[0], 1))
        Y3d = np.zeros((gridShape[1], gridShape[0], 1))
        Z3d = np.zeros((gridShape[1], gridShape[0], 1))
        X3d[:,:,0] = X
        Y3d[:,:,0] = Y
        Z3d[:,:,0] = Z

        # Add mesh to the grid
        myGrid.add_block(X3d, Y3d, Z3d)

        # Add initial curve as a separate zone
        curve = self.curve.get_points().T
        myGrid.add_block(np.array([[curve[:,0]]]),
                         np.array([[curve[:,1]]]),
                         np.array([[curve[:,2]]]))

        # Export grid
        plot3d_interface.export_plot3d(myGrid, filename)

    #================================================================
    #================================================================
    # OPTIONS METHODS

    def _getDefaultOptions(self):
        """ Define default options and pass back a dict. """

        self.optionsDict = {

            'bc1' : 'splay',
            'bc2' : 'splay',
            'dStart' : 1.e-2,
            'numLayers' : 17,
            'epsE0' : 1.0,
            'theta' : 0.0,
            'alphaP0' : 0.25,
            'numSmoothingPasses' : 0,
            'nuArea' : 0.16,
            'numAreaPasses' : 0,
            'sigmaSplay' : 0.3,
            'cMax' : 3.0,
            'ratioGuess' : 20,
            'plotQuality' : False,
            'guideCurves' : [],
            'remesh' : False,
            }

    def _applyUserOptions(self, options):
        # Override default options with user options
        for userKey in options.keys():

            # Treat special options first
            if userKey.lower() == 'extension':

                # Store extension as the relevant marching parameter
                self.marchParameter = options['extension']

                # Check if we did not flag growth ratio as a parameter previously
                if self.extension_given is None:
                    self.extension_given = True
                else:
                    error('Cannot define extension AND growthRatio parameters. Please select only one to include in the options dictionary.')

            elif userKey.lower() == 'growthratio':

                # Store extension as the relevant marching parameter
                self.marchParameter = options['growthratio']

                # Check if we did not flag growth ratio as a parameter previously
                if self.extension_given is None:
                    self.extension_given = False
                else:
                    error('Cannot define extension AND growthRatio parameters. Please select only one to include in the options dictionary.')

            # All other options are here
            else:

                unusedOption = True
                for defaultKey in self.optionsDict.keys():
                    if userKey.lower() == defaultKey.lower():
                        unusedOption = False
                        self.optionsDict[defaultKey] = options[userKey]
                        break
                if unusedOption:
                    message = "{} key not in default options dictionary.".format(userKey)
                    warn(message)


    #================================================================
    #================================================================
    # TESTING METHODS

    def test_all(self):

        '''
        This method executes all possible tests
        '''

        print ''
        print '#===================================================#'
        print 'Initiating hypSurf tests'
        print ''

        self.test_internal_fortran_subroutines()
        self.test_python_fortran_consistency(plotError=True)
        self.test_forwardAD_FD(plotError=True)
        self.test_forwardAD_reverseAD()

        print ''
        print 'Finalized hypSurf tests'
        print '#===================================================#'
        print ''

    def test_internal_fortran_subroutines(self, fixedSeed=True):

        '''
        This method test the underlying Fortran subroutines used by the
        AD code
        '''

        print ''
        print '#===================================================#'
        print 'Checking internal Fortran subroutines'
        print ''

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # Make sure variables have the correct type
        dStart = float(self.optionsDict['dStart'])
        cMax = float(self.optionsDict['cMax'])
        ratioGuess = float(self.optionsDict['ratioGuess'])
        theta = float(self.optionsDict['theta'])
        bc1 = self.optionsDict['bc1']
        bc2 = self.optionsDict['bc2']
        sigmaSplay = float(self.optionsDict['sigmaSplay'])
        numLayers = int(self.optionsDict['numLayers'])
        epsE0 = float(self.optionsDict['epsE0'])
        theta = float(self.optionsDict['theta'])
        alphaP0 = float(self.optionsDict['alphaP0'])
        numSmoothingPasses = int(self.optionsDict['numSmoothingPasses'])
        numAreaPasses = int(self.optionsDict['numAreaPasses'])
        nuArea = float(self.optionsDict['nuArea'])

        # First we need to solve the current case to get a baseline mesh
        self.createMesh(fortran_flag=True)
        mesh = self.mesh.copy()

        ### DOT PRODUCT TEST FOR SMOOTHING
        coor = mesh[:,:,3].T.flatten()#geom.curves[curve].coor.reshape(-1, order='F')
        rOut0 = hypsurfAPI.hypsurfapi.smoothing(coor, 5., alphaP0, 1, numLayers)

        # FORWARD MODE
        coord = np.random.random_sample(coor.shape)
        coord = coord/np.sqrt(np.sum(coord**2))
        coord_copy = coord.copy()
        rOut, rOutd = hypsurfAPI.hypsurfapi.smoothing_d(coor, coord, 5., alphaP0, 1, numLayers)

        # REVERSE MODE
        rOutb = np.random.random_sample(coor.shape)
        rOutb_copy = rOutb.copy()
        coorb, rOut = hypsurfAPI.hypsurfapi.smoothing_b(coor, 5., alphaP0, 1, numLayers, rOutb)

        print
        print 'Dot product test for Smoothing (this should be around 1e-14):'
        dotprod = np.sum(coord_copy*coorb) - np.sum(rOutd*rOutb_copy)
        print dotprod
        print

        # FINITE DIFFERENCE
        # Choose step size
        stepSize = 1e-7

        # Perturb input variables
        coor = coor + stepSize*coord

        # Run code
        rOut = hypsurfAPI.hypsurfapi.smoothing(coor, 5., alphaP0, 1, numLayers)

        # Compute derivatives
        rOutd_FD = (rOut - rOut0)/stepSize

        # Compare derivatives
        FD_check = np.max(np.abs(rOutd - rOutd_FD))

        # Print results
        print ''
        print 'Finite difference test for Smoothing (this should be around 1e-8):'
        print FD_check
        print ''

        ### DOT PRODUCT TEST FOR PROJECTION

        r0 = mesh[:,:,1].T.flatten()

        # Clean all projection dictionaries
        self.projDict = []
        self.curveProjDict1 = []
        self.curveProjDict2 = []
        for index in range(len(self.curveProjDictGuide)):
            self.curveProjDictGuide[index] = []

        # FORWARD PASS
        rNext, NNext = self.projection(r0)

        # Store projection dictionaries
        projDict = self.projDict[:]
        curveProjDict1 = self.curveProjDict1[:]
        curveProjDict2 = self.curveProjDict2[:]
        curveProjDictGuide = []
        for index in range(len(self.curveProjDictGuide)):
            curveProjDictGuide.append(self.curveProjDictGuide[index][:])

        # FORWARD DERIVATIVES
        r0d = np.array(np.random.rand(r0.shape[0]),order='F')
        r0d = r0d/np.sqrt(np.sum(r0d**2))

        # Set geomtry object seeds
        coord, curveCoord = self.ref_geom.set_randomADSeeds(mode='forward')

        # Propagate derivatives
        rNextd, NNextd = self.projection_d(r0, r0d, rNext, NNext, 0)

        # REVERSE DERIVATIVES
        rNextb = np.random.rand(rNext.shape[0])
        NNextb = np.random.rand(NNext.shape[0],NNext.shape[1])
        rNextb_copy = rNextb.copy()
        NNextb_copy = NNextb.copy()

        self.ref_geom.clean_reverseADSeeds()

        # Restore projection dictionaries
        self.projDict = projDict
        self.curveProjDict1 = curveProjDict1
        self.curveProjDict2 = curveProjDict2

        r0b = self.projection_b(r0, rNext, rNextb, NNext, NNextb, 0)
        rNextb = rNextb_copy
        NNextb = NNextb_copy

        coorb, curveCoorb = self.ref_geom.get_reverseADSeeds(clean=True)

        dotprod = 0.0
        dotprod = dotprod + np.sum(r0d*r0b)
        dotprod = dotprod + np.sum(coord*coorb)
        for curveName in curveCoord:
            dotprod = dotprod + np.sum(curveCoord[curveName]*curveCoorb[curveName])
        dotprod = dotprod - np.sum(rNextd*rNextb)
        dotprod = dotprod - np.sum(NNextd*NNextb)

        print ''
        print 'Dot product test for Projection (this should be around 1e-14):'
        print dotprod
        print ''

        # FORWARD DERIVATIVES WITH FINITE DIFFERENCING

        stepSize = 1e-7

        r0_FD = r0 + r0d*stepSize
        self.ref_geom.update(self.ref_geom.coor + np.array(coord*stepSize,order='F'))
        for curveName in self.ref_geom.curves:
            self.ref_geom.curves[curveName].set_points(self.ref_geom.curves[curveName].get_points() + curveCoord[curveName]*stepSize)

        rNext_FD, NNext_FD = self.projection(r0_FD)

        self.ref_geom.update(self.ref_geom.coor - np.array(coord*stepSize,order='F'))
        for curveName in self.ref_geom.curves:
            self.ref_geom.curves[curveName].set_points(self.ref_geom.curves[curveName].get_points() - curveCoord[curveName]*stepSize)

        rNextd_FD = (rNext_FD-rNext)/stepSize
        NNextd_FD = (NNext_FD-NNext)/stepSize

        FD_check = np.max([np.max(rNextd - rNextd_FD), np.max(NNextd - NNextd_FD)])

        print ''
        print 'Finite difference test for Projection (this should be around 1e-8):'
        print FD_check
        print ''

        ### DOT PRODUCT TEST FOR AREAFACTOR

        r0 = rNext
        d = self.optionsDict['dStart']
        layerindex = 1

        # FORWARD MODE
        r0_d = np.random.random_sample(r0.shape)
        r0_d = r0_d/np.sqrt(np.sum(r0_d**2))
        d_d = 1.0

        r0_d_copy = r0_d.copy()
        d_d_copy = d_d

        s0,s_d,maxstretch =  hypsurfAPI.hypsurfapi.areafactor_test_d(r0, r0_d, d, d_d, nuArea, numAreaPasses, bc1, bc2, self.guideIndices)

        # REVERSE MODE
        s_b = np.random.random_sample(s0.shape)

        s_b_copy = s_b.copy()

        r0_b,d_b = hypsurfAPI.hypsurfapi.areafactor_test_b(r0, d, nuArea, numAreaPasses, bc1, bc2, self.guideIndices, s0, s_b, maxstretch)


        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(r0_d_copy*r0_b)
        dotProd = dotProd + d_d_copy*d_b
        dotProd = dotProd - np.sum(s_d*s_b_copy)

        print ''
        print 'Dot product test for AreaFactor (this should be around 1e-14):'
        print dotProd
        print ''

        # FINITE DIFFERENCE
        # Choose step size
        stepSize = 1e-7

        # Perturb input variables
        r0 = r0 + stepSize*r0_d
        d = d + stepSize*d_d

        # Run code
        s,_,maxstretch =  hypsurfAPI.hypsurfapi.areafactor_test_d(r0, r0_d, d, d_d, nuArea, numAreaPasses, bc1, bc2, self.guideIndices)

        # Compute derivatives
        s_d_FD = (s - s0)/stepSize

        # Compare derivatives
        FD_check = np.max(np.abs(s_d - s_d_FD))

        # Print results
        print ''
        print 'Finite difference test for AreaFactor (this should be around 1e-8):'
        print FD_check
        print ''

        ### DOT PRODUCT TEST FOR DISSIPATIONCOEFFICIENTS

        r0 = rNext
        N0 = NNext
        layerIndex = 1
        iIndex= 3

        r0_xi = 0.5*(r0[3*(iIndex+1):3*(iIndex+1)+3] - r0[3*(iIndex-1):3*(iIndex-1)+3])

        import hypsurf_python as hp
        angle = hp.giveAngle(r0[3*(iIndex-1):3*(iIndex-1)+3],
                             r0[3*(iIndex+0):3*(iIndex+0)+3],
                             r0[3*(iIndex+1):3*(iIndex+1)+3],
                             N0[:,iIndex])

        B0 = np.zeros((3,3))
        B0[0,:] = r0_xi
        B0[1,:] = np.cross(N0[:,iIndex], r0_xi)
        B0[2,:] = N0[:,iIndex]
        Binv = np.linalg.inv(B0)
        r0_eta = Binv[:,2]*s0[iIndex]
        d_sensor = 1.0

        # FORWARD MODE
        r0_xid = np.array(np.random.random_sample(r0_xi.shape),order='F')
        r0_xid = r0_xid/np.sqrt(np.sum(r0_xid**2))
        r0_etad = np.array(np.random.random_sample(r0_eta.shape),order='F')
        r0_etad = r0_etad/np.sqrt(np.sum(r0_etad**2))
        d_sensord = 1.0
        angled = 1.0

        epsEstart, epsEd, epsIstart, epsId =  hypsurfAPI.hypsurfapi.dissipationcoefficients_d_test(layerIndex,
                                                                                                   r0_xi, r0_xid,
                                                                                                   r0_eta, r0_etad,
                                                                                                   d_sensor, d_sensord,
                                                                                                   angle, angled,
                                                                                                   numLayers, epsE0)

        # REVERSE MODE
        epsEb = 0.5
        epsIb = 0.3

        epsEb_copy = epsEb
        epsIb_copy = epsIb

        r0_xib, r0_etab, d_sensorb, angleb = hypsurfAPI.hypsurfapi.dissipationcoefficients_b_test(layerIndex,
                                                                                                  r0_xi, r0_eta,
                                                                                                  d_sensor, angle,
                                                                                                  numLayers,
                                                                                                  epsE0,
                                                                                                  epsEstart, epsEb,
                                                                                                  epsIstart, epsIb)


        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(r0_xid*r0_xib)
        dotProd = dotProd + np.sum(r0_etad*r0_etab)
        dotProd = dotProd + d_sensord*d_sensorb
        dotProd = dotProd + angled*angleb
        dotProd = dotProd - epsEd*epsEb_copy
        dotProd = dotProd - epsId*epsIb_copy

        print ''
        print 'Dot product test for DissipationCoefficients (this should be around 1e-14):'
        print dotProd
        print ''

        # FINITE DIFFERENCE
        # Choose step size
        stepSize = 1e-7

        # Perturb input variables
        r0_xi = r0_xi + stepSize*r0_xid
        r0_eta = r0_eta + stepSize*r0_etad
        d_sensor = d_sensor + stepSize*d_sensord
        angle = angle + stepSize*angled

        # Run code
        epsE, epsI =  hypsurfAPI.hypsurfapi.dissipationcoefficients_test(layerIndex,
                                                                         r0_xi,
                                                                         r0_eta,
                                                                         d_sensor,
                                                                         angle,
                                                                         numLayers, epsE0)

        # Compute derivatives
        epsEd_FD = (epsE - epsEstart)/stepSize
        epsId_FD = (epsI - epsIstart)/stepSize

        # Compare derivatives
        FD_check = np.max(np.abs([epsEd - epsEd_FD,
                                  epsId - epsId_FD]))

        # Print results
        print ''
        print 'Finite difference test for DissipationCoefficients (this should be around 1e-8):'
        print FD_check
        print ''

        ### DOT PRODUCT TEST FOR MATRIXBUILDER

        r0 = rNext[:]
        rm1 = rNext[:]
        n0 = NNext[:,:]
        s0 = s
        sm1 = s
        layerindex = 4
        iindex = 3

        numLayers = 30
        epsE0 = 5.5
        layerIndex = 2
        theta = 0.0
        numNodes = 5

        K = np.zeros((len(r0),len(r0)),order='F')
        f = np.zeros(len(r0),order='F')

        # Generate normalized random seeds
        r0_d = np.random.random_sample(r0.shape)
        r0_d = r0_d/np.sqrt(np.sum(r0_d**2))

        n0_d = np.random.random_sample(n0.shape)
        n0_d = n0_d/np.sqrt(np.sum(n0_d**2))

        s0_d = np.random.random_sample(s0.shape)
        s0_d = s0_d/np.sqrt(np.sum(s0_d**2))

        rm1_d = np.random.random_sample(rm1.shape)
        rm1_d = rm1_d/np.sqrt(np.sum(rm1_d**2))

        sm1_d = np.random.random_sample(sm1.shape)
        sm1_d = sm1_d/np.sqrt(np.sum(sm1_d**2))

        K_d = np.zeros((len(r0),len(r0)),order='F')
        f_d = np.zeros(len(r0),order='F')

        hypsurfAPI.hypsurfapi.matrixbuilder_test(iIndex, bc1, bc2, r0, rm1, n0, s0, sm1,
                                                 numLayers, epsE0, layerIndex, theta, K, f)

        K0 = K[:,:]
        f0 = f[:]

        hypsurfAPI.hypsurfapi.matrixbuilder_d_test(iIndex, bc1, bc2, r0, r0_d, rm1, rm1_d, n0, n0_d, s0, s0_d, sm1, sm1_d,
                                                   numLayers, epsE0, layerIndex, theta, K, K_d, f, f_d)

        print ''
        print 'matrixbuilder test'
        print np.max(np.abs(K-K0))
        print np.max(np.abs(f-f0))
        print ''

        ### DOT PRODUCT TEST FOR COMPUTEMATRICES
        r0 = rNext
        rm1 = rNext
        n0 = NNext
        s0 = s
        sm1 = s
        layerindex = 1

        # FORWARD MODE

        # Generate normalized random seeds
        r0_d = np.random.random_sample(r0.shape)
        r0_d = r0_d/np.sqrt(np.sum(r0_d**2))

        n0_d = np.random.random_sample(n0.shape)
        n0_d = n0_d/np.sqrt(np.sum(n0_d**2))

        s0_d = np.random.random_sample(s0.shape)
        s0_d = s0_d/np.sqrt(np.sum(s0_d**2))

        rm1_d = np.random.random_sample(rm1.shape)
        rm1_d = rm1_d/np.sqrt(np.sum(rm1_d**2))

        sm1_d = np.random.random_sample(sm1.shape)
        sm1_d = sm1_d/np.sqrt(np.sum(sm1_d**2))

        # Copy values to avoid any modification done in Fortran
        r0_d_copy = r0_d.copy()
        n0_d_copy = n0_d.copy()
        s0_d_copy = s0_d.copy()
        rm1_d_copy = rm1_d.copy()
        sm1_d_copy = sm1_d.copy()

        retainSpacing = self.optionsDict['remesh']

        # Run the forward AD
        rnext0, rnext_d = hypsurfAPI.hypsurfapi.computematrices_d(r0, r0_d, n0, n0_d, s0, s0_d, rm1, rm1_d, sm1, sm1_d, layerindex, theta, sigmaSplay, bc1, bc2, self.guideIndices, retainSpacing, numLayers, epsE0)

        # REVERSE MODE
        rnext_b = np.random.random_sample(r0.shape)

        rnext_b_copy = rnext_b.copy()

        r0_b, n0_b, s0_b, rm1_b, sm1_b, rnext = hypsurfAPI.hypsurfapi.computematrices_b(r0, n0, s0, rm1, sm1, layerindex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, self.guideIndices, retainSpacing, rnext_b)

        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(rnext_d*rnext_b_copy)
        dotProd = dotProd - np.sum(r0_b*r0_d_copy)
        dotProd = dotProd - np.sum(n0_b*n0_d_copy)
        dotProd = dotProd - np.sum(s0_b*s0_d_copy)
        dotProd = dotProd - np.sum(rm1_b*rm1_d_copy)
        dotProd = dotProd - np.sum(sm1_b*sm1_d_copy)

        print ''
        print 'Dot product test for ComputeMatrices (this should be around 1e-14):'
        print dotProd
        print ''

        # FINITE DIFFERENCE
        # Define step size
        stepSize = 1e-7

        # Perturb input variables
        r0 = r0 + r0_d*stepSize
        n0 = n0 + n0_d*stepSize
        s0 = s0 + s0_d*stepSize
        rm1 = rm1 + rm1_d*stepSize
        sm1 = sm1 + sm1_d*stepSize

        # Run code
        rnext = hypsurfAPI.hypsurfapi.computematrices(r0, n0, s0, rm1, sm1, layerindex, theta, sigmaSplay, bc1, bc2, self.guideIndices, retainSpacing, numLayers, epsE0)

        # Compute derivatives
        rnext_d_FD = (rnext - rnext0)/stepSize

        # Compute maximum error in the FD code
        FD_check = np.max(np.abs(rnext_d - rnext_d_FD))

        print ''
        print 'Finite difference test for ComputeMatrices (this should be around 1e-8):'
        print FD_check
        print ''

        print '#===================================================#'
        print ''

        #quit()

        # Final cleanup
        self.clean()

    def test_python_fortran_consistency(self, plotError=False):

        from time import time

        print ''
        print '#===================================================#'
        print 'Checking consistency between the Python and Fortran versions of hypsurf'
        print ''

        # Call the Fortran version
        st = time()
        self.createMesh(fortran_flag=True)
        self.exportPlot3d('output_fortran.xyz')
        mesh_fortran = self.mesh.copy()
        print 'Fortran version took: ',time() - st, 'secs'

        # Call the Python version
        st = time()
        self.createMesh(fortran_flag=False)
        self.exportPlot3d('output_python.xyz')
        mesh_python = self.mesh.copy()

        print 'Python version took: ',time() - st, 'secs'

        # Compare differences
        differences = np.abs(mesh_fortran-mesh_python)
        error = np.max(differences)

        print
        print 'Maximum difference between Fortran and Python (should be around 1e-14):'
        print error
        print

        # Show matrix of errors if max error is too large
        if plotError or error > 1e-9:

            # Reshape 3D array to 2D so we can plot it
            differences = differences.T.reshape((self.optionsDict['numLayers'],-1))

            # Plot error matrix
            view_mat(differences)

        print '#===================================================#'
        print ''

    def test_forwardAD_FD(self, stepSize=1e-7, fixedSeed=True, plotError=False):

        '''
        This method compares derivatives computed with automatic differentiation (AD)
        with derivatives computed by finite differencing (FD).

        fixedSeed: boolean -> Determines if we set a fixed seed to the random number
        generator so that we always get the same derivatives. This helps debugging as
        values will not change every time you run the same code.

        Ney Secco 2017-01
        '''

        # INITIALIZATION

        print ''
        print '#===================================================#'
        print 'Checking forward AD derivatives with Finite Differences'
        print ''

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # REFERENCE POINT

        # Run the initial mesh in Fortran
        self.createMesh(fortran_flag=True)

        # Store the initial mesh
        mesh0 = self.mesh[:,:,:]

        # DERIVATIVE SEEDS

        # Triangulated surface nodes seeds
        coord, curveCoord = self.ref_geom.set_randomADSeeds(mode='forward')

        # Initial curve seeds
        initCurveCoord = self.curve.set_randomADSeeds(mode='forward')
        initCurveCoor0 = self.curve.get_points()

        # AD VERSION

        # Propagate derivatives using AD
        self.compute_forwardAD()

        # Get output derivatives
        meshd_AD = self.get_forwardADSeeds()

        # FD VERSION

        # Perturb nodes based on derivative seeds
        self.ref_geom.update(self.ref_geom.coor + coord*stepSize)
        for curveName in self.ref_geom.curves:
            curve = self.ref_geom.curves[curveName]
            currPts = curve.get_points()
            newPts = currPts + curveCoord[curveName]*stepSize
            curve.set_points(newPts)

        # Perturb initial curve
        self.curve.set_points(initCurveCoor0 + stepSize*initCurveCoord)

        # Run the perturbed mesh in Fortran
        self.createMesh(fortran_flag=True)

        # Store the perturbed mesh
        mesh = self.mesh[:,:,:]

        # Restore initial nodes
        self.ref_geom.update(self.ref_geom.coor - coord*stepSize)
        for curveName in self.ref_geom.curves:
            curve = self.ref_geom.curves[curveName]
            currPts = curve.get_points()
            newPts = currPts - curveCoord[curveName]*stepSize
            curve.set_points(newPts)

        self.curve.set_points(initCurveCoor0)

        # Compute derivatives with Finite Differences
        meshd_FD = (mesh-mesh0)/stepSize

        # Get differences between both versions
        differences = np.zeros((self.optionsDict['numLayers'], self.numNodes*3))
        differences[:,::3] = (meshd_AD[0,:,:] - meshd_FD[0,:,:]).T
        differences[:,1::3] = (meshd_AD[1,:,:] - meshd_FD[1,:,:]).T
        differences[:,2::3] = (meshd_AD[2,:,:] - meshd_FD[2,:,:]).T
        error = np.max(np.abs(differences))

        # Print result
        print ''
        print 'Maximum discrepancy between AD and FD derivatives (should be around 1e-7):'
        print error
        print ''

        # Plot errors if requested by the user or if error is too large
        if plotError or error > 1e-4:
            view_mat(np.abs(differences))

        print '#===================================================#'
        print ''

    def test_forwardAD_reverseAD(self, fixedSeed=True):

        '''
        This method uses the dot product test to verify if the forward and
        reverse AD versions are consistent.

        Ney Secco 2017-01
        '''

        # INITIALIZATION

        print ''
        print '#===================================================#'
        print 'Checking forward and reverse AD with the dot product test'
        print ''

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # REFERENCE POINT

        # Run the initial mesh in Fortran
        self.createMesh(fortran_flag=True)

        # Store the initial mesh
        mesh0 = self.mesh[:,:,:]

        # DERIVATIVE SEEDS

        # Generate a set of random derivative seeds

        # Triangulated surface nodes seeds
        coord, curveCoord = self.ref_geom.set_randomADSeeds(mode='forward')

        # Initial curve seeds
        initCurveCoord = self.curve.set_randomADSeeds(mode='forward')

        # Surface mesh seeds
        meshb = self.set_randomADSeeds(mode='reverse')

        # FORWARD AD VERSION

        # Propagate derivatives using forward AD
        self.compute_forwardAD()

        # Get output derivatives
        meshd = self.get_forwardADSeeds()

        # REVERSE AD VERSION

        # Back-propagate derivatives using reverse AD
        self.compute_reverseAD()

        # Get input derivatives
        coorb, curveCoorb = self.ref_geom.get_reverseADSeeds(clean=False)

        initCurveCoorb = self.curve.get_reverseADSeeds(clean=False)

        # Only store the initial curve seed if it does not belong to the geometry object.
        # Otherwise, the initial curve seeds are already stored inside the geometry object
        for geomCurve in self.ref_geom.curves.itervalues():
            if self.curve is geomCurve:
                initCurveCoorb = 0.0

        # DOT PRODUCT TEST

        dotProduct = 0.0
        dotProduct = dotProduct + np.sum(coord*coorb)
        for curveName in self.ref_geom.curves:
            dotProduct = dotProduct + np.sum(curveCoord[curveName]*curveCoorb[curveName])
        dotProduct = dotProduct + np.sum(initCurveCoord*initCurveCoorb)
        dotProduct = dotProduct - np.sum(meshd*meshb)

        print ''
        print ' Marching dot product test (this should be around 1e-14):'
        print dotProduct
        print ''
        print '#===================================================#'
        print ''


'''
==============================================
MORE AUXILIARY FUNCTIONS
==============================================
'''

#=============================================
#=============================================

def plotGrid(X,Y,Z,show=False):

    # This function plots the grid

    # IMPORTS
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Get number of nodes and number of layers
    num_xi, num_eta = X.shape

    # Initialize figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot eta lines
    for i in range(num_eta):
        ax.plot(X[:,i], Y[:,i], Z[:,i],'k')

    # Plot xi lines
    for i in range(num_xi):
        ax.plot(X[i,:], Y[i,:], Z[i,:],'k')

    # Plot the initial curve
    ax.plot(X[0,:],Y[0,:],Z[0,:],'r',linewidth=2)

    # Add labels
    plt.xlabel('x')
    plt.ylabel('y')

    # Adjust axis
    ax.axis('equal')
    ax.set_axis_off()

    if show:
        # Show
        plt.show()

    # Return the figure handle
    return fig

def warn(message):
    print ''
    print 'WARNING!:'
    print message
    print ''

def error(message):
    print ''
    print 'ERROR!:'
    print message
    print ''

#=============================================
#=============================================

def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

def closest_node(guideCurve, curve):
    """ Find closest node from a list of node coordinates. """
    curve = np.asarray(curve.get_points()).T

    # Get number of points in the seed curve
    nPoints = curve.shape[0]

    # Initialize arrays to call projection function
    dist2 = np.ones(nPoints)*1e10
    xyzProj = np.zeros((nPoints,3))
    tangents = np.zeros((nPoints,3))
    elemIDs = np.zeros((nPoints),dtype='int32')

    # Call projection function to find the distance between every node of the
    # seed curve to the guide curve
    guideCurve.project(curve, dist2, xyzProj, tangents, elemIDs)

    # Find closest point
    ind = np.argmin(dist2)

    return ind
