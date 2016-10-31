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

fortran_flag = True
deriv_check = True

'''
TO DO

- add subiteration loop
- check smoothing
- set up scaling factor based on distance
- blend the angle-based dissipation coefficient
'''


class HypSurfMesh(object):
    """
    Hyperbolic surface mesh created by projecting points extruded from a
    curve onto a surface.

    """

    def __init__(self, curve, ref_geom, options={}):

        # Check to see if the marching extension (distance) is given
        try:
            _ = options['extension']
            self.extension_given = True
        except:
            self.extension_given = False
            pass

        # Check to see if the growth ratio is given
        try:
            _ = options['growthRatio']
            self.growthRatio_given = True
        except:
            self.growthRatio_given = False
            pass

        self._getDefaultOptions()
        self._applyUserOptions(options)

        if isinstance(curve, np.ndarray):
            self.curve = curve
        else:
            # We assume that curve is a string that defines a curve name
            self.curve = ref_geom.curves[curve].extract_points()

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

        # Get the number of nodes
        self.numNodes = self.curve.shape[0]
        self.mesh = np.zeros((3, self.numNodes, self.optionsDict['numLayers']))

    def createMesh(self):
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


        dStart = self.optionsDict['dStart']
        cMax = self.optionsDict['cMax']
        ratioGuess = self.optionsDict['ratioGuess']
        theta = self.optionsDict['theta']
        bc1 = self.optionsDict['bc1']
        bc2 = self.optionsDict['bc2']
        sigmaSplay = self.optionsDict['sigmaSplay']
        numLayers = self.optionsDict['numLayers']
        epsE0 = self.optionsDict['epsE0']
        theta = self.optionsDict['theta']
        alphaP0 = self.optionsDict['alphaP0']
        numSmoothingPasses = self.optionsDict['numSmoothingPasses']
        numAreaPasses = self.optionsDict['numAreaPasses']
        nuArea = self.optionsDict['nuArea']

        if self.growthRatio_given and self.extension_given:
            error('Cannot define extension and growthRatio parameters. Please select only one to include in the options dictionary.')

        if self.extension_given:
            extension = self.optionsDict['extension']
            marchParameter = extension
        else:
            growthRatio = self.optionsDict['growthRatio']
            marchParameter = growthRatio

        # Flatten the coordinates vector
        # We do this because the linear system is assembled assuming that
        # the coordinates vector is flattened
        rStart = self.curve.flatten().astype(float)

        if fortran_flag:

            # Perform the marching algorithm and output the results into R,
            # which contains the mesh.
            # fail is a flag set to true if the marching algo failed
            # ratios is the ratios of quality for the mesh
            R, fail, ratios, majorIndices = hypsurfAPI.hypsurfapi.march(self.projection, rStart, dStart, theta, sigmaSplay, bc1.lower(), bc2.lower(), epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, self.extension_given, numSmoothingPasses, numAreaPasses, numLayers)

            # Obtain the pseudomesh, or subiterations mesh from the three stages of marching.
            # These are used in the adjoint formulation.
            R_initial_march = np.array(hypsurfAPI.hypsurfapi.r_initial_march)
            R_smoothed = np.array(hypsurfAPI.hypsurfapi.r_smoothed)
            R_final = np.array(hypsurfAPI.hypsurfapi.r_final)
            S = np.array(hypsurfAPI.hypsurfapi.s)
            N = np.array(hypsurfAPI.hypsurfapi.n)

            # Derivative check
            if deriv_check:
                rStartd = np.random.random_sample(rStart.shape)

                rStartd_copy = rStartd.copy()
                R_, Rd, fail, ratios, _ = hypsurfAPI.hypsurfapi.march_d(self.projection, rStart, rStartd, dStart, theta, sigmaSplay, bc1.lower(), bc2.lower(), epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, self.extension_given, numSmoothingPasses, numAreaPasses, numLayers)

                Rb = np.random.random_sample(R.shape)
                Rb_copy = Rb.copy()

                rStartb, fail = hypsurfAPI.hypsurfapi.march_b(self.projection, rStart, R_initial_march, R_smoothed, R_final, N, majorIndices, dStart, theta, sigmaSplay, bc1.lower(), bc2.lower(), epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, self.extension_given, numSmoothingPasses, numAreaPasses, R, Rb, ratios, numLayers)

                print ' Marching dot product test, this should be zero:', np.sum(Rd*Rb_copy) - np.sum(rStartd_copy*rStartb), '(unless the surface is curved; don\'t have projections working)'
                print

            # Release the pseudomesh information from the hypsurfAPI instance
            hypsurfAPI.hypsurfapi.releasememory()

        else:

            # Initialize 2D array that will contains all the surface grid points in the end
            R = np.zeros((numLayers,len(rStart)))

            # Project onto the surface or curve (if applicable)
            rNext, NNext = self.projection(rStart)

            # Initialize step size and total marched distance
            d = dStart
            dTot = 0

            if self.extension_given:
                # Find the characteristic radius of the mesh
                radius = findRadius(rNext)

                # Find the desired marching distance
                dMax = radius*(extension-1)

                # Compute the growth ratio necessary to match this distance
                dGrowth = findRatio(dMax, dStart, numLayers, ratioGuess)
            else:
                dGrowth = growthRatio

            # Print growth ratio
            print 'Growth ratio: ',dGrowth

            # Store the initial curve
            R[0,:] = rNext

            # Issue a warning message if the projected points are far from the given points
            if max(abs(rNext - rStart)) > 1.0e-5:
                warn('The given points (rStart) might not belong to the given surface (surf)\nThese points were projected for the surface mesh generation')


            # We need a guess for the first-before-last curve in order to compute the grid distribution sensor
            # As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
            rm1 = rNext[:]

            #===========================================================

            '''
            The marching function actually begins here
            '''
            # Some functions require the area factors of the first-before-last curve
            # We will repeat the first curve areas for simplicity.
            # rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
            # But we still need to find Sm1
            Sm1, maxStretch = self.areaFactor(rNext, d)

            fail = False

            # MARCH!!!
            for layerIndex in range(numLayers-1):

                # Get the coordinates computed by the previous iteration
                r0 = rNext[:]

                # Compute the new area factor for the desired marching distance
                S0, maxStretch = self.areaFactor(r0, d)

                # The subiterations will use pseudo marching steps.
                # If the required marching step is too large, the hyperbolic marching might become
                # unstable. So we subdivide the current marching step into multiple smaller pseudo-steps

                # Compute the factor between the current stretching ratio and the allowed one.
                # If the current stretching ratio is smaller than cMax, the cFactor will be 1.0, and
                # The pseudo-step will be the same as the desired step.
                cFactor = int(np.ceil(maxStretch/cMax))

                # Constrain the marching distance if the stretching ratio is too high
                dPseudo = d/cFactor

                # Subiteration
                # The number of subiterations is the one required to meet the desired marching distance
                for indexSubIter in range(cFactor):

                    # Recompute areas with the pseudo-step
                    S0, maxStretch = self.areaFactor(r0, dPseudo)

                    # Update the Normals with the values computed in the last iteration.
                    # We do this because, the last iteration already projected the new
                    # points to the surface and also computed the normals. So we don't
                    # have to repeat the projection step
                    N0 = NNext[:,:]

                    # March using the pseudo-marching distance
                    eta = layerIndex+2
                    rNext, NNext = self.subIteration(r0, N0, S0, rm1, Sm1, layerIndex)

                    # Update Sm1 (Store the previous area factors)
                    Sm1 = S0[:]

                    # Update rm1
                    rm1 = r0[:]

                    # Update r0
                    r0 = rNext[:]

                # Store grid points
                R[layerIndex+1,:] = rNext

                # Check quality of the mesh
                if layerIndex > 1:
                    fail, ratios = self.qualityCheck(R[layerIndex-2:layerIndex+2, :], layerIndex)


                if fail:
                    # If the mesh is not valid, only save the mesh up until that point.
                    # Otherwise we'd see a bunch of points at 0,0,0.
                    R = R[:layerIndex+2, :]
                    break

                # Compute the total marched distance so far
                dTot = dTot + d

                # Update step size
                d = d*dGrowth

            if self.optionsDict['plotQuality']:
                fail, ratios = self.qualityCheck(R)

        print time() - st, 'secs'

        if self.optionsDict['plotQuality']:
            view_mat(ratios)


        # Convert to X, Y and Z
        X = R[:,::3]
        Y = R[:,1::3]
        Z = R[:,2::3]

        self.mesh = np.zeros((3, self.numNodes, numLayers))
        self.mesh[0, :, :] = X.T
        self.mesh[1, :, :] = Y.T
        self.mesh[2, :, :] = Z.T


    def subIteration(self, r0, N0, S0, rm1, Sm1, layerIndex):

        '''
        r0 -> 1D array (1 x 3*numNodes): flattened vector containing coordinates of the latest curve
              r0 = [x1 y1 z1 x2 y2 z2 x3 y3 z3 ... ]

        N0 -> 2D array (3 x numNodes): Surface normals at every node of the latest curve (r0)
              N0 = [ nx1 nx2 nx3 nx4 ... ]
                   [ ny1 ny2 ny3 ny4 ... ]
                   [ nz1 nz2 nz3 nz4 ... ]

        S0 -> 1D array (1 x numNodes): Area factor (related to marching distance) of every node of
                                       the latest curve (r0). The area is the marching distance times
                                       the distance between the node and its neighbors
              S0 = [ S1 S2 S3 S4 ... ]

        rm1 -> 1D array (1 x 3*numNodes): flattened vector containing coordinates of the
                                          first-before-last curve
               rm1 = [x1m1 y1m1 z1m1 x2m1 y2m1 z2m1 x3m1 y3m1 z3m1 ... ]

        Sm1 -> 1D array (1 x numNodes): Area factor (related to marching distance) of every node of
                                        the first-before-last curve (rm1)
              S0 = [ S1m1 S2m1 S3m1 S4m1 ... ]

        layerIndex -> integer : 0-based index corresponding to the current layer (the starting curve
                                 has layerIndex = 0, the first generated curve has layerIndex = 1,
                                 and so on...

        '''


        theta = self.optionsDict['theta']
        bc1 = self.optionsDict['bc1']
        bc2 = self.optionsDict['bc2']
        sigmaSplay = self.optionsDict['sigmaSplay']
        numLayers = self.optionsDict['numLayers']
        epsE0 = self.optionsDict['epsE0']
        theta = self.optionsDict['theta']
        alphaP0 = self.optionsDict['alphaP0']
        numSmoothingPasses = self.optionsDict['numSmoothingPasses']
        eta = layerIndex+2

        dr = self.computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex)


        # Update r
        rNext = r0 + dr

        # Smooth coordinates
        rNext_ = self.smoothing(rNext,layerIndex+2)

        rNext, NNext = self.projection(rNext_)

        # RETURNS
        return rNext, NNext

    def projection(self, rNext):

        # Initialize normals array
        NNext = np.zeros((3, self.numNodes))

        # Save endpoints
        node1 = rNext[:3]
        node2 = rNext[-3:]

        # Project onto surface and compute surface normals
        rNext, NNext, projDict = self.ref_geom.project_on_surface(rNext.reshape((self.numNodes, 3)))
        rNext = rNext.flatten()
        NNext = NNext.T

        # Replace end points if we use curve BCs
        if self.optionsDict['bc1'].lower().startswith('curve'):
            rNext[:3], NNextAux = self.ref_geom.project_on_curve(node1.reshape((1, 3)), curveCandidates=[self.ref_curve1])
            NNext[:, 0] = NNextAux.T[:, 0]

        if self.optionsDict['bc2'].lower().startswith('curve'):
            rNext[-3:], NNextAux = self.ref_geom.project_on_curve(node2.reshape((1, 3)), curveCandidates=[self.ref_curve2])
            NNext[:, -1] = NNextAux.T[:, 0]

        return rNext, NNext

    def areaFactor(self, r0, d):

        '''
        This function computes the area factor distribution for the current curve.
        The area factor is the marched distance times the distance between the node
        and its neighbors.
        This function also computes the stretching ratio of the desired marching distance
        (ratio between marching distance and the neighbor distance), so that we can avoid
        large steps

        INPUTS:
        r0 -> 1d array (1 x 3*numNodes) : coordinates of all points in a layer:
        r0 = [x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]

        d -> float : desired marching distance
        '''

        nuArea = self.optionsDict['nuArea']
        numAreaPasses = self.optionsDict['numAreaPasses']

        # Extrapolate the end points
        r0minus1 = 2*r0[0:3] - r0[3:6]
        r0plus1 = 2*r0[-3:] - r0[-6:-3]
        r0_extrap = np.hstack([r0minus1, r0, r0plus1])

        # Reshape so we have 3D array
        R0_extrap = np.reshape(r0_extrap,(-1,3)).T

        # Compute the distance of each node to its neighbors
        neighborDist = 0.5*(np.linalg.norm(R0_extrap[:,1:-1] - R0_extrap[:,:-2],axis=0) + np.linalg.norm(R0_extrap[:,2:] - R0_extrap[:,1:-1],axis=0))

        # Multiply distances by the step size to get the areas
        S = d*neighborDist

        # Divide the marching distance and the neighbor distance to get the stretch ratios
        stretchRatio = d/neighborDist

        # Get the maximum stretch ratio
        maxStretch = np.max(stretchRatio)

        # Do the requested number of averagings
        for index in xrange(numAreaPasses):
            # Store previous values
            Splus = S[1]
            Sminus = S[-2]
            # Do the averaging for the central nodes
            S[1:-1] = (1-nuArea)*S[1:-1] + nuArea/2*(S[:-2] + S[2:])
            # Average for the extremum nodes
            S[0] = (1-nuArea)*S[0] + nuArea*Splus
            S[-1] = (1-nuArea)*S[-1] + nuArea*Sminus

        # If we use curve boundary conditions, we need just the marching distance, and not area, for the end nodes
        if self.optionsDict['bc1'].lower().startswith('curve'):
            S[0] = d
        if self.optionsDict['bc2'].lower().startswith('curve'):
            S[-1] = d

        # RETURNS
        return S, maxStretch

    def computeMatrices(self, r0, N0, S0, rm1, Sm1, layerIndex):

        '''
        This function computes the derivatives r_zeta and r_eta for a group
        of coordinates given in r.
        It should be noted that r0 is a 1d array that contains x and y
        coordinates of all points in a layer:
        r0 = [x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]

        N0 is a 2D array with the surface normals at each r0 point
        N0 = [ nx1 nx2 nx3 nx4 ... ]
             [ ny1 ny2 ny3 ny4 ... ]
             [ nz1 nz2 nz3 nz4 ... ]
        '''

        # Get options
        theta = self.optionsDict['theta']
        bc1 = self.optionsDict['bc1']
        bc2 = self.optionsDict['bc2']
        sigmaSplay = self.optionsDict['sigmaSplay']
        numLayers = self.optionsDict['numLayers']
        epsE0 = self.optionsDict['epsE0']

        # Initialize arrays
        K = np.zeros((3*self.numNodes, 3*self.numNodes))
        f = np.zeros(3*self.numNodes)

        def matrixBuilder(curr_index):

            if curr_index == 0:  # forward case

                if bc1 != 'continuous':

                    neighbor1_index = 1
                    neighbor2_index = 2

                    # Using forward differencing for xi = 1
                    r0_xi = 0.5*(-3*r0[3*(curr_index):3*(curr_index)+3] + 4*r0[3*(neighbor1_index):3*(neighbor1_index)+3] - r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                    angle = np.pi

                else:

                    neighbor1_index = self.numNodes - 2
                    neighbor2_index = curr_index + 1

                    # Using central differencing for zeta = 2:numNodes-1
                    r0_xi = 0.5*(r0[3*(neighbor2_index):3*(neighbor2_index)+3] - r0[3*(neighbor1_index):3*(neighbor1_index)+3])

                    # Compute the local grid angle based on the neighbors
                    angle = giveAngle(r0[3*(neighbor1_index):3*(neighbor1_index)+3],
                                      r0[3*(curr_index):3*(curr_index)+3],
                                      r0[3*(neighbor2_index):3*(neighbor2_index)+3],
                                      N0[:,curr_index])


            elif curr_index == self.numNodes - 1:  # backward case
                neighbor1_index = curr_index - 1
                neighbor2_index = curr_index - 2

                # Using backward differencing for xi = numNodes
                r0_xi = 0.5*(3*r0[3*(curr_index):3*(curr_index)+3] - 4*r0[3*(neighbor1_index):3*(neighbor1_index)+3] + r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                angle = np.pi

            else:  # central case
                neighbor1_index = curr_index - 1
                neighbor2_index = curr_index + 1

                # Using central differencing for zeta = 2:numNodes-1
                r0_xi = 0.5*(r0[3*(neighbor2_index):3*(neighbor2_index)+3] - r0[3*(neighbor1_index):3*(neighbor1_index)+3])

                # Compute the local grid angle based on the neighbors
                angle = giveAngle(r0[3*(neighbor1_index):3*(neighbor1_index)+3],
                                  r0[3*(curr_index):3*(curr_index)+3],
                                  r0[3*(neighbor2_index):3*(neighbor2_index)+3],
                                  N0[:,curr_index])
                if np.isnan(angle):
                    pdb.set_trace()

            x0_xi = r0_xi[0]
            y0_xi = r0_xi[1]
            z0_xi = r0_xi[2]

            # Get current normal
            nx = N0[0,curr_index]
            ny = N0[1,curr_index]
            nz = N0[2,curr_index]

            # Assemble B0 matrix
            B0 = np.array([[x0_xi, y0_xi, z0_xi],
                        [ny*z0_xi-nz*y0_xi, nz*x0_xi-nx*z0_xi, nx*y0_xi-ny*x0_xi],
                        [nx, ny, nz]])

            # Invert B0
            B0inv = np.linalg.inv(B0)

            # Compute eta derivatives
            r0_eta = B0inv.dot(np.array([0, Sm1[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = np.array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = np.linalg.norm(rm1[3*(neighbor2_index):3*(neighbor2_index)+3]-rm1[3*(curr_index):3*(curr_index)+3]) + np.linalg.norm(rm1[3*(neighbor1_index):3*(neighbor1_index)+3]-rm1[3*(curr_index):3*(curr_index)+3])
            dden = np.linalg.norm(r0[3*(neighbor2_index):3*(neighbor2_index)+3]-r0[3*(curr_index):3*(curr_index)+3]) + np.linalg.norm(r0[3*(neighbor1_index):3*(neighbor1_index)+3]-r0[3*(curr_index):3*(curr_index)+3])
            dSensor = dnum/dden

            # Sharp convex corner detection
            if angle < 70*np.pi/180: # Corner detected

                # Populate matrix with Eq 8.3
                K[3*curr_index:3*curr_index+3,3*neighbor2_index:3*neighbor2_index+3] = -np.eye(3)
                K[3*curr_index:3*curr_index+3,3*curr_index:3*curr_index+3] = 2*np.eye(3)
                K[3*curr_index:3*curr_index+3,3*neighbor1_index:3*neighbor1_index+3] = -np.eye(3)
                f[3*curr_index:3*curr_index+3] = np.array([0,0,0])

            else:

                # Compute C0 = B0inv*A0
                C0 = B0inv.dot(A0)

                # Compute smoothing coefficients
                epsE, epsI = self.dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle)

                # Compute RHS components
                B0invg = B0inv.dot(np.array([0, S0[curr_index], 0]))

                if curr_index == 0:

                    if self.optionsDict['bc1'] != 'continuous':# forwards
                        De = epsE*(r0[3*(curr_index):3*(curr_index)+3] - 2*r0[3*(neighbor1_index):3*(neighbor1_index)+3] + r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                        # Compute block matrices
                        L_block = -0.5*(1+theta)*C0 - epsI*np.eye(3)
                        M_block = 2*(1+theta)*C0 + 2*epsI*np.eye(3)
                        N_block = -1.5*(1+theta)*C0 + (1-epsI)*np.eye(3)
                        f_block = B0invg + De

                        # Populate matrix
                        K[3*(curr_index):3*(curr_index)+3,3*(neighbor2_index):3*(neighbor2_index)+3] = L_block
                        K[3*(curr_index):3*(curr_index)+3,3*(neighbor1_index):3*(neighbor1_index)+3] = M_block
                        K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = N_block
                        f[3*(curr_index):3*(curr_index)+3] = f_block

                    else:

                        De = epsE*(r0[3*(neighbor1_index):3*(neighbor1_index)+3] - 2*r0[3*(curr_index):3*(curr_index)+3] + r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                        # Compute block matrices
                        L_block = -0.5*(1+theta)*C0 - epsI*np.eye(3)
                        M_block = (1 + 2*epsI)*np.eye(3)
                        N_block = 0.5*(1+theta)*C0 - epsI*np.eye(3)
                        f_block = B0invg + De

                        # Populate matrix
                        K[3*(curr_index):3*(curr_index)+3,3*(neighbor1_index):3*(neighbor1_index)+3] = L_block
                        K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = M_block
                        K[3*(curr_index):3*(curr_index)+3,3*(neighbor2_index):3*(neighbor2_index)+3] = N_block
                        f[3*(curr_index):3*(curr_index)+3] = f_block

                elif curr_index == self.numNodes - 1:  # backwards
                    De = epsE*(r0[3*(curr_index):3*(curr_index)+3] - 2*r0[3*(neighbor1_index):3*(neighbor1_index)+3] + r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                    # Compute block matrices
                    L_block = 0.5*(1+theta)*C0 - epsI*np.eye(3)
                    M_block = -2*(1+theta)*C0 + 2*epsI*np.eye(3)
                    N_block = 1.5*(1+theta)*C0 + (1-epsI)*np.eye(3)
                    f_block = B0invg + De

                    # Populate matrix
                    K[3*(curr_index):3*(curr_index)+3,3*(neighbor2_index):3*(neighbor2_index)+3] = L_block
                    K[3*(curr_index):3*(curr_index)+3,3*(neighbor1_index):3*(neighbor1_index)+3] = M_block
                    K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = N_block
                    f[3*(curr_index):3*(curr_index)+3] = f_block

                else:  # central
                    De = epsE*(r0[3*(neighbor1_index):3*(neighbor1_index)+3] - 2*r0[3*(curr_index):3*(curr_index)+3] + r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                    # Compute block matrices
                    L_block = -0.5*(1+theta)*C0 - epsI*np.eye(3)
                    M_block = (1 + 2*epsI)*np.eye(3)
                    N_block = 0.5*(1+theta)*C0 - epsI*np.eye(3)
                    f_block = B0invg + De

                    # Populate matrix
                    K[3*(curr_index):3*(curr_index)+3,3*(neighbor1_index):3*(neighbor1_index)+3] = L_block
                    K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = M_block
                    K[3*(curr_index):3*(curr_index)+3,3*(neighbor2_index):3*(neighbor2_index)+3] = N_block
                    f[3*(curr_index):3*(curr_index)+3] = f_block

        # Now loop over each node

        for index in [0]:

            if bc1 is 'splay':

                # Get coordinates

                r_curr = r0[3*(index):3*(index)+3]
                r_next = r0[3*(index+1):3*(index+1)+3]

                # Get vector that connects r_next to r_curr
                d_vec = r_next - r_curr

                # Get marching direction vector (orthogonal to the curve and to the surface normal)
                d_vec_rot = np.cross(N0[:,index],d_vec)

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.array([[d_vec_rot[0], d_vec_rot[1], d_vec_rot[2]],
                                                                [N0[0,index], N0[1,index], N0[2,index]],
                                                                [d_vec[0], d_vec[1], d_vec[2]]])
                f[3*index:3*index+3] = np.array([S0[index]*(1-sigmaSplay), 0, 0])


            elif bc1 is 'constX':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constY':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constZ':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1.lower().startswith('curve'):

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] = S0[index] * N0[:,index]

            else:

                # Call assembly routine
                matrixBuilder(index)

        for index in xrange(1,self.numNodes-1):
            # Call assembly routine
            matrixBuilder(index)

        for index in [self.numNodes-1]:

            if bc2 is 'continuous':

                # Populate matrix (use same displacements of first node)
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                K[3*index:3*index+3,:3] = -np.eye(3)
                f[3*index:3*index+3] = [0, 0, 0]

            elif bc2 is 'splay':

                # Get coordinates
                r_curr = r0[3*(index):3*(index)+3]
                r_prev = r0[3*(index-1):3*(index-1)+3]

                # Get vector that connects r_next to r_curr
                d_vec = r_curr - r_prev

                # Get marching direction vector (orthogonal to the curve and to the surface normal)
                d_vec_rot = np.cross(N0[:,index],d_vec)

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.array([[d_vec_rot[0], d_vec_rot[1], d_vec_rot[2]],
                                                                [N0[0,index], N0[1,index], N0[2,index]],
                                                                [d_vec[0], d_vec[1], d_vec[2]]])
                f[3*index:3*index+3] = np.array([S0[index]*(1-sigmaSplay), 0, 0])

            elif bc2 is 'constX':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc2 is 'constY':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc2 is 'constZ':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif self.optionsDict['bc2'].lower().startswith('curve'):

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
                f[3*index:3*index+3] = S0[index] * N0[:,index]

            else:

                # Call assembly routine
                matrixBuilder(index)

        # Solve the linear system
        dr = np.linalg.solve(K,f)

        # RETURNS
        return dr

    def smoothing(self, r, eta):

        alphaP0 = self.optionsDict['alphaP0']
        numSmoothingPasses = self.optionsDict['numSmoothingPasses']
        numLayers = self.optionsDict['numLayers']

        # This function does the grid smoothing

        # Loop over the desired number of smoothing passes
        for index_pass in range(numSmoothingPasses):

            # Initialize array of smoothed coordinates
            r_smooth = np.zeros(3*self.numNodes)

            # Copy the edge nodes
            r_smooth[:3] = r[:3]
            r_smooth[-3:] = r[-3:]

            # Smooth every node
            for index in xrange(1,self.numNodes-1):

                # Get coordinates
                r_curr = r[3*(index):3*(index)+3]
                r_next = r[3*(index+1):3*(index+1)+3]
                r_prev = r[3*(index-1):3*(index-1)+3]

                # Compute distances
                lp = np.linalg.norm(r_next - r_curr)
                lm = np.linalg.norm(r_curr - r_prev)

                # Compute alpha'
                alphaP = min(alphaP0, alphaP0*(eta-2)/numLayers)

                # Compute smoothed coordinates
                r_smooth[3*index:3*index+3] = (1-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)

            # Copy coordinates to allow next pass
            r = r_smooth[:]

        # RETURNS
        return r

    def dissipationCoefficients(self, layerIndex, r0_xi, r0_eta, dSensor, angle):

        # Get options
        numLayers = self.optionsDict['numLayers']
        epsE0 = self.optionsDict['epsE0']

        # Compute N (Eq. 6.3)
        N = np.linalg.norm(r0_eta)/np.linalg.norm(r0_xi)

        # Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layerIndex+2
        ltrans = int(3/4*numLayers)

        if l <= ltrans:
            Sl = np.sqrt((l-1)/(ltrans-1))
        else:
            Sl = np.sqrt((ltrans-1)/(numLayers-1))

        # Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max([dSensor**(2/Sl), 0.1])

        # Compute a (Eq 6.12 adjusted for entire angle (angle=2*alpha))
        if angle <= np.pi: # Convex corner
            a = 1.0
        else:
            a = 1.0/(1.0 - np.cos(angle/2)*np.cos(angle/2))

        # Compute auxiliary variable R (Eq. 6.4)
        R = Sl*dbar*a

        # Compute the dissipation coefficients
        epsE = epsE0*R*N
        epsI = 2*epsE

        # RETURNS
        return epsE, epsI

    def qualityCheck(self, R, layerIndex=None):
        '''
        This function checks the quality of a given mesh interval. This is done
        by constructing the Jacobians at each node, taking its determinant,
        and finding the ratio of the min over the max of each of the four
        nodal Jacobians for each face. Values near 1 within ratios are
        desirable while values near 0 point to bad mesh quality. Values below
        0 mean that the mesh is no longer valid and we stop marching.
        '''

        # Convert the flattened array R into a 3 x nw1 x nw2 array.
        # nw1 -> number of nodes in the direction of marching
        # nw2 -> number of nodes in direction of curve
        XYZ = np.array([R[:, ::3], R[:, 1::3], R[:, 2::3]])
        nw1, nw2 = XYZ.shape[1:]

        # Setup nodal normals
        nodalNormals = np.empty((3, nw1, nw2))

        # Get the panel normals from the interior points of the mesh.
        # Here we take the cross product of the diagonals of each face
        vec1 = XYZ[:, 1:, 1:] - XYZ[:, :-1, :-1]
        vec2 = XYZ[:, 1:, :-1] - XYZ[:, :-1, 1:]
        panelNormals = np.cross(vec2, vec1, axis=0)
        panelNormals = panelNormals / np.linalg.norm(panelNormals, axis=0)

        # Set the interior normals using an average of the panel normals
        vec1 = panelNormals[:, 1:, 1:] + panelNormals[:, :-1, :-1]
        vec2 = panelNormals[:, 1:, :-1] + panelNormals[:, :-1, 1:]
        normals = vec1 + vec2
        nodalNormals[:, 1:-1, 1:-1] = normals / np.linalg.norm(normals, axis=0)

        # Set the boundary normals
        nodalNormals[:, 1:, 0] = panelNormals[:, :, 0]
        nodalNormals[:, 0, :-1] = panelNormals[:, 0, :]
        nodalNormals[:, :-1, -1] = panelNormals[:, :, -1]
        nodalNormals[:, -1, 1:] = panelNormals[:, -1, :]

        # Setup nodal derivatives
        nodalDerivs = np.zeros((3, 2, nw1, nw2))

        # Compute interior derivatives using 2nd order central differencing
        nodalDerivs[:, 0, 1:-1, 1:-1] = (XYZ[:, 2:, 1:-1] - XYZ[:, :-2, 1:-1]) / 2.
        nodalDerivs[:, 1, 1:-1, 1:-1] = (XYZ[:, 1:-1, 2:] - XYZ[:, 1:-1, :-2]) / 2.

        # Compute i derivatives using 1st order differencing
        nodalDerivs[:, 0, 0, :] = XYZ[:, 1, :] - XYZ[:, 0, :]
        nodalDerivs[:, 0, -1, :] = XYZ[:, -1, :] - XYZ[:, -2, :]

        nodalDerivs[:, 0, 1:-1, 0] = (XYZ[:, 2:, 0] - XYZ[:, :-2, 0]) / 2.
        nodalDerivs[:, 0, 1:-1, -1] = (XYZ[:, 2:, -1] - XYZ[:, :-2, -1]) / 2.

        # Compute j derivatives using 1st order differencing
        nodalDerivs[:, 1, :, 0] = XYZ[:, :, 1] - XYZ[:, :, 0]
        nodalDerivs[:, 1, :, -1] = XYZ[:, :, -1] - XYZ[:, :, -2]

        nodalDerivs[:, 1, 0, 1:-1] = (XYZ[:, 0, 2:] - XYZ[:, 0, :-2]) / 2.
        nodalDerivs[:, 1, -1, 1:-1] = (XYZ[:, -1, 2:] - XYZ[:, -1, :-2]) / 2.

        # Assemble nodal Jacobians
        nodalJacs = np.zeros((3, 3, nw1, nw2))
        nodalJacs[0, :, :, :] = nodalDerivs[:, 0, :, :]
        nodalJacs[1, :, :, :] = nodalDerivs[:, 1, :, :]
        nodalJacs[2, :, :, :] = nodalNormals[:, :, :]

        # Compute determinants of Jacobians and find ratio of min to max per face
        ratios = np.zeros((nw1-1, nw2-1))

        # Compute the determinants of each nodal Jacobian
        det = np.linalg.det(np.swapaxes(np.swapaxes(nodalJacs, 1, 3), 0, 2))

        # Find the ratio of the minimum valued determinant to the maximum
        # valued determinant.
        # This is a measure of quality, with 1 being desirable and anything
        # less than 0 meaning the mesh is no longer valid.
        for i in range(nw1-1):
            for j in range(nw2-1):
                ratios[i, j] = np.min(det[i:i+2, j:j+2]) / np.max(det[i:i+2, j:j+2])

        fail = False
        # Throw an error and set the failure flag if the mesh is not valid
        if (np.any(np.isnan(ratios)) or np.min(ratios) <= 0.) and layerIndex:
            error("The mesh is not valid after step {}.".format(layerIndex+1))
            fail = True

        # Throw a warning if the mesh is low quality
        elif np.min(ratios) <= .2 and layerIndex:
            warn("The mesh may be low quality after step {}.".format(layerIndex+1))

        return fail, ratios

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
        myGrid.add_block(np.array([[self.curve[:,0]]]),
                         np.array([[self.curve[:,1]]]),
                         np.array([[self.curve[:,2]]]))

        # Export grid
        plot3d_interface.export_plot3d(myGrid, filename)

    def _getDefaultOptions(self):
        """ Define default options and pass back a dict. """

        self.optionsDict = {

            'bc1' : 'splay',
            'bc2' : 'splay',
            'dStart' : 1.e-2,
            'numLayers' : 17,
            'extension' : 2.,
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
            'growthRatio' : 1.2,
            }

    def _applyUserOptions(self, options):
        # Override default options with user options
        for userKeys in options.keys():
            unusedOption = True
            for defaultKeys in self.optionsDict.keys():
                if userKeys.lower() == defaultKeys.lower():
                    unusedOption = False
                    self.optionsDict[defaultKeys] = options[userKeys]
                    break
            if unusedOption:
                message = "{} key not in default options dictionary.".format(userKeys)
                warn(message)


if __name__ == '__main__':
    st = time()
    mesh = HypSurfMesh()
    mesh.createMesh()
    mesh.exportMesh('output.xyz')


'''
==============================================
MORE AUXILIARY FUNCTIONS
==============================================
'''

#=============================================
#=============================================

def giveAngle(r0,r1,r2,N1):

    '''
    This function gives the angle between the vectors joining
    r0 to r1 and r1 to r2. We assume that the body is to the right
    of the vector, while the mesh is propagated to the left

          r0
          |
     body | mesh
          |
          V angle
          r1--------->r2
              body

    N1 is the surface normal at the point r1

    angles > pi indicate convex corners, while angles < pi
    indicate concave corners
    '''

    dr1 = r1 - r0
    dr2 = r2 - r1

    dr1dotdr2 = dr1.dot(dr2) # dot product
    dr1crossdr2 = np.cross(dr1,dr2) # cross product

    # Compute acute angle and ensure it's <= 1.0
    arccos_inside = dr1dotdr2/np.linalg.norm(dr1)/np.linalg.norm(dr2)
    angle = np.arccos(np.min([arccos_inside, 1.0]))

    # If the cross product points in the same direction of the surface
    # normal, we have an acute corner
    if dr1crossdr2.dot(N1) > 0:
        angle = np.pi + angle
    else:
        angle = np.pi - angle

    return angle

#=============================================
#=============================================

def findRatio(dMax, d0, numLayers, ratioGuess):

    '''
    This function returns the geometrical progression ratio that satisfies
    the farfield distance and the number of cells. Newton search is used
    INPUTS
    dMax: distance that should be reached
    s0: cell edge length at the wall
    numLayers: number of cells used to reach farfield
    '''

    # Extra parameters
    nIters = 200 # Maximum number of iterations for Newton search
    q0 = ratioGuess # Initial ratio

    # Initialize ratio
    q = q0

    # Newton search loop
    for i in xrange(nIters):
       # Residual function
       R = d0*(1-q**(numLayers-1)) - dMax*(1-q)

       # Residual derivative
       Rdot = -(numLayers-1)*d0*q**(numLayers-2) + dMax

       # Update ratio with Newton search
       q = q - R/Rdot

    # Check if we got a reasonable value
    if (q <= 1) or (q >= q0):
        error('Ratio may be too large...\nIncrease number of cells or reduce extension')
        from sys import exit
        exit()

    # RETURNS
    return q

#=============================================
#=============================================

def findRadius(r):

    '''
    This function find the largest radius of the bounding box
    that encompasses the given mesh. This length will be used
    to calculate the maching distance
    '''

    # IMPORTS
    from numpy import max, min

    # Split coordinates
    x = r[::3]
    y = r[1::3]
    z = r[2::3]

    # Find bounds
    minX = np.min(x)
    maxX = np.max(x)
    minY = np.min(y)
    maxY = np.max(y)
    minZ = np.min(z)
    maxZ = np.max(z)

    # Find longest radius (we give only half of the largest side to be considered as radius)
    radius = np.max([maxX-minX, maxY-minY, maxZ-minZ])/2

    # RETURNS
    return radius

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
