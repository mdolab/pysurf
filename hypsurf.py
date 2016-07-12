# 3D Hyperbolic Surface Grid generator
# Written by:
# Ney Rafael Secco
# John Jasa
# neysecco@umich.edu
# July 2016

from __future__ import division

'''
TO DO

- add subiteration loop
- check smoothing
- set up scaling factor based on distance
- blend the angle-based dissipation coefficient
'''

def createMesh(rStart, bc1, bc2, dStart, numLayers, extension, surf,
               epsE0 = 1.0, theta = 0.0,
               alphaP0 = 0.25, numSmoothingPasses = 0,
               nuArea = 0.16, numAreaPasses = 0,
               sigmaSplay = 0.3,
               cMax = 3.0,
               ratioGuess = 20):

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

    # IMPORTS
    from numpy import zeros, copy, ceil

    # Flatten the coordinates vector
    # We do this because the linear system is assembled assuming that
    # the coordinates vector is flattened
    rStart = rStart.flatten()

    # Get the number of nodes
    numNodes = int(len(rStart)/3)

    # Initialize 2D array that will contains all the surface grid points in the end
    R = zeros((numLayers,len(rStart)))

    # Store the initial curve
    R[0,:] = rStart

    # Initialize step size and total marched distance
    d = dStart
    dTot = 0

    # Find the characteristic radius of the mesh
    radius = findRadius(rStart)

    # Find the desired marching distance
    dMax = radius*(extension-1)

    # Compute the growth ratio necessary to match this distance
    dGrowth = findRatio(dMax, dStart, numLayers, ratioGuess)

    # Print growth ratio
    print 'Growth ratio: ',dGrowth

    # Project the given points to the surface
    # (just to make sure the points are actually one the surface)
    surf.inverse_evaluate(rStart.reshape((numNodes, 3)))
    rNext = surf.get_points(None).flatten()

    # Compute surface normals at the projected points that will be used in the first iteration
    NNext = surf.get_normals(None).T

    # Issue a warning message if the projected points are far from the given points
    if max(abs(rNext - rStart)) > 1.0e-4:
        print ''
        print 'WARNING!:' 
        print 'The given points (rStart) might not belong to the given surface (surf)'
        print 'These points were projected for the surface mesh generation'
        print ''

    # We need a guess for the first-before-last curve in order to compute the grid distribution sensor
    # As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
    rm1 = rStart[:]

    '''
    ==========================================
    Now let's define the functions that do the heavy work.
    These functions are the actual hyperbolic generation
    calculations.
    I define them here so that the several optional parameters
    are readily available to them.
    ==========================================
    '''

    def subIteration(r0, N0, S0, rm1, Sm1, layerIndex):

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

        # IMPORTS
        from numpy.linalg import solve

        # Generate matrices of the linear system
        K,f = computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex)

        # Solve the linear system
        dr = solve(K,f)

        # Get maximum residual
        res = K.dot(dr) - f
        maxRes = max(abs(res))

        # Update r
        rNext = r0 + dr

        # Smooth coordinates
        rNext = smoothing(rNext,layerIndex+2)

        # Project onto surface and compute surface normals
        surf.inverse_evaluate(rNext.reshape((numNodes, 3)))
        rNext = surf.get_points(None).flatten()
        NNext = surf.get_normals(None).T

        # RETURNS
        return rNext, NNext, maxRes

    def areaFactor(r0,d):

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

        # IMPORT
        from numpy import hstack, reshape, max
        from numpy.linalg import norm

        # Extrapolate the end points
        r0minus1 = 2*r0[0:3] - r0[3:6]
        r0plus1 = 2*r0[-3:] - r0[-6:-3]
        r0_extrap = hstack([r0minus1, r0, r0plus1])

        # Reshape so we have 3D array
        R0_extrap = reshape(r0_extrap,(-1,3)).T

        # Compute the distance of each node to its neighbors
        neighborDist = 0.5*(norm(R0_extrap[:,1:-1] - R0_extrap[:,:-2],axis=0) + norm(R0_extrap[:,2:] - R0_extrap[:,1:-1],axis=0))

        # Multiply distances by the step size to get the areas
        S = d*neighborDist

        # Divide the marching distance and the neighbor distance to get the stretch ratios
        stretchRatio = d/neighborDist

        # Get the maximum stretch ratio
        maxStretch = max(stretchRatio)

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

        # RETURNS
        return S, maxStretch

    def computeMatrices(r0,N0,S0,rm1,Sm1,layerIndex):

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

        # IMPORTS
        from numpy import zeros, sqrt, eye, array, arctan2, pi, cross
        from numpy.linalg import norm, inv

        # Find number of nodes
        numNodes = int(len(r0)/3)

        # Initialize arrays
        K = zeros((3*numNodes,3*numNodes))
        f = zeros(3*numNodes)

        # Define assembly routines
        def centralMatrices(prev_index,curr_index,next_index):

            # Using central differencing for zeta = 2:numNodes-1
            r0_xi = 0.5*(r0[3*(next_index):3*(next_index)+3] - r0[3*(prev_index):3*(prev_index)+3])
            x0_xi = r0_xi[0]
            y0_xi = r0_xi[1]
            z0_xi = r0_xi[2]

            # Get current normal
            nx = N0[0,curr_index]
            ny = N0[1,curr_index]
            nz = N0[2,curr_index]

            # Assemble B0 matrix
            B0 = array([[x0_xi, y0_xi, z0_xi],
                        [ny*z0_xi-nz*y0_xi, nz*x0_xi-nx*z0_xi, nx*y0_xi-ny*x0_xi],
                        [nx, ny, nz]])

            # Invert B0
            B0inv = inv(B0)

            # Compute eta derivatives
            r0_eta = B0inv.dot(array([0, Sm1[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(next_index):3*(next_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(prev_index):3*(prev_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(next_index):3*(next_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(prev_index):3*(prev_index)+3]-r0[3*(index):3*(index)+3])
            dSensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = giveAngle(r0[3*(prev_index):3*(prev_index)+3],
                              r0[3*(curr_index):3*(curr_index)+3],
                              r0[3*(next_index):3*(next_index)+3],
                              N0[:,curr_index])

            # Sharp convex corner detection
            if angle > 240*pi/180: # Corner detected

                # Populate matrix with Eq 8.3
                K[3*curr_index:3*curr_index+3,3*next_index:3*next_index+3] = -eye(3)
                K[3*curr_index:3*curr_index+3,3*curr_index:3*curr_index+3] = 2*eye(3)
                K[3*curr_index:3*curr_index+3,3*prev_index:3*prev_index+3] = -eye(3)
                f[3*curr_index:3*curr_index+3] = array([0,0,0])

            else:

                # Compute C0 = B0inv*A0
                C0 = B0inv.dot(A0)

                # Compute smoothing coefficients
                epsE, epsI = dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle)

                # Compute RHS components
                B0invg = B0inv.dot(array([0, S0[curr_index], 0]))
                De = epsE*(r0[3*(prev_index):3*(prev_index)+3] - 2*r0[3*(index):3*(index)+3] + r0[3*(next_index):3*(next_index)+3])

                # Compute block matrices
                L_block = -0.5*(1+theta)*C0 - epsI*eye(3)
                M_block = (1 + 2*epsI)*eye(3)
                N_block = 0.5*(1+theta)*C0 - epsI*eye(3)
                f_block = B0invg + De

                # Populate matrix
                K[3*(curr_index):3*(curr_index)+3,3*(prev_index):3*(prev_index)+3] = L_block
                K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = M_block
                K[3*(curr_index):3*(curr_index)+3,3*(next_index):3*(next_index)+3] = N_block
                f[3*(curr_index):3*(curr_index)+3] = f_block

        def forwardMatrices(curr_index,next_index,nextnext_index):

            # Using forward differencing for xi = 1
            r0_xi = 0.5*(-3*r0[3*(curr_index):3*(curr_index)+3] + 4*r0[3*(next_index):3*(next_index)+3] - r0[3*(nextnext_index):3*(nextnext_index)+3])
            x0_xi = r0_xi[0]
            y0_xi = r0_xi[1]
            z0_xi = r0_xi[2]

            # Get current normal
            nx = N0[0,curr_index]
            ny = N0[1,curr_index]
            nz = N0[2,curr_index]

            # Assemble B0 matrix
            B0 = array([[x0_xi, y0_xi, z0_xi],
                        [ny*z0_xi-nz*y0_xi, nz*x0_xi-nx*z0_xi, nx*y0_xi-ny*x0_xi],
                        [nx, ny, nz]])

            # Invert B0
            B0inv = inv(B0)

            # Compute eta derivatives
            r0_eta = B0inv.dot(array([0, Sm1[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(next_index):3*(next_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(nextnext_index):3*(nextnext_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(next_index):3*(next_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(nextnext_index):3*(nextnext_index)+3]-r0[3*(index):3*(index)+3])
            dSensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

            # Compute C0 = B0inv*A0
            C0 = B0inv.dot(A0)

            # Compute smoothing coefficients
            epsE, epsI = dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle)

            # Compute RHS components
            B0invg = B0inv.dot(array([0, S0[curr_index], 0]))
            De = epsE*(r0[3*(nextnext_index):3*(nextnext_index)+3] - 2*r0[3*(next_index):3*(next_index)+3] + r0[3*(curr_index):3*(curr_index)+3])

            # Compute block matrices
            L_block = -0.5*(1+theta)*C0 - epsI*eye(3)
            M_block = 2*(1+theta)*C0 + 2*epsI*eye(3)
            N_block = -1.5*(1+theta)*C0 + (1-epsI)*eye(3)
            f_block = B0invg + De

            # Populate matrix
            K[3*(curr_index):3*(curr_index)+3,3*(nextnext_index):3*(nextnext_index)+3] = L_block
            K[3*(curr_index):3*(curr_index)+3,3*(next_index):3*(next_index)+3] = M_block
            K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = N_block
            f[3*(curr_index):3*(curr_index)+3] = f_block

        def backwardMatrices(curr_index,prev_index,prevprev_index):

            # Using backward differencing for xi = numNodes
            r0_xi = 0.5*(3*r0[3*(curr_index):3*(curr_index)+3] - 4*r0[3*(prev_index):3*(prev_index)+3] + r0[3*(prevprev_index):3*(prevprev_index)+3])
            x0_xi = r0_xi[0]
            y0_xi = r0_xi[1]
            z0_xi = r0_xi[2]

            # Get current normal
            nx = N0[0,curr_index]
            ny = N0[1,curr_index]
            nz = N0[2,curr_index]

            # Assemble B0 matrix
            B0 = array([[x0_xi, y0_xi, z0_xi],
                        [ny*z0_xi-nz*y0_xi, nz*x0_xi-nx*z0_xi, nx*y0_xi-ny*x0_xi],
                        [nx, ny, nz]])

            # Invert B0
            B0inv = inv(B0)

            # Compute eta derivatives
            r0_eta = B0inv.dot(array([0, Sm1[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(prev_index):3*(prev_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(prevprev_index):3*(prevprev_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(prev_index):3*(prev_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(prevprev_index):3*(prevprev_index)+3]-r0[3*(index):3*(index)+3])
            dSensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

             # Compute C0 = B0inv*A0
            C0 = B0inv.dot(A0)

            # Compute smoothing coefficients
            epsE, epsI = dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle)

            # Compute RHS components
            B0invg = B0inv.dot(array([0, S0[curr_index], 0]))
            De = epsE*(r0[3*(prevprev_index):3*(prevprev_index)+3] - 2*r0[3*(prev_index):3*(prev_index)+3] + r0[3*(curr_index):3*(curr_index)+3])

            # Compute block matrices
            L_block = 0.5*(1+theta)*C0 - epsI*eye(3)
            M_block = -2*(1+theta)*C0 + 2*epsI*eye(3)
            N_block = 1.5*(1+theta)*C0 + (1-epsI)*eye(3)
            f_block = B0invg + De

            # Populate matrix
            K[3*(curr_index):3*(curr_index)+3,3*(prevprev_index):3*(prevprev_index)+3] = L_block
            K[3*(curr_index):3*(curr_index)+3,3*(prev_index):3*(prev_index)+3] = M_block
            K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = N_block
            f[3*(curr_index):3*(curr_index)+3] = f_block

        # Now loop over each node

        for index in [0]:

            if bc1 is 'free':

                # Name indexes
                curr_index = index
                next_index = index+1
                nextnext_index = index+2

                # Call assembly routine
                forwardMatrices(curr_index,next_index,nextnext_index)

            elif bc1 is 'continuous':

                # Name indexes
                prev_index = numNodes-2
                curr_index = index
                next_index = index+1

                # Call assembly routine
                centralMatrices(prev_index,curr_index,next_index)

            elif bc1 is 'splay':

                # Get coordinates
                r_curr = r0[3*(index):3*(index)+3]
                r_next = r0[3*(index+1):3*(index+1)+3]

                # Get vector that connects r_next to r_curr
                d_vec = r_next - r_curr

                # Get marching direction vector (orthogonal to the curve and to the surface normal)
                d_vec_rot = cross(N0[:,index],d_vec)

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = array([[d_vec_rot[0], d_vec_rot[1], d_vec_rot[2]],
                                                                [N0[0,index], N0[1,index], N0[2,index]],
                                                                [d_vec[0], d_vec[1], d_vec[2]]])
                f[3*index:3*index+3] = array([S0[index]*(1-sigmaSplay), 0, 0])

            elif bc1 is 'constX':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constY':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constZ':

                # Populate matrix
                K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

        for index in xrange(1,numNodes-1):

            # Name indexes
            prev_index = index-1
            curr_index = index
            next_index = index+1

            # Call assembly routine
            centralMatrices(prev_index,curr_index,next_index)

        for index in [numNodes-1]:

            if bc2 is 'free':

                # Name indexes
                curr_index = index
                prev_index = index-1
                prevprev_index = index-2

                # Call assembly routine
                backwardMatrices(curr_index,prev_index,prevprev_index)

            elif bc2 is 'continuous':

                # Populate matrix (use same displacements of first node)
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                K[3*index:3*index+3,:3] = -eye(3)
                f[3*index:3*index+3] = [0, 0, 0]

            elif bc2 is 'splay':

                # Get coordinates
                r_curr = r0[3*(index):3*(index)+3]
                r_prev = r0[3*(index-1):3*(index-1)+3]

                # Get vector that connects r_next to r_curr
                d_vec = r_curr - r_prev

                # Get marching direction vector (orthogonal to the curve and to the surface normal)
                d_vec_rot = cross(N0[:,index],d_vec)

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = array([[d_vec_rot[0], d_vec_rot[1], d_vec_rot[2]],
                                                                [N0[0,index], N0[1,index], N0[2,index]],
                                                                [d_vec[0], d_vec[1], d_vec[2]]])
                f[3*index:3*index+3] = array([S0[index]*(1-sigmaSplay), 0, 0])

            elif bc2 is 'constX':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc2 is 'constY':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc2 is 'constZ':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

        # RETURNS
        return K,f

    def smoothing(r,eta):

        # This function does the grid smoothing

        # IMPORTS
        from numpy import copy, zeros
        from numpy.linalg import norm

        # Loop over the desired number of smoothing passes
        for index_pass in range(numSmoothingPasses):

            # Initialize array of smoothed coordinates
            r_smooth = zeros(3*numNodes)

            # Copy the edge nodes
            r_smooth[:3] = r[:3]
            r_smooth[-3:] = r[-3:]

            # Smooth every node
            for index in xrange(1,numNodes-1):

                # Get coordinates
                r_curr = r[3*(index):3*(index)+3]
                r_next = r[3*(index+1):3*(index+1)+3]
                r_prev = r[3*(index-1):3*(index-1)+3]

                # Compute distances
                lp = norm(r_next - r_curr)
                lm = norm(r_curr - r_prev)

                # Compute alpha'
                alphaP = min(alphaP0, alphaP0*(eta-2)/numLayers)

                # Compute smoothed coordinates
                r_smooth[3*index:3*index+3] = (1-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)

            # Copy coordinates to allow next pass
            r = copy(r_smooth)

        # RETURNS
        return r

    def dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle):

        # IMPORTS
        from numpy import sqrt, cos, pi
        from numpy.linalg import norm

        # Compute N (Eq. 6.3)
        N = norm(r0_eta)/norm(r0_xi)

        # Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layerIndex+2
        ltrans = int(3/4*numLayers)
        if l <= ltrans:
            Sl = sqrt((l-1)/(ltrans-1))
        else:
            Sl = sqrt((ltrans-1)/(numLayers-1))

        # Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max([dSensor**(2/Sl), 0.1])

        # Compute a (Eq 6.12 adjusted for entire angle (angle=2*alpha))
        if angle >= pi: # Convex corner
            a = 1.0
        else:
            a = 1.0/(1.0 - cos(angle/2)*cos(angle/2))

        # Compute auxiliary variable R (Eq. 6.4)
        R = Sl*dbar*a

        # Compute the dissipation coefficients
        epsE = epsE0*R*N
        epsI = 2*epsE

        # RETURNS
        return epsE, epsI

    #===========================================================

    '''
    The marching function actually begins here
    '''
    # Some functions require the area factors of the first-before-last curve
    # We will repeat the first curve areas for simplicity.
    # rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
    # But we still need to find Sm1
    Sm1, maxStretch = areaFactor(rStart,d)

    maxRes = 0

    # MARCH!!!
    for layerIndex in range(numLayers-1):

        # Get the coordinates computed by the previous iteration
        r0 = rNext[:]

        # Compute the new area factor for the desired marching distance
        S0, maxStretch = areaFactor(r0,d)

        # The subiterations will use pseudo marching steps.
        # If the required marching step is too large, the hyperbolic marching might become
        # unstable. So we subdivide the current marching step into multiple smaller pseudo-steps

        # Compute the factor between the current stretching ratio and the allowed one.
        # If the current stretching ratio is smaller than cMax, the cFactor will be 1.0, and
        # The pseudo-step will be the same as the desired step.
        cFactor = int(ceil(maxStretch/cMax))

        # Constrain the marching distance if the stretching ratio is too high
        dPseudo = d/cFactor

        print cFactor, maxRes

        # Subiteration
        # The number of subiterations is the one required to meet the desired marching distance
        for indexSubIter in range(cFactor):

            # Recompute areas with the pseudo-step
            S0, maxStretch = areaFactor(r0,dPseudo)

            # Update the Normals with the values computed in the last iteration.
            # We do this because, the last iteration already projected the new
            # points to the surface and also computed the normals. So we don't
            # have to repeat the projection step
            N0 = NNext[:,:]

            # March using the pseudo-marching distance
            rNext, NNext, maxRes = subIteration(r0, N0, S0, rm1, Sm1, layerIndex)

            # Update Sm1 (Store the previous area factors)
            Sm1 = S0[:]

            # Update rm1
            rm1 = r0[:]

            # Update r0
            r0 = rNext[:]

        # Store grid points
        R[layerIndex+1,:] = rNext

        # Compute the total marched distance so far
        dTot = dTot + d

        # Update step size
        d = d*dGrowth

    # Convert to X, Y and Z
    X = R[:,::3]
    Y = R[:,1::3]
    Z = R[:,2::3]

    # RETURNS
    return X,Y,Z

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

    # IMPORTS
    from numpy import pi, cross, arccos
    from numpy.linalg import norm

    dr1 = r1 - r0
    dr2 = r2 - r1

    dr1dotdr2 = dr1.dot(dr2) # dot product
    dr1crossdr2 = cross(dr1,dr2) # cross product

    # Compute acute angle
    angle = arccos(dr1dotdr2/norm(dr1)/norm(dr2))

    # If the cross product points in the same direction of the surface
    # normal, we have an acute corner
    if dr1crossdr2.dot(N1) > 0:
        angle = pi + angle
    else:
        angle = pi - angle

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
       print 'Ratio may be too large...'
       print 'Increase number of cells or reduce extension'
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
    minX = min(x)
    maxX = max(x)
    minY = min(y)
    maxY = max(y)
    minZ = min(z)
    maxZ = max(z)

    # Find longest radius (we give only half of the largest side to be considered as radius)
    radius = max([maxX-minX, maxY-minY, maxZ-minZ])/2

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

#=============================================
#=============================================

def exportPlot3d(X,Y,Z,filename):

    '''
    This function exports a 3D mesh in plot3d format.
    The user specifies the span and number of nodes in the z direction.
    '''

    # IMPORTS
    from plot3d_interface import Grid, export_plot3d
    from numpy import array, copy, ones, linspace

    # Initialize grid object
    myGrid = Grid()

    # Expand the coordinate matrices
    X3d = array([X])
    Y3d = array([Y])
    Z3d = array([Z])

    # Add block to the grid
    myGrid.add_block(X3d, Y3d, Z3d)

    # Export grid
    export_plot3d(myGrid, filename)
