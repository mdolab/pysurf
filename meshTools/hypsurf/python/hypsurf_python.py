'''
This script contains a Python version of the surface mesh generation.
This version cannot compute derivatives, but we still keep it because
it is easier to test new features in a Python enviroment and then
pass it to Fortran later on.

Ney Secco 2017-01
'''

from __future__ import division
import numpy as np

def march(projection_func, rStart, dStart, theta, sigmaSplay, bc1, bc2,
          epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax,
          extensionGiven, guideIndices, remeshFlag, numSmoothingPasses,
          numAreaPasses, numLayers):

    '''
    This is the main function that executes the marching.
    I kept all these multiple inputs so that the Python and Fortran versions had the
    same interface.
    '''

    # INITIALIZATION

    # Make sure BCs are in lowercase
    bc1 = bc1.lower()
    bc2 = bc2.lower()

    # Get number of nodes
    numNodes = int(len(rStart)/3)

    # Create dictionary to make data passing easier
    extraData = {}
    extraData['theta'] = theta
    extraData['bc1'] = bc1
    extraData['bc2'] = bc2
    extraData['sigmaSplay'] = sigmaSplay
    extraData['numLayers'] = numLayers
    extraData['epsE0'] = epsE0
    extraData['theta'] = theta
    extraData['alphaP0'] = alphaP0
    extraData['numSmoothingPasses'] = numSmoothingPasses
    extraData['numNodes'] = numNodes
    extraData['remeshFlag'] = remeshFlag

    # COMPUTE NORMALIZED ARC-LENGTHS
    if remeshFlag:

        # If we detect guide curves, we will record separate arc-lengths for
        # each subinterval defined by these guide curves.

        # Initialize list of arc-lengths.
        # Each entry of this list will correspond to an arc-length subinterval
        # defined by the guide curves.
        arcLength = []

        # Now we define a list of nodes that define the boundaries of the subintervals.
        # These are the nodes that follow guide curves (and also the first and last nodes)
        breakNodes = [0] + guideIndices.tolist() + [numNodes-1]
        breakNodes = sorted(set(breakNodes)) # Remove duplicates

        # Compute the number of subinterval based on the number of break nodes
        numIntervals = len(breakNodes) - 1

        # Store the normalized arc-lengths of each subinterval
        for intervalID in range(numIntervals):

            # Get indices of the first and last nodes of the interval
            node1 = breakNodes[intervalID]
            node2 = breakNodes[intervalID+1]

            # Take the slice of the nodal coordinates that corresponds to the
            # current subinterval
            rInterval = rStart[3*node1:3*(node2+1)]

            # Compute the normalized arc-lengths of this interval
            currArcLength = compute_arc_length(rInterval)

            # Store the arc-lengths of the current interval
            arcLength.append(currArcLength)

    else:
        arcLength = []

    # Initialize 2D array that will contain all surface grid points in the end
    R = np.zeros((numLayers,len(rStart)))

    # Project initial curve onto the surface or curve (if applicable)
    rNext, NNext = projection_func(rStart)

    # Initialize step size and total marched distance
    d = dStart
    dTot = 0

    # Find the characteristic radius of the mesh
    radius = findRadius(rNext)

    if extensionGiven:

        # Save extension
        extension = marchParameter

        # Find the desired marching distance
        dMax = radius*(extension-1.0)

        # Compute the growth ratio necessary to match this distance
        dGrowth = findRatio(dMax, dStart, numLayers, ratioGuess)
    else:
        # Remember that growthRatio = self.marchParameter
        dGrowth = marchParameter

        # Recompute desired marching distance
        dMax = dStart*(1.0 - dGrowth**(numLayers-1))/(1.0 - dGrowth)

        # Recompute equivalent extension with the given growth ratio
        extension = dMax/radius + 1.0

    # Print marching information
    if extensionGiven:
        parameterName = 'extension'
    else:
        parameterName = 'growth ratio'
    print ''
    print 'Running with user provided '+parameterName
    print 'Extension: ',extension
    print 'Growth ratio: ',dGrowth
    print ''

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

    fail = False

    # MARCH!!!
    for layerIndex in range(numLayers-1):

        # Get the coordinates computed by the previous iteration
        r0 = rNext[:]
        N0 = NNext[:,:]

        # Compute the new area factor for the desired marching distance
        S0, maxStretch = areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, guideIndices)

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
            Sm1, maxStretch = areaFactor(rm1, dPseudo, nuArea, numAreaPasses, bc1, bc2, guideIndices)
            S0, maxStretch = areaFactor(r0, dPseudo, nuArea, numAreaPasses, bc1, bc2, guideIndices)

            # March using the pseudo-marching distance
            rNext, NNext = subIteration(r0, N0, S0, rm1, Sm1, layerIndex, guideIndices, arcLength, extraData, projection_func)

            # Update Sm1 (Store the previous area factors)
            Sm1 = S0[:]

            # Update rm1
            rm1 = r0[:]

            # Update r0
            r0 = rNext[:]

            # Update the Normals with the values computed in the last iteration.
            # We do this because, the last iteration already projected the new
            # points to the surface and also computed the normals. So we don't
            # have to repeat the projection step
            N0 = NNext[:,:]

        # Store grid points
        R[layerIndex+1,:] = rNext

        # Check quality of the mesh
        if layerIndex > 1:
            fail, ratios = qualityCheck(R[layerIndex-2:layerIndex+2, :], layerIndex)


        if fail:
            # If the mesh is not valid, only save the mesh up until that point.
            # Otherwise we'd see a bunch of points at 0,0,0.
            R = R[:layerIndex+2, :]
            break

        # Compute the total marched distance so far
        dTot = dTot + d

        # Update step size
        d = d*dGrowth

    # RETURNS
    return R, fail, ratios

#========================================================
# AUXILIARY FUNCTIONS

def subIteration(r0, N0, S0, rm1, Sm1, layerIndex, guideIndices, arcLength, extraData, projection_func):

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

    theta = extraData['theta']
    bc1 = extraData['bc1']
    bc2 = extraData['bc2']
    sigmaSplay = extraData['sigmaSplay']
    numLayers = extraData['numLayers']
    epsE0 = extraData['epsE0']
    theta = extraData['theta']
    alphaP0 = extraData['alphaP0']
    numSmoothingPasses = extraData['numSmoothingPasses']
    numNodes = extraData['numNodes']
    remeshFlag = extraData['remeshFlag']
    eta = layerIndex+2

    dr = computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, guideIndices, extraData)

    # Update r
    rNext = r0 + dr

    # Smooth coordinates
    rNext_ = smoothing(rNext, eta, alphaP0, numSmoothingPasses, numLayers, numNodes)

    rNext, NNext = projection_func(rNext_)

    # Remesh curve with initial spacing if chosen by the user
    if remeshFlag:

        # Now we will redistribute the nodes based on the original curve arc-length.
        # Remember that we should do this for each subinterval defined by guide curves.

        # Initialize array for redistributed coordinates
        rNext_ = np.zeros(rNext.shape)

        # Initialize offset index for the initial node of the subinterval
        nodeID_offset = 0

        # Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
        for intervalID in range(len(arcLength)):

            # Get the number of nodes in the current subinterval
            numNodes = len(arcLength[intervalID])

            # Slice the coordinate array to get only nodes of the current subinterval
            rInterval = rNext[nodeID_offset*3:(nodeID_offset+numNodes)*3]

            # Redistribute the nodes
            rNew = redistribute_nodes_by_arc_length(rInterval,arcLength[intervalID])

            # Save the redistributed nodes in the original array
            rNext_[nodeID_offset*3:(nodeID_offset+numNodes)*3] = rNew

            # Update the offset
            nodeID_offset = nodeID_offset + numNodes - 1

        # Project the remeshed curve back onto the surface
        rNext, NNext = projection_func(rNext_)

    # RETURNS
    return rNext, NNext

#========================================================

def areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, guideIndices):

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
    if bc1.startswith('curve'):
        S[0] = d
    if bc2.startswith('curve'):
        S[-1] = d
    # Set guideCurve marching distances
    if guideIndices.any():
        for index in guideIndices:
            S[index] = d

    # RETURNS
    return S, maxStretch

#========================================================

def computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, guideIndices, extraData):

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
    theta = extraData['theta']
    bc1 = extraData['bc1']
    bc2 = extraData['bc2']
    sigmaSplay = extraData['sigmaSplay']
    numLayers = extraData['numLayers']
    epsE0 = extraData['epsE0']
    numNodes = extraData['numNodes']

    # Initialize arrays
    K = np.zeros((3*numNodes, 3*numNodes))
    f = np.zeros(3*numNodes)

    #####################################
    # HELPER FUNCTION TO BUILD MATRICES #
    #####################################

    def matrixBuilder(curr_index):

        if curr_index == 0:  # forward case

            if bc1 != 'continuous':

                neighbor1_index = 1
                neighbor2_index = 2

                # Using forward differencing for xi = 1
                r0_xi = 0.5*(-3*r0[3*(curr_index):3*(curr_index)+3] + 4*r0[3*(neighbor1_index):3*(neighbor1_index)+3] - r0[3*(neighbor2_index):3*(neighbor2_index)+3])

                angle = np.pi

            else:

                neighbor1_index = numNodes - 2
                neighbor2_index = curr_index + 1

                # Using central differencing for zeta = 2:numNodes-1
                r0_xi = 0.5*(r0[3*(neighbor2_index):3*(neighbor2_index)+3] - r0[3*(neighbor1_index):3*(neighbor1_index)+3])

                # Compute the local grid angle based on the neighbors
                angle = giveAngle(r0[3*(neighbor1_index):3*(neighbor1_index)+3],
                                  r0[3*(curr_index):3*(curr_index)+3],
                                  r0[3*(neighbor2_index):3*(neighbor2_index)+3],
                                  N0[:,curr_index])


        elif curr_index == numNodes - 1:  # backward case
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
            epsE, epsI = dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0)

            # Compute RHS components
            B0invg = B0inv.dot(np.array([0, S0[curr_index], 0]))

            if curr_index == 0:

                if bc1 != 'continuous':# forwards
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

            elif curr_index == numNodes - 1:  # backwards
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

    #####################################
    # END OF HELPER FUNCTION            #
    #####################################

    # Now loop over each node

    for index in [0]:

        if bc1 == 'splay':

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


        elif bc1 == 'constX':

            # Populate matrix
            K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc1 == 'constY':

            # Populate matrix
            K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc1 == 'constZ':

            # Populate matrix
            K[3*index:3*index+3,3*(index+1):3*(index+1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc1.startswith('curve'):

            # Populate matrix
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] = S0[index] * N0[:,index]

        else:

            # Call assembly routine
            matrixBuilder(index)

    for index in xrange(1,numNodes-1):

        # Set guideCurve matrix contributions
        if guideIndices.any():
            if index in guideIndices:

                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)

                f[3*index:3*index+3] = S0[index] * N0[:,index]

                '''
                # Neighbor average
                # Populate matrix
                K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)/2
                K[3*index:3*index+3,3*index-3:3*index-3+3] = np.eye(3)/4
                K[3*index:3*index+3,3*index+3:3*index+3+3] = np.eye(3)/4

                f[3*index:3*index+3] = S0[index] * N0[:,index]
                '''

                '''
                # Dissipation coefficients

                # Call assembly routine
                matrixBuilder(index)

                # Get dissipation factor from original matrix
                epsI_local = (K[3*index,3*index] - 1.0)/2.0
                epsE_local = epsI_local/2.0

                # Here we keep the same diagonal matrix, but we modify
                # the neighbors' matrices to keep only dissipation terms
                K[3*index:3*index+3,3*index-3:3*index-3+3] = -epsI_local*np.eye(3)
                K[3*index:3*index+3,3*index+3:3*index+3+3] = -epsI_local*np.eye(3)

                # Compute dissipation term for RHS
                r_prev = r0[3*(index-1):3*(index-1)+3]
                r_curr = r0[3*(index):3*(index)+3]
                r_next = r0[3*(index+1):3*(index+1)+3]
                De = epsE_local*(r_prev - 2*r_curr + r_next)

                # Replace RHS to keep marching direction and dissipation term
                f[3*index:3*index+3] = S0[index] * N0[:,index] + De
                '''
            else:
                # Call assembly routine
                matrixBuilder(index)
        else:
            # Call assembly routine
            matrixBuilder(index)

    for index in [numNodes-1]:

        if bc2 == 'continuous':

            # Populate matrix (use same displacements of first node)
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            K[3*index:3*index+3,:3] = -np.eye(3)
            f[3*index:3*index+3] = [0, 0, 0]

        elif bc2 == 'splay':

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

        elif bc2 == 'constX':

            # Populate matrix
            K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc2 == 'constY':

            # Populate matrix
            K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc2 == 'constZ':

            # Populate matrix
            K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
            K[3*index:3*index+3,3*index:3*index+3] = np.eye(3)
            f[3*index:3*index+3] =  [0, 0, 0]

        elif bc2.startswith('curve'):

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

#========================================================

def smoothing(r, eta, alphaP0, numSmoothingPasses, numLayers, numNodes):

    '''
    This function does the grid smoothing
    '''

    # Loop over the desired number of smoothing passes
    for index_pass in range(numSmoothingPasses):

        # Initialize array of smoothed coordinates
        r_smooth = np.zeros(3*numNodes)

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

#========================================================

def dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0):

    # Compute N (Eq. 6.3)
    N = np.linalg.norm(r0_eta)/np.linalg.norm(r0_xi)

    # Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
    l = layerIndex+2
    ltrans = int(3/4*numLayers)

    if l <= ltrans:
        Sl = np.sqrt((l-1)/(numLayers-1))
    else:
        Sl = np.sqrt((ltrans-1)/(numLayers-1))

    # Compute adjusted grid distribution sensor (Eq. 6.7)
    dbar = max([dSensor**(2/Sl), 0.1])

    # Compute a (Eq 6.12 adjusted for entire angle (angle=2*alpha))
    if angle <= np.pi: # Convex corner
        a = 1.0
    else:
        a = 1.0/(1.0 - np.cos(angle/2)**2)

    # Compute auxiliary variable R (Eq. 6.4)
    R = Sl*dbar*a

    # Compute the dissipation coefficients
    epsE = epsE0*R*N
    epsI = 2*epsE

    # RETURNS
    return epsE, epsI

#========================================================

def qualityCheck(R, layerIndex=None):
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

#========================================================

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

#========================================================

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

def findRadius(r):

    '''
    This function find the largest radius of the bounding box
    that encompasses the given mesh. This length will be used
    to calculate the maching distance
    '''

    numNodes = int(len(r)/3)

    # Find the average of all nodes
    # We know that this is not the proper centroid, but it is ok since
    # the radius is already an approximation
    r_centroid = np.zeros(3)

    r_centroid[0] = np.sum(r[::3])/numNodes
    r_centroid[1] = np.sum(r[1::3])/numNodes
    r_centroid[2] = np.sum(r[2::3])/numNodes

    # Now we need a way to estimate the maximum distance to this centroid.
    # We will use a KS-function approach to estimate the maximum radius
    # Remember that KS = (1/rho)*log(sum(exp(rho*g_i)))
    # In our case, g_i are the distances of each node to the centroid.
    # According to the KS function property: KS >= max(g_i) for rho > 0,
    # So it can give an estimate of maximum radius

    # Define KS-function constant (if you change rho, remember to run Makefile_tapenade)
    rho = 60.0

    # Compute the sum of exponentials of the distances
    sumExp = 0.0

    for ii in range(numNodes):
           
        # Get coordinates of the current node
        r_node = r[3*ii:3*ii+3]

        # Compute distance
        distVec = r_node - r_centroid
        dist = np.sqrt(np.sum(distVec**2))

        # Add contribution to the sum
        sumExp = sumExp + np.exp(rho*dist)

    # Use the KS function to estimate maximum radius
    radius = np.log(sumExp)/rho

    # RETURNS
    return radius

#========================================================

def compute_arc_length(r):

    '''
    This function will compute normalized arclengths of
    the coordinates given in r.
    Remember that r is a flat vector containing the x, y, and z coordinates
    of each node in the following order:
    r = [x1, y1, z1, x2, y2, x2, x3, y3, z3, ... ]

    The first node will be at arc-length of 0.0
    The last node will be at arc-length of 1.0

    Ney Secco 2016-11
    '''

    # Get the number of nodes
    nNodes = int(len(r)/3)

    # Initialize an array that will store the arclength of each node.
    # That is, the distance, along the curve from the current node to
    # the first node of the curve.
    arcLength = np.zeros(nNodes)

    # Store coordinates of the first node (the other nodes will be covered in the loop)
    node1 = r[0:3]

    # Loop over each element to increment arcLength
    for nodeID in range(1,nNodes):

        # Get coordinates of the next node
        node2 = r[3*nodeID:3*nodeID+3]

        # Compute distance between nodes
        dist = np.linalg.norm(node1 - node2)

        # Store nodal arc-length
        arcLength[nodeID] = arcLength[nodeID-1] + dist

        # Store coordinates for the next loop
        node1 = node2

    # Normalize the arc-lengths
    arcLength = arcLength/arcLength[-1]

    # Return arc-lengths
    return arcLength

#========================================================

def redistribute_nodes_by_arc_length(r,arcLength):

    '''
    This function will receive a set of nodal coordinates defined in r and
    redistribute them along the same curve using the normalized arc-lengths
    provided in arcLength.

    Remember that r is a flat vector containing the x, y, and z coordinates
    of each node in the following order:
    r = [x1, y1, z1, x2, y2, x2, x3, y3, z3, ... ]

    The first and last nodes will remain at the same place.

    Ney Secco 2016-11
    '''

    # Import interpolation function
    from scipy.interpolate import interp1d

    # Get the number of nodes
    nNodes = int(len(r)/3)

    # Compute arc-lengths of the current curve
    origArcLength = compute_arc_length(r)

    # INTERPOLATE NEW NODES

    # Now we sample the new coordinates based on the interpolation method given by the user

    # Create interpolants for x, y, and z
    fX = interp1d(origArcLength, r[0::3])
    fY = interp1d(origArcLength, r[1::3])
    fZ = interp1d(origArcLength, r[2::3])

    # Initialize array of new coordinates
    rNew = np.zeros(r.shape)

    # Sample new points using the interpolation functions
    rNew[0::3] = fX(arcLength)
    rNew[1::3] = fY(arcLength)
    rNew[2::3] = fZ(arcLength)

    # Return the remeshed curve
    return rNew

#=============================================

def warn(message):
    print ''
    print 'WARNING!:'
    print message
    print ''

#=============================================

def error(message):
    print ''
    print 'ERROR!:'
    print message
    print ''
