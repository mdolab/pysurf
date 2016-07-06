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
- check 3D metric correction
- check smoothing
- update plot3d export function
- add compute normal and projection steps

-set up scaling factor based on distance
-blend the angle-based dissipation coefficient
'''

def create_mesh(rBaseline, NBaseline, bc1, bc2, sBaseline, numLayers, extension, \
                epsE0 = 1.0, theta = 0.0, \
                alphaP0 = 0.25, num_smoothing_passes = 0, \
                nuArea = 0.16, num_area_passes = 0, \
                sigmaSplay = 0.3, \
                metricCorrection = False, \
                ratioGuess = 20):

    # This function will march N steps

    # IMPORTS
    from numpy import zeros, copy

    # Initialize 2D array that contains all the grid points
    R = zeros((numLayers,len(rBaseline)))
    R[0,:] = rBaseline
    rm1 = rBaseline

    # Initialize step size and total marched distance
    s = sBaseline
    s_tot = 0

    # Find the characteristic radius of the mesh
    radius = findRadius(rBaseline)

    # Find the desired marching distance
    s_max = radius*(extension-1)

    # Compute the growth ratio necessary to match this distance
    growth = findRatio(s_max, sBaseline, numLayers, ratioGuess)

    # Print growth ratio
    print 'Growth ratio: ',growth

    '''
    ==========================================
    Now let's define the functions that do the heavy work.
    These functions are the actual hyperbolic generation
    calculations.
    I define them here so that the several optional parameters
    are readily available to them.
    ==========================================
    '''

    def sub_iteration(r0, rm1, Aprev, s, s_tot, layer_index):

        # IMPORTS
        from numpy.linalg import solve

        # Smooth coordinates
        #r0 = smoothing(r0,layer_index+2)

        # Compute area factor
        A = area_factor(r0,s)

        # Repeat areas for the first iteration
        # As we don't have a previous layer, we just repeat areas
        if isinstance(Aprev,int):
            Aprev = A[:]

        '''
        Add function to compute normals HERE
        '''
        # Temporary normals
        N0 = NBaseline

        # Generate matrices of the linear system
        K,f = compute_matrices(r0, rm1, N0, Aprev, A, layer_index, s_tot)

        # Store areas to be used by the next layer
        Aprev = A[:]

        # Solve the linear system
        dr = solve(K,f)

        # Update r
        r_next = r0 + dr

        # RETURNS
        return r_next, A

    def area_factor(r0,s):

        # This function computes the area factor distribution

        # IMPORT
        from numpy import hstack, reshape
        from numpy.linalg import norm

        # Extrapolate the end points
        r0minus1 = 2*r0[0:3] - r0[3:6]
        r0plus1 = 2*r0[-3:] - r0[-6:-3]
        r0_extrap = hstack([r0minus1, r0, r0plus1])

        # Reshape so we have 3D array
        R0_extrap = reshape(r0_extrap,(-1,3)).T

        # Compute the distance of each node to its neighbors
        d = 0.5*(norm(R0_extrap[:,1:-1] - R0_extrap[:,:-2],axis=0) + norm(R0_extrap[:,2:] - R0_extrap[:,1:-1],axis=0))

        # Multiply distances by the step size to get the areas
        A = d*s

        # Do the requested number of averagings
        for index in xrange(num_area_passes):
            # Store previous values
            Aplus = A[1]
            Aminus = A[-2]
            # Do the averaging
            A[1:-1] = (1-nuArea)*A[1:-1] + nuArea/2*(A[:-2] + A[2:])
            A[0] = (1-nuArea)*A[0] + nuArea*Aplus
            A[-1] = (1-nuArea)*A[-1] + nuArea*Aminus

        # RETURNS
        return A

    def metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index):

        # This function applies the metric correction algorithm to update x0_eta and y0_eta

        '''

        # Compute distance vectors (Eq. 6.9)
        r_plus = r0_next - r0_curr
        r_minus = r0_prev - r0_curr

        # Compute corrected zeta derivatives (Eq. 7.2)
        r_zeta = 0.25*(abs(r_plus) + abs(r_minus))*(r_plus/abs(r_plus) - r_minus/abs(r_minus))
        xp_zeta = r_zeta[0]
        yp_zeta = r_zeta[1]

        # Find the adjusted eta derivatives (Eq. 7.1)
        gamma_p = xp_zeta*xp_zeta + yp_zeta*yp_zeta
        xp_eta = -dA_curr/gamma_p*yp_zeta
        yp_eta = dA_curr/gamma_p*xp_zeta

        # Compute metric correction factor
        nu_mc = 2**(-layer_index)

        # Apply metric correction (Eq. 7.3)
        x0_eta = (1-nu_mc)*x0_eta + nu_mc*xp_eta
        y0_eta = (1-nu_mc)*y0_eta + nu_mc*yp_eta

        '''

        # RETURNS
        return x0_eta, y0_eta

    def compute_matrices(r0,rm1,N0,Aprev,A,layer_index,s_tot):

        # This function computes the derivatives r_zeta and r_eta for a group
        # of coordinates given in r.
        # It shoud be noted that r is a 1d array that contains x and y
        # coordinates of all points in a layer:
        # r0 = [x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]
        #
        # N0 is a 2D array with the surface normals at each r0 point
        # N0 = nx1 nx2 nx3 nx4 ...
        #      ny1 ny2 ny3 ny4 ...
        #      nz1 nz2 nz3 nz4

        # IMPORTS
        from numpy import zeros, sqrt, eye, array, arctan2, pi, cross
        from numpy.linalg import norm, inv

        # Find number of nodes
        numNodes = int(len(r0)/3)

        # Initialize arrays
        K = zeros((3*numNodes,3*numNodes))
        f = zeros(3*numNodes)

        # Define assembly routines
        def central_matrices(prev_index,curr_index,next_index):

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
            r0_eta = B0inv.dot(array([0, Aprev[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            '''
            if metricCorrection:
                r0_prev = r0[2*(prev_index):2*(prev_index)+2]
                r0_curr = r0[2*(curr_index):2*(curr_index)+2]
                r0_next = r0[2*(next_index):2*(next_index)+2]
                dA_curr = dA[curr_index]
                x0_eta, y0_eta = metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index)
            '''

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(next_index):3*(next_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(prev_index):3*(prev_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(next_index):3*(next_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(prev_index):3*(prev_index)+3]-r0[3*(index):3*(index)+3])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = give_angle(r0[3*(prev_index):3*(prev_index)+3],
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
                epsE, epsI = dissipation_coefficients(layer_index, r0_xi, r0_eta, d_sensor, angle)

                # Compute RHS components
                B0invg = B0inv.dot(array([0, A[curr_index], 0]))
                De = epsE*(r0[3*(prev_index):3*(prev_index)+3] - 2*r0[3*(index):3*(index)+3] + r0[3*(next_index):3*(next_index)+3])

                # Compute block matrices
                L_block = -0.5*(1+theta)*C0 - epsI*eye(3)
                M_block = (1 + 2*epsI)*eye(3)
                N_block = 0.5*(1+theta)*C0 - epsI*eye(3)
                f_block = B0invg - De

                # Populate matrix
                K[3*(curr_index):3*(curr_index)+3,3*(prev_index):3*(prev_index)+3] = L_block
                K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = M_block
                K[3*(curr_index):3*(curr_index)+3,3*(next_index):3*(next_index)+3] = N_block
                f[3*(curr_index):3*(curr_index)+3] = f_block

        def forward_matrices(curr_index,next_index,nextnext_index):

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
            r0_eta = B0inv.dot(array([0, Aprev[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            '''
            if metricCorrection:
                r0_prev = r0[2*(prev_index):2*(prev_index)+2]
                r0_curr = r0[2*(curr_index):2*(curr_index)+2]
                r0_next = r0[2*(next_index):2*(next_index)+2]
                dA_curr = dA[curr_index]
                x0_eta, y0_eta = metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index)
            '''

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(next_index):3*(next_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(nextnext_index):3*(nextnext_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(next_index):3*(next_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(nextnext_index):3*(nextnext_index)+3]-r0[3*(index):3*(index)+3])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

            # Compute C0 = B0inv*A0
            C0 = B0inv.dot(A0)

            # Compute smoothing coefficients
            epsE, epsI = dissipation_coefficients(layer_index, r0_xi, r0_eta, d_sensor, angle)

            # Compute RHS components
            B0invg = B0inv.dot(array([0, A[curr_index], 0]))
            De = epsE*(r0[3*(nextnext_index):3*(nextnext_index)+3] - 2*r0[3*(next_index):3*(next_index)+3] + r0[3*(curr_index):3*(curr_index)+3])

            # Compute block matrices
            L_block = -0.5*(1+theta)*C0 - epsI*eye(3)
            M_block = 2*(1+theta)*C0 + 2*epsI*eye(3)
            N_block = -1.5*(1+theta)*C0 + (1-epsI)*eye(3)
            f_block = B0invg - De

            # Populate matrix
            K[3*(curr_index):3*(curr_index)+3,3*(nextnext_index):3*(nextnext_index)+3] = L_block
            K[3*(curr_index):3*(curr_index)+3,3*(next_index):3*(next_index)+3] = M_block
            K[3*(curr_index):3*(curr_index)+3,3*(curr_index):3*(curr_index)+3] = N_block
            f[3*(curr_index):3*(curr_index)+3] = f_block

        def backward_matrices(curr_index,prev_index,prevprev_index):

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
            r0_eta = B0inv.dot(array([0, Aprev[curr_index], 0]))
            x0_eta = r0_eta[0]
            y0_eta = r0_eta[1]
            z0_eta = r0_eta[2]

            # Assemble A0 matrix
            A0 = array([[x0_eta, y0_eta, z0_eta],
                        [ny*z0_eta-nz*y0_eta, nz*x0_eta-nx*z0_eta, nx*y0_eta-ny*x0_eta],
                        [0, 0, 0]])

            '''
            if metricCorrection:
                r0_prev = r0[2*(prev_index):2*(prev_index)+2]
                r0_curr = r0[2*(curr_index):2*(curr_index)+2]
                r0_next = r0[2*(next_index):2*(next_index)+2]
                dA_curr = dA[curr_index]
                x0_eta, y0_eta = metric_correction(x0_eta, y0_eta, dA_curr, r0_prev, r0_curr, r0_next, layer_index)
            '''

            # Compute grid distribution sensor (Eq. 6.8a)
            dnum = norm(rm1[3*(prev_index):3*(prev_index)+3]-rm1[3*(index):3*(index)+3]) + norm(rm1[3*(prevprev_index):3*(prevprev_index)+3]-rm1[3*(index):3*(index)+3])
            dden = norm(r0[3*(prev_index):3*(prev_index)+3]-r0[3*(index):3*(index)+3]) + norm(r0[3*(prevprev_index):3*(prevprev_index)+3]-r0[3*(index):3*(index)+3])
            d_sensor = dnum/dden

            # Compute the local grid angle based on the neighbors
            angle = pi

             # Compute C0 = B0inv*A0
            C0 = B0inv.dot(A0)

            # Compute smoothing coefficients
            epsE, epsI = dissipation_coefficients(layer_index, r0_xi, r0_eta, d_sensor, angle)

            # Compute RHS components
            B0invg = B0inv.dot(array([0, A[curr_index], 0]))
            De = epsE*(r0[3*(prevprev_index):3*(prevprev_index)+3] - 2*r0[3*(prev_index):3*(prev_index)+3] + r0[3*(curr_index):3*(curr_index)+3])

            # Compute block matrices
            L_block = 0.5*(1+theta)*C0 - epsI*eye(3)
            M_block = -2*(1+theta)*C0 + 2*epsI*eye(3)
            N_block = 1.5*(1+theta)*C0 + (1-epsI)*eye(3)
            f_block = B0invg - De

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
                forward_matrices(curr_index,next_index,nextnext_index)

            elif bc1 is 'continuous':

                # Name indexes
                prev_index = numNodes-2
                curr_index = index
                next_index = index+1

                # Call assembly routine
                central_matrices(prev_index,curr_index,next_index)

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
                f[3*index:3*index+3] = array([A[index]*(1-sigmaSplay), 0, 0])
                #f[3*index:3*index+3] = array([0, s*(1-sigmaSplay)*norm(d_vec_rot)])

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
            central_matrices(prev_index,curr_index,next_index)

        for index in [numNodes-1]:

            if bc2 is 'free':

                # Name indexes
                curr_index = index
                prev_index = index-1
                prevprev_index = index-2

                # Call assembly routine
                backward_matrices(curr_index,prev_index,prevprev_index)

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
                f[3*index:3*index+3] = array([A[index]*(1-sigmaSplay), 0, 0])

            elif bc1 is 'constX':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[0, 0, 0],[0, -1, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constY':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, 0, 0],[0, 0, -1]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

            elif bc1 is 'constZ':

                # Populate matrix
                K[3*index:3*index+3,3*(index-1):3*(index-1)+3] = [[-1, 0, 0],[0, -1, 0],[0, 0, 0]]
                K[3*index:3*index+3,3*index:3*index+3] = eye(3)
                f[3*index:3*index+3] =  [0, 0, 0]

        # RETURNS
        return K,f

    def smoothing(r,eta):

        '''

        # This function does the grid smoothing

        # IMPORTS
        from numpy import copy, zeros
        from numpy.linalg import norm

        # Find number of nodes
        numNodes = int(len(r)/2)

        # Loop over the desired number of smoothing passes
        for index_pass in range(num_smoothing_passes):

            # Initialize array of smoothed coordinates
            r_smooth = zeros(2*numNodes)

            # Copy the edge nodes
            r_smooth[:2] = r[:2]
            r_smooth[-2:] = r[-2:]

            # Smooth every node
            for index in xrange(1,numNodes-1):

                # Get coordinates
                r_curr = r[2*(index):2*(index)+2]
                r_next = r[2*(index+1):2*(index+1)+2]
                r_prev = r[2*(index-1):2*(index-1)+2]

                # Compute distances
                lp = norm(r_next - r_curr)
                lm = norm(r_curr - r_prev)

                # Compute alpha'
                alphaP = min(alphaP0, alphaP0*(eta-2)/numLayers)

                # Compute smoothed coordinates
                r_smooth[2*index:2*index+2] = (1-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)

            # Copy coordinates to allow next pass
            r = copy(r_smooth)

        '''

        # RETURNS
        return r

    def dissipation_coefficients(layer_index, r0_xi, r0_eta, d_sensor, angle):

        # IMPORTS
        from numpy import sqrt, cos, pi
        from numpy.linalg import norm

        # Compute N (Eq. 6.3)
        N = norm(r0_eta)/norm(r0_xi)

        # Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layer_index+2
        ltrans = int(3/4*numLayers)
        if l <= ltrans:
            Sl = sqrt((l-1)/(ltrans-1))
        else:
            Sl = sqrt((ltrans-1)/(numLayers-1))

        # Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max([d_sensor**(2/Sl), 0.1])

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

    # MARCH!!!
    for layer_index in range(numLayers-1):
        # Get the previous coordinates
        r0 = R[layer_index,:]

        # Compute the total marched distance
        s_tot = s_tot + s

        if layer_index == 0:
            Aprev = 0

        # Subiter
        for index_sub in range(1):
            r_next, Aprev = sub_iteration(r0, rm1, Aprev, s, s_tot, layer_index)

        # Store grid points
        R[layer_index+1,:] = r_next

        # Update step size
        s = s*growth

        # Update rm1
        rm1 = R[layer_index,:]

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

def give_angle(r0,r1,r2,N1):

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

def findRatio(s_max, s0, numLayers, ratioGuess):

    '''
    This function returns the geometrical progression ratio that satisfies
    the farfield distance and the number of cells. Newton search is used
    INPUTS
    s_max: distance that should be reached
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
       R = s0*(1-q**(numLayers-1)) - s_max*(1-q)

       # Residual derivative
       Rdot = -(numLayers-1)*s0*q**(numLayers-2) + s_max

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

def plot_grid(X,Y,Z,show=False):

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

def export_plot3d(X,Y,filename,zSpan=1,zNodes=2):

    '''
    This function exports a 3D mesh in plot3d format.
    The user specifies the span and number of nodes in the z direction.
    '''

    # IMPORTS
    from plot3d_writer import Grid, export_plot3d
    from numpy import array, copy, ones, linspace

    # Initialize grid object
    myGrid = Grid()

    # Expand the coordinate matrices
    X3d = array([copy(X) for dummy in range(zNodes)])
    Y3d = array([copy(Y) for dummy in range(zNodes)])

    # Create the span-wise coordinates
    Z3d = array([z*ones(X.shape) for z in linspace(0,zSpan,zNodes)])

    # Add block to the grid
    myGrid.add_block(X3d, Y3d, Z3d)

    # Export grid
    export_plot3d(myGrid, filename)
