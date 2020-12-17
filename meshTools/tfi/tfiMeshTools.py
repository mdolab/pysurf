'''
tfiMeshTools

Ney Secco
2016-08

This module contains tools to ease mesh generation using TFI method.

This is basically  high-level interface to use tfiLib in mesh generation problems.

The TFI geometry used in this work is as follows:

  +--> j,v
  |  +-------------------+
  V  |        top        |
 i,u |                   |
     |left          right|
     |                   |
     |       bottom      |
     +-------------------+

'''

# GENERAL IMPORTS

import numpy as np
from . import tfiLib

def computeMesh(X, Y, Z, derivDict=None):

    '''
    This is the main function that will generate the TFI mesh.
    If no derivatives are provided, the code will use linear interpolation.

    This function has no explicit outputs. It modifies X, Y, and Z instead.
    '''

    # Get problem size
    nu, nv = X.shape

    # Detect interpolation order based on the presence of derivatives.
    # Define interpolation functions based on the order.

    if derivDict is None: 

        # Let's use linear interpolation
        alpha00 = lambda u : 1-u
        alpha01 = lambda u : 0.0
        alpha10 = lambda u : u
        alpha11 = lambda u : 0.0
        beta00  = lambda v : 1-v
        beta01  = lambda v : 0.0
        beta10  = lambda v : v
        beta11  = lambda v : 0.0

    else:

        # Use cubic Hermite interpolation
        alpha00 = lambda u : 2*u**3 - 3*u**2 + 1
        alpha01 = lambda u : u**3 - 2*u**2 + u
        alpha10 = lambda u : -2*u**3 + 3*u**2
        alpha11 = lambda u : u**3 - u**2
        beta00  = lambda v : 2*v**3 - 3*v**2 + 1
        beta01  = lambda v : v**3 - 2*v**2 + v
        beta10  = lambda v : -2*v**3 + 3*v**2
        beta11  = lambda v : v**3 - v**2

    # Initialize TFI object
    tfi = tfiLib.tfi(nu, nv,
                     alpha00, alpha01, alpha10, alpha11,
                     beta00,  beta01,  beta10,  beta11)

    # Run the interpolation for each coordinate
    for coorTuple in [(X,'X'), (Y,'Y'), (Z,'Z')]:

        # Get pointer and string to current variable
        coor = coorTuple[0]
        coorName = coorTuple[1]

        # Extract derivatives
        if derivDict is None:

            dcoordu = np.zeros((2,nv))
            dcoordv = np.zeros((nu,2))
            ddcoordudv = np.zeros((2,2))

        else:

            dcoordu = derivDict['d'+coorName+'du']
            dcoordv = derivDict['d'+coorName+'dv']
            ddcoordudv = derivDict['dd'+coorName+'dudv']

        # Run interpolation
        tfi.interpolate(coor, dcoordu, dcoordv, ddcoordudv)

#======================================
#======================================
#======================================

def computeDerivatives(X, Y, Z, s0u, s0v, N0u, N0v, blendingFactor=0.6):

    '''
    This functions computes derivatives necessary for cubic Hermite TFI based on surface normals
    and user-provided grid spacings.

    INPUTS:
    X, Y, Z: array[nu,nv] --> Mesh nodal coordinates. Just the external values should be filled. The 
                              interior values in these matrices do not matter.

    s0u: array[nu,3] --> Desired marching distance for the nodes on the left (s0u[:,0]) and on the
                         right (s0u[:,1]) curves. The marching distance is the desired distance between
                         the boundary node and the first interior node.

    s0u: array[3,nv] --> Desired marching distance for the nodes on the top (s0v[0,:]) and on the
                         bottom (s0v[1,:]) curves. The marching distance is the desired distance between
                         the boundary node and the first interior node.

    N0u: array[nu,6] --> Reference surface normals (nx,ny,nz) computed at the nodes on the left (N0u[:,0:3])
                         and on the right (N0u[:,3:6]) curves.
                         N0u = [nx(0,0)  ny(0,0)  nz(0,0)  |  nx(0,-1)  ny(0,-1)  nz(0,-1)]
                               [nx(1,0)  ny(1,0)  nz(1,0)  |  nx(1,-1)  ny(1,-1)  nz(1,-1)]
                               [nx(2,0)  ny(2,0)  nz(2,0)  |  nx(2,-1)  ny(2,-1)  nz(2,-1)]
                                   :        :        :            :        :          :
                               [nx(nu,0) ny(nu,0) nz(nu,0) |  nx(nu,-1) ny(nu,-1) nz(nu,-1)]


    N0v: array[6,nv] --> Reference surface normals (nx,ny,nz) computed at the nodes on the top (N0v[0:3,:])
                         and on the bottom (N0v[3:6,:]) curves. It structure is similar to N0u, but transposed.

    blendingFactor: float (0<float<1) --> 0 -- Do not blend any information from the boundaries.
                                          1 -- Only boundary values will be used.

    '''

    #=================#
    # HELPER FUNCTION #
    #=================#

    '''
    First we define a function that computes derivatives perpendicular to an arbitrary line.
    We will call the same function to compute derivatives along the bottom, right, top, and left lines.
    '''

    def computeDerivativesOnLine(Xl, Yl, Zl,
                                 Xm0, Ym0, Zm0,
                                 Xm1, Ym1, Zm1,
                                 nm,
                                 s0l, N0l,
                                 blendingFactor,
                                 flipBCdir=False):
        
        '''
        This function computes derivatives perpendicular to the curve direction.
        In this particular case, we use l to represent the curve direction (which
        may be either u or v, depending on the case) and m to represent the direction
        perpendicular to the curve, pointing inwards.

        We also need the first nodes of the boundary curves, as shown below. They are necessary
        to blend derivatives.

        m0                   m1
        |                    |
        l0--l1--l2--l3--...--ln -> This is the current curve (Xl, Yl, Zl)


        INPUTS:

        Xl, Yl, Zl: array[nl] --> Coordinates at the nodes along the line. We assume
                                  that the nodes are linearly spaced in parametric coordinates. 

        Xm0, Ym0, Zm0: float --> X,Y, and Z coordinates of the first node of the boundary curve that starts at l=0.

        Xm1, Ym1, Zm1: float --> X,Y, and Z coordinates of the first node of the boundary curve that starts at l=-1.

        nm: integer --> number of nodes on the perpendicular curve

        s0l: array[nl] --> Desired marching distance for the first nodes along the line

        N0l: array[nl,3] --> Surface normals computed on each node of the line

        blendingFactor: float (0<float<1) --> 0 -- Do not blend any information from the boundaries.
                                              1 -- Only boundary values will be used.

        flipBCdir: logical --> should be True for bottom and right curves, so we can flip the BC derivative sign
                               when computing delta(X)/delta(m).
        '''

        # PRE-PROCESSING

        # Get number of nodes along the l curve
        nl = len(Xl)

        # Concatenate coordinates
        Fl = np.vstack([Xl, Yl, Zl])

        # COMPUTING DERIVATIVES AT THE BOUNDARIES

        '''
        In the TFI approach, we need to specify all boundaries of the domain. Consequently, some derivatives
        at the corner nodes are set implicitly as we already give the node spacing. We will compute those derivatives
        here so we can merge them with the user-specified derivatives later on.
        '''

        # Derivatives at l=0 (beginning of the line). This is just delta(X)/delta(m)
        dXdm0 = (Xm0-Xl[0])*(nm-1)
        dYdm0 = (Ym0-Yl[0])*(nm-1)
        dZdm0 = (Zm0-Zl[0])*(nm-1)

        # Derivatives at l=-1 (end of the line). This is just delta(X)/delta(m)
        dXdm1 = (Xm1-Xl[-1])*(nm-1)
        dYdm1 = (Ym1-Yl[-1])*(nm-1)
        dZdm1 = (Zm1-Zl[-1])*(nm-1)

        # Flip values according to the curve orientation
        if flipBCdir:
            dXdm0 = -dXdm0
            dYdm0 = -dYdm0
            dZdm0 = -dZdm0
            dXdm1 = -dXdm1
            dYdm1 = -dYdm1
            dZdm1 = -dZdm1

        # Now use linear interpolation to get boundary-constrained derivatives along the l curve
        dXdmBC = np.linspace(dXdm0,dXdm1,nl)
        dYdmBC = np.linspace(dYdm0,dYdm1,nl)
        dZdmBC = np.linspace(dZdm0,dZdm1,nl)

        '''
        Now we will computer another set of derivatives using user-supplied data, such
        as desired marching distance s0l
        '''

        # ESTIMATING d/dl DERIVATIVES ALONG THE LINE
        # Second order differencing

        # Initialize matrix with derivatives
        dFdl = np.zeros((3,nl))

        # Central differencing on interior nodes
        dFdl[:,1:-1] = 0.5*(Fl[:,2:] - Fl[:,:-2])*(nl-1)

        # One-sided differencing on end nodes
        dFdl[:,0] = 0.5*(-3.0*Fl[:,0] + 4.0*Fl[:,1] - Fl[:,2])*(nl-1)
        dFdl[:,-1] = 0.5*(Fl[:,-3] - 4.0*Fl[:,-2] + 3.0*Fl[:,-1])*(nl-1)

        # AREA FACTOR

        # Initialize array that will have average distances between consecutive nodes
        dNodes = np.zeros(nl)

        # Compute distances between nodes
        dd = np.linalg.norm(Fl[:,1:] - Fl[:,:-1],axis=0)

        # Fill average distances for interior nodes
        dNodes[1:-1] = 0.5*(dd[1:] + dd[:-1])

        # Compute average distances for end nodes (assuming we have ghost nodes with same distance
        # as the other neighbor)
        dNodes[0] = dd[0]
        dNodes[-1] = dd[-1]

        # Compute area factor
        dS = dNodes*s0l

        # COMPUTE PERPENDICULAR DERIVATIVES

        # Initialize matrix of perpendicular derivatives
        dFdm = np.zeros((3,nl))

        # Solve linear systems on each node to get m derivatives
        # (that is, derivatives perpendicular to the line)
        for ii in range(nl):

            # Initialize matrix
            A = np.zeros((3,3))

            # Assemble matrix
            A[0,:] = dFdl[:,ii]
            A[1,:] = np.cross(N0l[ii,:], A[0,:])
            A[2,:] = N0l[ii,:]

            # Invert and compute solution
            dFdm[:,ii] = np.linalg.inv(A).dot([0, dS[ii], 0])

        # Split derivatives
        dXdm = dFdm[0,:]
        dYdm = dFdm[1,:]
        dZdm = dFdm[2,:]

        # BLENDING DERIVATIVES

        # Now we blend the user-supplied and the BC-constrained derivatives using a beta distribution

        # Get beta-distribution coefficients based on blending factor
        if blendingFactor < 1.0:
            alpha = -1.0/(blendingFactor-1)
            beta = alpha
        else:
            # This will use only BC information
            alpha = None
            beta = None
            print('Only BC information will be used to estimate derivatives.')

        # Get weights with beta distribution
        weights = betaDist(nl,alpha,beta)

        # Now blend derivatives
        dXdm = weights*dXdm + (1-weights)*dXdmBC
        dYdm = weights*dYdm + (1-weights)*dYdmBC
        dZdm = weights*dZdm + (1-weights)*dZdmBC

        # Return derivatives
        return dXdm, dYdm, dZdm

    #===========#
    # EXECUTION #
    #===========#

    '''
    Now we will call the same derivative function for all lines
    '''

    # Get problem size
    nu, nv = X.shape

    # Check if the patch normal is consistent with the surface normal.
    # If we assume a well behaved shape, we only need to check the first cell.
    # The vector u x v of the first cell should point in the same direction
    # of the surface normal in this cell.
    uDir = [X[1,0]-X[0,0], Y[1,0]-Y[0,0], Z[1,0]-Z[0,0]]
    vDir = [X[0,1]-X[0,0], Y[0,1]-Y[0,0], Z[0,1]-Z[0,0]]
    uvDir = np.cross(uDir,vDir)
    normProj = min(uvDir.dot(N0u[0,:3]), uvDir.dot(N0u[0,:3])) # Even though the two normals should be the same, we still check both
    if normProj < 0:
        print('TFI WARNING:')
        print('The patch normal is not consistent with the surface normal.')
        print('The code will continue but the mesh may step in the wrong direction')
        print('and extrapolate the patch boundaries.')
        print('The easiest way to fix this is by flipping u and v directions by')
        print('transposing the coordinates and normal matrices.')

    # We need to multiply the marching distances to take into account
    # the normalization of u and v between 0 and 1, instead of 0 and nu or 0 and
    # nv usuallly used in hyperbolic marching
    s0u = s0u*(nu-1)*(nv-1)
    s0v = s0v*(nu-1)*(nv-1)

    # Initialize arrays with derivatives
    dXdu = np.zeros((2,nv))
    dXdv = np.zeros((nu,2))
    ddXdudv = np.zeros((2,2))
    dYdu = np.zeros((2,nv))
    dYdv = np.zeros((nu,2))
    ddYdudv = np.zeros((2,2))
    dZdu = np.zeros((2,nv))
    dZdv = np.zeros((nu,2))
    ddZdudv = np.zeros((2,2))

    # Sometimes we need to flip the node ordering so that N0 x Xu points in the Xv direction
    # and vice-versa.

    # Left line derivatives
    dXdv[:,0], dYdv[:,0], dZdv[:,0] = computeDerivativesOnLine(X[:,0], Y[:,0], Z[:,0],
                                                               X[0,1], Y[0,1], Z[0,1],
                                                               X[-1,1], Y[-1,1], Z[-1,1],
                                                               nv,
                                                               s0u[:,0], N0u[:,:3],
                                                               blendingFactor)

    # Bottom line derivatives
    dXdu[-1,::-1], dYdu[-1,::-1], dZdu[-1,::-1] = computeDerivativesOnLine(X[-1,::-1], Y[-1,::-1], Z[-1,::-1],
                                                                           X[-2,-1], Y[-2,-1], Z[-2,-1],
                                                                           X[-2,0], Y[-2,0], Z[-2,0],
                                                                           nu,
                                                                           s0v[-1,::-1], N0v[3:,::-1].T,
                                                                           blendingFactor,
                                                                           flipBCdir=True)

    # Right line derivatives
    dXdv[:,-1], dYdv[:,-1], dZdv[:,-1] = computeDerivativesOnLine(X[:,-1], Y[:,-1], Z[:,-1],
                                                                  X[0,-2], Y[0,-2], Z[0,-2],
                                                                  X[-1,-2], Y[-1,-2], Z[-1,-2],
                                                                  nv,
                                                                  s0u[:,-1], N0u[:,3:],
                                                                  blendingFactor,
                                                                  flipBCdir=True)

    # Top line derivatives
    dXdu[0,::-1], dYdu[0,::-1], dZdu[0,::-1] = computeDerivativesOnLine(X[0,::-1], Y[0,::-1], Z[0,::-1],
                                                                        X[1,-1], Y[1,-1], Z[1,-1],
                                                                        X[1,0], Y[1,0], Z[1,0],
                                                                        nu,
                                                                        s0v[0,::-1], N0v[:3,::-1].T,
                                                                        blendingFactor)
    
    # Now we estimate the cross derivatives at the corners using finite diferencing.
    # We average the estimates in each orthogonal direction

    # Get indices for "middle" nodes
    um = int((nu-1)/2)
    vm = int((nv-1)/2)

    # NODE (0,0)
    '''
    ddXdudv_u = (dXdv[um,0] - dXdv[0,0])*(nu-1)/um
    ddXdudv_v = (dXdu[0,vm] - dXdu[0,0])*(nv-1)/vm
    ddXdudv[0,0] = 0.5*(ddXdudv_u + ddXdudv_v)

    ddYdudv_u = (dYdv[um,0] - dYdv[0,0])*(nu-1)/um
    ddYdudv_v = (dYdu[0,vm] - dYdu[0,0])*(nv-1)/vm
    ddYdudv[0,0] = 0.5*(ddYdudv_u + ddYdudv_v)

    ddZdudv_u = (dZdv[um,0] - dZdv[0,0])*(nu-1)/um
    ddZdudv_v = (dZdu[0,vm] - dZdu[0,0])*(nv-1)/vm
    ddZdudv[0,0] = 0.5*(ddZdudv_u + ddZdudv_v)

    # NODE (0,1)

    ddXdudv_u = (dXdv[um,-1] - dXdv[0,-1])*(nu-1)/um
    ddXdudv_v = (dXdu[0,-1] - dXdu[0,vm])*(nv-1)/vm
    ddXdudv[0,1] = 0.5*(ddXdudv_u + ddXdudv_v)

    ddYdudv_u = (dYdv[um,-1] - dYdv[0,-1])*(nu-1)/um
    ddYdudv_v = (dYdu[0,-1] - dYdu[0,vm])*(nv-1)/vm
    ddYdudv[0,1] = 0.5*(ddYdudv_u + ddYdudv_v)

    ddZdudv_u = (dZdv[um,-1] - dZdv[0,-1])*(nu-1)/um
    ddZdudv_v = (dZdu[0,-1] - dZdu[0,vm])*(nv-1)/vm
    ddZdudv[0,1] = 0.5*(ddZdudv_u + ddZdudv_v)

    # NODE (1,0)

    ddXdudv_u = (dXdv[-1,0] - dXdv[um,0])*(nu-1)/um
    ddXdudv_v = (dXdu[-1,vm] - dXdu[-1,0])*(nv-1)/vm
    ddXdudv[1,0] = 0.5*(ddXdudv_u + ddXdudv_v)

    ddYdudv_u = (dYdv[-1,0] - dYdv[um,0])*(nu-1)/um
    ddYdudv_v = (dYdu[-1,vm] - dYdu[-1,0])*(nv-1)/vm
    ddYdudv[1,0] = 0.5*(ddYdudv_u + ddYdudv_v)

    ddZdudv_u = (dZdv[-1,0] - dZdv[um,0])*(nu-1)/um
    ddZdudv_v = (dZdu[-1,vm] - dZdu[-1,0])*(nv-1)/vm
    ddZdudv[1,0] = 0.5*(ddZdudv_u + ddZdudv_v)

    # NODE (1,1)

    ddXdudv_u = (dXdv[-1,-1] - dXdv[um,-1])*(nu-1)/um
    ddXdudv_v = (dXdu[-1,-1] - dXdu[-1,vm])*(nv-1)/vm
    ddXdudv[1,1] = 0.5*(ddXdudv_u + ddXdudv_v)

    ddYdudv_u = (dYdv[-1,-1] - dYdv[um,-1])*(nu-1)/um
    ddYdudv_v = (dYdu[-1,-1] - dYdu[-1,vm])*(nv-1)/vm
    ddYdudv[1,1] = 0.5*(ddYdudv_u + ddYdudv_v)

    ddZdudv_u = (dZdv[-1,-1] - dZdv[um,-1])*(nu-1)/um
    ddZdudv_v = (dZdu[-1,-1] - dZdu[-1,vm])*(nv-1)/vm
    ddZdudv[1,1] = 0.5*(ddZdudv_u + ddZdudv_v)
    '''
    # Assemble derivatives in a dictionary
    derivDict = {'dXdu':dXdu,
                 'dXdv':dXdv,
                 'ddXdudv':ddXdudv,
                 'dYdu':dYdu,
                 'dYdv':dYdv,
                 'ddYdudv':ddYdudv,
                 'dZdu':dZdu,
                 'dZdv':dZdv,
                 'ddZdudv':ddZdudv}

    # Return derivatives
    return derivDict

#======================================
#======================================
#======================================

def computeArcLengths(X,Y,Z):

    '''
    This function will compute the normalized arc-lengths
    '''

#======================================
#======================================
#======================================

def exportPlot3d(X, Y, Z, filename):

    '''
    This function exports a 3D mesh in plot3d format.
    '''

    # IMPORTS
    from plot3d_interface import Grid, export_plot3d

    # Initialize grid object
    myGrid = Grid()

    # Expand the coordinate matrices
    X3d = np.array([X])
    Y3d = np.array([Y])
    Z3d = np.array([Z])

    # Add mesh to the grid
    myGrid.add_block(X3d, Y3d, Z3d)

    # Export grid
    export_plot3d(myGrid, filename)

#======================================
#======================================
#======================================

def betaDist(n,alpha,beta):

    '''
    This functions uses the beta distribution to compute n weighing
    factors ranging from 0 to 1.

    alpha and beta should be >= 1

    If alpha=beta=1.0, we return all weights as ones
    If alpha=beta=None, we return all weights as zeros
    '''

    # Check if we have Nones
    if alpha is None and beta is None:
        return np.zeros(n)

    # Show warning for bad inputs
    if alpha < 1.0 or beta < 1.0:
        print('Warning: the beta distribution function in this code only')
        print('works properly for alpha >= 1.0 and beta >= 1.0.')

    # If alpha=beta=1, all weights should be 1.0
    if alpha <= 1.0 and beta <= 1.0:
        return np.ones(n)

    # Compute normalization factors so that the maximum output is 1.0
    xopt = (1-alpha)/(2-alpha-beta)
    K = 1/(xopt**(alpha-1)*(1-xopt)**(beta-1))

    # Now compute weights
    x = np.linspace(0,1,n)
    w = K*x**(alpha-1)*(1-x)**(beta-1)

    # Return weights
    return w
