"""
tfiLib.py

Ney Secco
2016-08

This file defines the tfi library

The TFI geometry used in this work is as follows:

  +--> j,v
  |  +-------------------+
  V  |        top        |
 i,u |                   |
     |left          right|
     |                   |
     |       bottom      |
     +-------------------+


"""

# IMPORTS

import numpy as np


class tfi(object):

    """
    This is the main tfi object

    The TFI creates an interpolant g for function f based on the
    following expression

    g(u,v) = f(u0,v)*alpha0(u) + f(u1,v)*alpha1(u) +
           + f(u,v0)*beta0(v)  + f(u,v1)*beta1(v)  -
           - f(u0,v0)*alpha0(u)*beta0(v) - f(u1,v0)*alpha1(u)*beta0(v) -
           - f(u0,v1)*alpha0(u)*beta1(v) - f(u1,v1)*alpha1(u)*beta1(v)
    """

    def __init__(
        self, nu, nv, alpha00, alpha01, alpha10, alpha11, beta00, beta01, beta10, beta11, u="uniform", v="uniform"
    ):

        """
        INPUTS
        nu:     integer --> number of nodes in u direction (including end points)
        nv:     integer --> number of nodes in v direction (including end points)
        alpha00: function alpha00(u) --> univariate function used to interpolate functions in u
                                         alpha00(0) = 1.0 and alpha00(1) = 0.0
        alpha01: function alpha01(u) --> univariate function used to interpolate derivatives in u
        alpha10: function alpha10(u) --> univariate function used to interpolate functions in u
                                         alpha10(0) = 0.0 and alpha10(1) = 1.0
        alpha11: function alpha11(u) --> univariate function used to interpolate derivatives in u
        beta00:  function beta00(v) --> univariate function used to interpolate functions in v
                                        beta0(0) = 1.0 and beta0(1) = 0.0
        beta01:  function beta01(v) --> univariate function used to interpolate derivatives in v
        beta10:  function beta10(v) --> univariate function used to interpolate functions in v
                                        beta1(0) = 0.0 and beta1(1) = 1.0
        beta11:  function beta11(v) --> univariate function used to interpolate derivatives in v

        CHANGES
        what self.setTFI changes
        """

        # Store the function handles
        # but do not compute matrices yet because we did not set the size
        self.setTFI(nu, nv, alpha00, alpha01, alpha10, alpha11, beta00, beta01, beta10, beta11, u, v)

    # ===============================================

    def setTFI(
        self,
        nu=None,
        nv=None,
        alpha00=None,
        alpha01=None,
        alpha10=None,
        alpha11=None,
        beta00=None,
        beta01=None,
        beta10=None,
        beta11=None,
        u=None,
        v=None,
    ):

        """
        This function will update the TFI size and interpolation functions.
        The interpolation matrices will also be computed.

        Only the "non-None" parameters will be updated.

        INPUTS
        nu:     integer --> number of nodes in u direction (including end points)

        nv:     integer --> number of nodes in v direction (including end points)

        alpha00: function alpha00(u) --> univariate function used to interpolate functions in u
                                         alpha00(0) = 1.0 and alpha00(1) = 0.0

        alpha01: function alpha01(u) --> univariate function used to interpolate derivatives in u

        alpha10: function alpha10(u) --> univariate function used to interpolate functions in u
                                         alpha10(0) = 0.0 and alpha10(1) = 1.0

        alpha11: function alpha11(u) --> univariate function used to interpolate derivatives in u

        beta00:  function beta00(v) --> univariate function used to interpolate functions in v
                                        beta0(0) = 1.0 and beta0(1) = 0.0

        beta01:  function beta01(v) --> univariate function used to interpolate derivatives in v

        beta10:  function beta10(v) --> univariate function used to interpolate functions in v
                                        beta1(0) = 0.0 and beta1(1) = 1.0

        beta11:  function beta11(v) --> univariate function used to interpolate derivatives in v

        u: array[nu] --> Parametric coordinates (ranging from 0 to 1) that correspond to each u-layer.

        v: array[nv] --> Parametric coordinates (ranging from 0 to 1) that correspond to each v-layer.

        CHANGES
        self.__nu
        self.__nv
        self.__u
        self.__v
        self.__baseFuncs.alpha00
        self.__baseFuncs.alpha01
        self.__baseFuncs.alpha10
        self.__baseFuncs.alpha11
        self.__baseFuncs.beta00
        self.__baseFuncs.beta01
        self.__baseFuncs.beta10
        self.__baseFuncs.beta11
        what self.__computeInterpolationArrays changes
        """

        # Values that have None will not be updated

        # Update current discretization
        if nu is not None:
            self.__nu = nu
        if nv is not None:
            self.__nv = nv

        # Update function handles
        if alpha00 is not None:
            self.__alpha00 = alpha00
        if alpha01 is not None:
            self.__alpha01 = alpha01
        if alpha10 is not None:
            self.__alpha10 = alpha10
        if alpha11 is not None:
            self.__alpha11 = alpha11

        if beta00 is not None:
            self.__beta00 = beta00
        if beta01 is not None:
            self.__beta01 = beta01
        if beta10 is not None:
            self.__beta10 = beta10
        if beta11 is not None:
            self.__beta11 = beta11

        # Update parametric coordinates
        if u is not None:
            self.__u = u
        if v is not None:
            self.__v = v

        # Update interpolation arrays
        self.__computeInterpolationArrays()

    # ===============================================

    def interpolate(self, F, dFdu, dFdv, ddFdudv):

        """
        This function will perform the interpolation

        INPUTS
        F: array[nu,nv] --> This array should have the boundary values for
                            the interpolation. The values in the interior
                            of the array will not be used. They will be replaced
                            by the interpolated values.
               F = [ f(u_0,v_0)  f(u_0,v_1)  f(u_0,v_2)  ... f(u_0,v_nv)  ]
                   [ f(u_1,v_0)     0.0         0.0      ... f(u_1,v_nv)  ]
                   [ f(u_2,v_0)     0.0         0.0      ... f(u_2,v_nv)  ]
                   [     ...        ...         ...      ...      ...     ]
                   [ f(u_nu,v_0) f(u_nu,v_1) f(u_nu,v_2) ... f(u_nu,v_nv) ]

        INPUTS
        F: array[nu,nv] --> array containing the values f(u,v). Only the boundary
                            values will be used in this function. We will assume
                            that the values in f are in homogeneous parametric spacing

        dFdu: array[2,nv] --> Array with dFdu for nodes on the top [0,:]
                              and bottom nodes [1,:]

        dFdv: array[nu,2] --> Array with dFdv for nodes on the left [:,0]
                              and right nodes [:,1]

        ddFdudv: array[2,2] --> Array with ddFdudv for the corner nodes
                                ddFdudv[0,0]: top left corner
                                ddFdudv[1,0]: bottom left corner
                                ddFdudv[0,1]: top right corner
                                ddFdudv[1,1]: bottom right corner

        RETURNS
        This function has no explicit outputs. It modifies one of the inputs instead.

        MODS
        F: array[nu,nv] --> Updated array with interpolated values in the interior
                            elements.
               F = [ f(u_0,v_0)  f(u_0,v_1)  ... f(u_0,v_nv)  ]
                   [ f(u_1,v_0)  f(u_1,v_1)  ... f(u_1,v_nv)  ]
                   [     ...       ...             ...        ]
                   [ f(u_nu,v_0) f(u_nu,v_1) ... f(u_nu,v_nv) ]

        ATTENTION: F cannot be an array of integers because when we assign float to
                   an integer array, the floats becomes an integer as well!
        """

        # First check if the given array and the current TFI object have compatible sizes
        if F.shape != (self.__nu, self.__nv):

            # Print warnings
            print("WARNING:")
            print("The TFI object is set to a size different from the given matrix.")
            print("Reinitializing TFI with the correct size.")
            print("You can also do this using the tfi.setTFI function.")

            # Update tfi
            self.setTFI(nu=F.shape[0], nv=F.shape[1])

        # Now check if F is an array of integers. If this is the case, we must
        # stop the program because Python will transform floats back to integers...
        if type(F[0, 0]) == np.int64:

            # Print error message
            print("ERROR:")
            print("F should be an array of floats. Check if you are providing")
            print("integer values to be interpolated and convert them to floats.")

            # Finish program
            quit()

        # Get reference values arrays
        Fu, Fv, Fuv = self.__computeReferenceArrays(F, dFdu, dFdv, ddFdudv)

        # Do the interpolation
        F[:, :] = 1.0 * self.__A.T.dot(Fu) + 1.0 * Fv.T.dot(self.__B) - 1.0 * self.__A.T.dot(Fuv.dot(self.__B))

    # ===============================================

    def __computeInterpolationArrays(self):

        """
        This function will compute arrays of size (4 x nu) and (4 x nv) that will store
        the interpolation function values throughout the domain.
        This way, we can use the the same arrays to interpolate multiple functions.

        We can also provide non-uniform values for u and v. The literature recommends
        setting parametric coordinates based on the actual arc-length of the curve.

        CHANGES
        self.__A
        self.__B
        """

        # Get TFI size
        nu = self.__nu
        nv = self.__nv

        # Get parametric coordinates
        u = self.__u
        v = self.__v

        # Initialize parametric coordinates (if None is provided) or check
        # their size
        if u == "uniform":
            u = np.linspace(0, 1, nu)
        else:
            if len(u) is not nu:
                print("")
                print("TFI error: __computeInterpolation Arrays")
                print(" The number of parametric coordinates provided does not match")
                print(" the number of points.")
                print("")

        if v == "uniform":
            v = np.linspace(0, 1, nv)
        else:
            if len(v) is not nv:
                print("")
                print("TFI error: __computeInterpolation Arrays")
                print(" The number of parametric coordinates provided does not match")
                print(" the number of points.")
                print("")

        # Initialize matrices
        self.__A = np.zeros((4, nu))
        self.__B = np.zeros((4, nv))

        # Fill matrices
        for ii in range(nu):
            self.__A[0, ii] = self.__alpha00(u[ii])
            self.__A[1, ii] = self.__alpha01(u[ii])
            self.__A[2, ii] = self.__alpha10(u[ii])
            self.__A[3, ii] = self.__alpha11(u[ii])

        for jj in range(nv):
            self.__B[0, jj] = self.__beta00(v[jj])
            self.__B[1, jj] = self.__beta01(v[jj])
            self.__B[2, jj] = self.__beta10(v[jj])
            self.__B[3, jj] = self.__beta11(v[jj])

    # ===============================================

    def __computeReferenceArrays(self, F, dFdu, dFdv, ddFdudv):

        """
        This function will create arrays with the function reference values
        so that the entire interpolation computation is vectorized

        INPUTS
        F: array[nu,nv] --> array containing the values f(u,v). Only the boundary
                            values will be used in this function. We will assume
                            that the values in f are in homogeneous parametric spacing

        dFdu: array[2,nv] --> Array with dFdu for nodes on the top [0,:]
                              and bottom nodes [1,:]

        dFdv: array[nu,2] --> Array with dFdv for nodes on the left [:,0]
                              and right nodes [:,1]

        ddFdudv: array[2,2] --> Array with ddFdudv for the corner nodes
                                ddFdudv[0,0]: top left corner
                                ddFdudv[1,0]: bottom left corner
                                ddFdudv[0,1]: top right corner
                                ddFdudv[1,1]: bottom right corner

        RETURNS
        Fu0v: array[nu,nv] --> array with reference values
        Fu1v: array[nu,nv] --> array with reference values
        Fuv0: array[nu,nv] --> array with reference values
        Fuv1: array[nu,nv] --> array with reference values
        Fu0v0: float --> float with corner values
        Fu0v1: float --> float with corner values
        Fu1v0: float --> float with corner values
        Fu1v1: float --> float with corner values
        """

        # Get TFI size
        nu = self.__nu
        nv = self.__nv

        # Initialize matrices
        Fu = np.zeros((4, nv))  # Values of constant u: [f(u0,v), f'(u0,v), f(u1,v), f'(u1,v)]
        Fv = np.zeros((4, nu))  # Values of constant v: [f(u,v0), f'(u,v0), f(u,v1), f'(u,v1)]
        Fuv = np.zeros((4, 4))

        # Fill matrices
        for ii in range(nu):
            Fv[0, ii] = F[ii, 0]
            Fv[1, ii] = dFdv[ii, 0]
            Fv[2, ii] = F[ii, -1]
            Fv[3, ii] = dFdv[ii, -1]

        for jj in range(nv):
            Fu[0, jj] = F[0, jj]
            Fu[1, jj] = dFdu[0, jj]
            Fu[2, jj] = F[-1, jj]
            Fu[3, jj] = dFdu[-1, jj]

        # Get corner values to assemble Fuv
        Fuv[0, 0] = F[0, 0]
        Fuv[1, 0] = dFdu[0, 0]
        Fuv[2, 0] = F[-1, 0]
        Fuv[3, 0] = dFdu[-1, 0]

        Fuv[0, 1] = dFdv[0, 0]
        Fuv[1, 1] = ddFdudv[0, 0]
        Fuv[2, 1] = dFdv[-1, 0]
        Fuv[3, 1] = ddFdudv[-1, 0]

        Fuv[0, 2] = F[0, -1]
        Fuv[1, 2] = dFdu[0, -1]
        Fuv[2, 2] = F[-1, -1]
        Fuv[3, 2] = dFdu[-1, -1]

        Fuv[0, 3] = dFdv[0, -1]
        Fuv[1, 3] = ddFdudv[0, -1]
        Fuv[2, 3] = dFdv[-1, -1]
        Fuv[3, 3] = ddFdudv[-1, -1]

        # RETURNS
        return Fu, Fv, Fuv


# ===============================================
# ===============================================
# ===============================================

# AUXILIARY FUNCTIONS
