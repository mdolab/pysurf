'''
This function defines the child mesh object for TFI, and also other auxiliary functions

Ney Secco 2016-08
'''


import numpy as np
from mpi4py import MPI
from . import tfiMeshTools
from pysurf import plot3d_interface

# OPTIONS

# Tolerance to state that two nodes are equal
tol = 1e-2

# Create a child class from the parent Mesh class defined in meshClass.py
class TFIMesh(object):

    def __init__(self, leftCurve, bottomCurve, rightCurve, topCurve, refComp):

        '''
        This function will receive four curve objects/array and
        one component object to initilize the TFI mesh.

        The curves should be ordered and directed as shown in the scheme below:

        +--> j,v
        |   +-----------<-------+
        V   |1       top       4|
        i,u |                   ^
            |left          right|
            V                   |
            |2      bottom     3|
            +------->-----------+

        You can also use the function linkcurves defind in this same file to automatically
        order any set of curves acording to this criteria. Note however that you will not be able
        to amp derivatives anymore.

        Instead of Curve objects, the user can also provide
        arrays of size [numNodes,3] for each curve.
        '''

        # CHECK INPUTS

        # Check if curves are connected
        badcurves = False

        if np.linalg.norm(leftCurve[-1,:] - bottomCurve[0,:]) > tol:
            badcurves = True
        elif np.linalg.norm(bottomCurve[-1,:] - rightCurve[0,:]) > tol:
            badcurves = True
        elif np.linalg.norm(rightCurve[-1,:] - topCurve[0,:]) > tol:
            badcurves = True
        elif np.linalg.norm(topCurve[-1,:] - leftCurve[0,:]) > tol:
            badcurves = True

        if badcurves:
            print('TFI Mesh ERROR:')
            print('curves are disconnected or are not following the TFIMesh convention.')
            print('Check the documentation in the tfi_mesh.py file, or use the function')
            print('linkcurves (defined in the same file) to get appropriate curves.')

        # INITIALIZE CLASS

        self.leftCurve = leftCurve
        self.bottomCurve = bottomCurve
        self.rightCurve = rightCurve
        self.topCurve = topCurve
        self.refComp = refComp

    def generate_mesh(self):
        '''
        This function actually computes the mesh.

        OUTPUTS:
        This function has no explicit outputs. It updates self.coor instead.
        '''

        # EXTRACT COORDINATES

        # Left curve
        if type(self.leftCurve) is np.ndarray:
            Xleft = self.leftCurve[:,0]
            Yleft = self.leftCurve[:,1]
            Zleft = self.leftCurve[:,2]
        else: # We have a curve object
            coor = self.leftCurve.extract_points()
            Xleft = coor[:,0]
            Yleft = coor[:,1]
            Zleft = coor[:,2]

        # Bottom curve
        if type(self.bottomCurve) is np.ndarray:
            Xbottom = self.bottomCurve[:,0]
            Ybottom = self.bottomCurve[:,1]
            Zbottom = self.bottomCurve[:,2]
        else: # We have a curve object
            coor = self.bottomCurve.extract_points()
            Xbottom = coor[:,0]
            Ybottom = coor[:,1]
            Zbottom = coor[:,2]

        # Right curve
        if type(self.rightCurve) is np.ndarray:
            Xright = self.rightCurve[:,0]
            Yright = self.rightCurve[:,1]
            Zright = self.rightCurve[:,2]
        else: # We have a curve object
            coor = self.rightCurve.extract_points()
            Xright = coor[:,0]
            Yright = coor[:,1]
            Zright = coor[:,2]

        # Top curve
        if type(self.topCurve) is np.ndarray:
            Xtop = self.topCurve[:,0]
            Ytop = self.topCurve[:,1]
            Ztop = self.topCurve[:,2]
        else: # We have a curve object
            coor = self.topCurve.extract_points()
            Xtop = coor[:,0]
            Ytop = coor[:,1]
            Ztop = coor[:,2]

        # FILL COORDINATES MATRIX

        # Now we need to fill all [nu x nv] arrays that contains the mesh coordinates.
        # But first we get the mesh size
        nu = len(Xleft)
        nv = len(Xbottom)

        print('nu')
        print(nu)
        print('nv')
        print(nv)
        print(Xright.shape)

        # Initialize matrices
        X = np.zeros((nu, nv))
        Y = np.zeros((nu, nv))
        Z = np.zeros((nu, nv))

        # Now we need to fill the border values

        # Left curve
        X[:,0] = Xleft
        Y[:,0] = Yleft
        Z[:,0] = Zleft

        # Bottom curve
        X[-1,:] = Xbottom
        Y[-1,:] = Ybottom
        Z[-1,:] = Zbottom

        # Right curve
        X[:,-1] = Xright[::-1]
        Y[:,-1] = Yright[::-1]
        Z[:,-1] = Zright[::-1]

        # Top curve
        X[0,:] = Xtop[::-1]
        Y[0,:] = Ytop[::-1]
        Z[0,:] = Ztop[::-1]

        # GENERATE MESH USING TFI
        # This will call the TFI functions and fill the interior values of X, Y, and Z
        tfiMeshTools.computeMesh(X,Y,Z)

        # PROJECTION
        # Now we need to project all points back to the surface
        self._project_mesh(X,Y,Z)

        # UPDATE COORDINATES
        self.coor = np.zeros((3,nu,nv))
        self.coor[0,:,:] = X
        self.coor[1,:,:] = Y
        self.coor[2,:,:] = Z

    def export(self, filename, addcurves=False):

        '''
        This function exports a 3D mesh in plot3d format.
        The user specifies the span and number of nodes in the z direction.

        INPUTS
        addcurves: logical -> Export separate zones for bounding curves if True
        '''

        # Get coordiantes
        X = self.coor[0,:,:].T
        Y = self.coor[1,:,:].T
        Z = self.coor[2,:,:].T

        # Initialize grid object
        myGrid = plot3d_interface.Grid()

        # Expand the coordinate matrices
        X3d = np.array([X])
        Y3d = np.array([Y])
        Z3d = np.array([Z])

        # Add mesh to the grid
        myGrid.add_block(X3d, Y3d, Z3d)

        if addcurves:
            # Add initial curves as a separate zone
            myGrid.add_block(np.array([[self.leftCurve[:,0]]]),
                             np.array([[self.leftCurve[:,1]]]),
                             np.array([[self.leftCurve[:,2]]]))
            myGrid.add_block(np.array([[self.bottomCurve[:,0]]]),
                             np.array([[self.bottomCurve[:,1]]]),
                             np.array([[self.bottomCurve[:,2]]]))
            myGrid.add_block(np.array([[self.rightCurve[:,0]]]),
                             np.array([[self.rightCurve[:,1]]]),
                             np.array([[self.rightCurve[:,2]]]))
            myGrid.add_block(np.array([[self.topCurve[:,0]]]),
                             np.array([[self.topCurve[:,1]]]),
                             np.array([[self.topCurve[:,2]]]))

        # Export grid
        plot3d_interface.export_plot3d(myGrid, filename)

    # ==================================================================

    def _project_mesh(self,X,Y,Z):

        '''
        This function will project the mesh coordinates into refComp.

        OUTPUTS:
        This function has no explicit outputs. It modifies X, Y, and Z instead.

        Ney Secco 2016-08
        '''

        # Get matrix sizes
        nu, nv = X.shape

        # Flatten matrices to get an array of points
        coor = np.vstack([X.flatten(),
                          Y.flatten(),
                          Z.flatten()]).T

        # Project all points
        coor, normals, projDict = self.refComp.project_on_surface(coor)

        # Now we replace the old coordinates
        X[:,:] = coor[:,0].reshape((nu,nv))
        Y[:,:] = coor[:,1].reshape((nu,nv))
        Z[:,:] = coor[:,2].reshape((nu,nv))

# ==================================================================
# ==================================================================
# ==================================================================

# ==================================================================
# AUXILIARY FUNCTIONS
# ==================================================================

def link_curves(curveList):

    '''
    This function will receive four curve objects/arrays and
    check if they can be concatenated. Then it returns
    the another set of four curves which can be readily used to initialize
    a TFI mesh object.

    Instead of Curve objects, the user can also provide
    arrays of size [numNodes,3] for each curve.

    INPUTS:
    curveList: list[4] -> List that should contain exactly for curve objects.
            Each element can also be an array[numNodes,3] specifying nodal coordinates.

    OUTPUTS:
    leftCurve: array[numNodes,3] -> coordinates of the left curve according to
               the TFIMesh Class orientation.
    bottomCurve: array[numNodes,3] -> coordinates of the bottom curve according to
                 the TFIMesh Class orientation.
    rightCurve: array[numNodes,3] -> coordinates of the right curve according to
                the TFIMesh Class orientation.
    topCurve: array[numNodes,3] -> coordinates of the top curve according to
              the TFIMesh Class orientation.

    Ney Secco 2016-08
    '''

    # CHECKING INPUTS

    # Initialize list that will have nodal coordinates for all 4 curves
    allcurves = []

    # Extract coordinates if the user provided curve Classes
    # rather than points
    for curve in curveList:

        # Check if we need to extract nodal coordinates
        if type(curve) is not np.ndarray:
            curve = curve.extract_points()

        # Store coordinates
        allcurves.append(curve)

    # We need to check which curves can be concatenated so that we can
    # form an ordered set of curves just as shown below:
    #
    #  +--> j,v
    #  |  +-------------------+
    #  V  |1       top       4|
    # i,u |                   |
    #     |left          right|
    #     |                   |
    #     |2      bottom     3|
    #     +-------------------+
    #
    # The numbers indicate the corner IDs.
    #
    # We will check the corner nodes for that.

    # Initialize curve variables so we can check if curves are not connected later on
    leftCurve = None
    bottomCurve = None
    rightCurve = None
    topCurve = None

    # We arbitralily select the first point of the first curve as the first corner,
    # and the last point of the first curve as the second corner. Thus the first curve
    # corresponds to the left curve.
    corner1 = allcurves[0][0,:]
    corner2 = allcurves[0][-1,:]
    leftCurve = allcurves.pop(0) # This will remove curve1 from allcurves

    # Now we find which curve is connected to corner2
    for curveID in range(len(allcurves)):

        # Get first and last points of the current curve:
        curve = allcurves[curveID]
        startPt = curve[0,:]
        endPt = curve[-1,:]

        # Check if the curve is connected and flip its
        # direction if necessary.
        if np.linalg.norm(corner2 - startPt) < tol:
            corner3 = endPt
            bottomCurve = allcurves.pop(curveID)
            break
        elif np.linalg.norm(corner2 - endPt) < tol:
            corner3 = startPt
            bottomCurve = allcurves.pop(curveID)[::-1]
            break

    # Stop if we did not find any curve
    if bottomCurve is None:
        print('TFI Mesh ERROR:')
        print('curves are not connected.')
        quit()

    # Now we find which curve is connected to corner3
    for curveID in range(len(allcurves)):

        # Get first and last points of the current curve:
        curve = allcurves[curveID]
        startPt = curve[0,:]
        endPt = curve[-1,:]

        # Check if the curve is connected and flip its
        # direction if necessary.
        if np.linalg.norm(corner3 - startPt) < tol:
            corner4 = endPt
            rightCurve = allcurves.pop(curveID)
            break
        elif np.linalg.norm(corner3 - endPt) < tol:
            corner4 = startPt
            rightCurve = allcurves.pop(curveID)[::-1]
            break

    # Stop if we did not find any curve
    if rightCurve is None:
        print('TFI Mesh ERROR:')
        print('curves are not connected.')
        quit()

    # The last curve should be connected to corner 4
    curve = allcurves[0]
    startPt = curve[0,:]
    endPt = curve[-1,:]

    # Check if the curve is connected and flip its
    # direction if necessary.
    if max(np.linalg.norm(corner4 - startPt), np.linalg.norm(corner1 - endPt)) < tol:
        topCurve = allcurves.pop(0)
    elif max(np.linalg.norm(corner4 - endPt), np.linalg.norm(corner1 - startPt)) < tol:
        topCurve = allcurves.pop(0)[::-1]
    else:
        print('TFI Mesh ERROR:')
        print('The four given curves do not define a closed quad.')
        quit()

    return leftCurve, bottomCurve, rightCurve, topCurve

# ==================================================================
