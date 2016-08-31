from __future__ import division
import os
import numpy as np
from mpi4py import MPI
from ...baseClasses import Component, Curve
from ....utilities import plot3d_interface
import utilitiesAPI, curveSearch
import tsurf_tools as tst
import adtAPI
import copy

class TSurfComponent(Component):

    def _initialize(self, *arg):

        '''
        This function initializes and TSurfComponent.
        It is called by the __init__ method of the parent Component
        class defined in classes.py.

        The expected arguments for the initialization function are:
        TSurfComponent(fileName, sectionsList, comm)

        REQUIRED INPUTS:
        fileName: string -> Name of the CGNS file that contains the
                  triangulated surface definition.

        OPTIONAL INPUTS:
        sectionsList: list of strings -> List of strings containing
                  the names of the sections in the CGNS file that
                  should be included in the current ADTComponent.
                  If nothing is provided, or if sectionList is None
                  or an empty list, then all sections will be included.

        comm: MPI communicator -> An MPI communicator, such as
              MPI.COMM_WORLD
        '''

        # Set dummy value to filename, so we can check later on if
        # the user actually gave an input name.
        filename = None

        # Set default values in case we have no additional arguments
        selectedSections = None # We will update this later if necessary
        self.comm = MPI.COMM_WORLD # Communicate to all processors

        # Check the optional arguments and do the necessary changes
        for optarg in arg:

            if type(optarg) == MPI.Intracomm: # We have an MPI Communicator
                self.comm = optarg

            elif type(optarg) == list:
                selectedSections = optarg

            elif optarg in (None, []):
                print 'Reading all CGNS sections in ADTComponent assigment.'

            elif type(optarg) == str:
                # Get filename
                filename = optarg

        # Check if the user provided no input file
        if filename == None:
            print ' ERROR: Cannot initialize TSurf Component as no input file'
            print ' was specified.'
            quit()

        # Read CGNS file
        self.coor, sectionDict = tst.getCGNSsections(filename, self.comm)
        self.name = os.path.splitext(os.path.basename(filename))[0]

        # Select all section names in case the user provided none
        if selectedSections is None:
            selectedSections = sectionDict.keys()
        else:
            self.name = self.name + "__" + "_".join(selectedSections)

        # Now we call an auxiliary function to merge selected surface sections in a single
        # connectivity array
        self.triaConn, self.quadsConn = tst.merge_surface_sections(sectionDict, selectedSections)

        # Initialize curves
        tst.initialize_curves(self, sectionDict, selectedSections)

        # Now we remove unused points
        self.coor = tst.remove_unused_points(self.coor, triaConn=self.triaConn, quadsConn=self.quadsConn)

        # Create ADT for the current surface, now that coor, triaConn, and quadsConn are correct
        tst.initialize_surface(self)

    def add_curve(self, name, curve):

        '''
        Adds a given curve instance to the self.Curves dictionary.
        '''

        self.Curves[name] = copy.deepcopy(curve)

    def update(self, coor):

        '''
        This function updates the nodal coordinates used by both surface and curve objects.
        '''

        # Update coordinates
        self.coor = coor

        # Update surface definition
        tst.update_surface(self)


    def translate(self, x, y, z):
        tst.translate(self, x, y, z)
        tst.update_surface(self)
        for curve in self.Curves.itervalues():
            curve.translate(x, y, z)

    def scale(self, factor):
        tst.scale(self, factor)
        tst.update_surface(self)

    def rotate(self, angle, axis):
        tst.rotate(self, angle, axis)
        tst.update_surface(self)

    def project_on_surface(self, xyz):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        OUTPUTS:
        xyzProj -> float[numPts,3] : Coordinates of the projected points

        normProj -> float[numPts,3] : Surface normal at projected points
        '''

        '''
        Explanation of reference values:


        xyzProj -> float[numPts, 3] : If the user already have previous projection candidates, they
                                      should be provided in this array. The code will only replace
                                      values if it finds a better candidade. If the user has no
                                      previous candidate, initialize all elements to zero and also
                                      set large values to dist2.

        normProj -> float[numPts, 3] : Surface normals computed at the user-provided projection candidates.
                                       If the user has no previous candidate, initialize all elements
                                       to zero and also set large values to dist2.

        dist2 -> float[numPts] : distance**2 of the user-provided projection candidates. The code
                                 will use this distance value to check for the best candidate and
                                 then replace values in xyzProj, and normProj acordingly. If no previous
                                 candidates are available, set all elements to a large number (1e10) so that
                                 all information in xyzProj and normProj is replaced.
        '''

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        normProj = np.zeros((numPts,3))

        # Call projection function
        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(xyz.T, self.name,
                                                                                 dist2, xyzProj.T,
                                                                                 self.nodal_normals, normProj.T)

        # Return projections
        return xyzProj, normProj

    def project_on_curve(self, xyz, curveCandidates=None):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        OUTPUTS:
        xyzProj -> float[numPts,3] : Coordinates of the projected points

        tanProj -> float[numPts,3] : Curve tangent at projected points
        '''

        '''
        Explanation of reference values:

        componentsList -> list of strings : Names of the surfaces components on which we should
                                            look for projections. If nothing is provided, the
                                            code will use all available surfaces.

        xyzProj -> float[numPts, 3] : If the user already have previous projection candidates, they
                                      should be provided in this array. The code will only replace
                                      values if it finds a better candidade. If the user has no
                                      previous candidate, initialize all elements to zero and also
                                      set large values to dist2.

        tanProj -> float[numPts, 3] : Curve tangent computed at the user-provided projection candidates.
                                      If the user has no previous candidate, initialize all elements
                                      to zero and also set large values to dist2.

        dist2 -> float[numPts] : distance**2 of the user-provided projection candidates. The code
                                 will use this distance value to check for the best candidate and
                                 then replace values in xyzProj, and normProj acordingly. If no previous
                                 candidates are available, set all elements to a large number (1e10) so that
                                 all information in xyzProj and normProj is replaced.
        '''

        # Use all curves if None is provided by the user
        if curveCandidates is None:
            curveCandidates = self.Curves.keys()

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        tanProj = np.zeros((numPts,3))

        # Check if the candidates are actually defined
        curveKeys = self.Curves.keys()
        for curve in curveCandidates:
            if curve not in curveKeys:
                print 'ERROR: Curve',curve,'is not defined. Check the curve names in your CGNS file.'
                print '       Also check if you included this curve name when selecting CGNS sections.'
                quit()

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.Curves:
            if curveName in curveCandidates:
                self.Curves[curveName].project(xyz, dist2, xyzProj, tanProj)

        # Return projections
        return xyzProj, tanProj

    def extract_curves(self, feature='sharpness'):
        '''
        This function will define new curves in this component based on
        surface features. The new curves will be stored in self.Curves
        '''
        tst.extract_curves_from_surface(self, feature)

    def intersect(self, otherComponent, distTol=1e-7, flip=False):
        '''
        This method will intersect the current component (self) with the provided
        TSurfComponent. The new curves will be stored in each component.

        INPUTS:
        otherComponent: Component object -> Other object that we want to intersect with.

        distTol: float -> Tolerance used to merge close nodes when joining bar elements generated
                 by the intersection.

        OUTPUTS:
        This function has no explicit outputs. However, the new curves will be added to
        self.Curves and otherComponent.Curves.
        '''

        # Call the intersection function defined in tsurf_tools.py
        Intersection = tst._compute_pair_intersection(self, otherComponent, distTol)

        # Add new curves to both components
        for curve in Intersection:

            # Flip the curve direction if necessary
            if flip:
                curve.flip()

            # Add curve to the current component
            self.add_curve(curve.name,curve)
            
            # Flip the curve
            curve.flip()

            # Add flipped curve to the second component
            # We do this so that the mesh grows in the correct direction when
            # calling hypsurf.
            otherComponent.add_curve(curve.name,curve)

        # Print the number of intersections
        print 'Number of intersections between',self.name,'and',otherComponent.name,'is',len(Intersection)

#=============================================================

class TSurfCurve(Curve):

    '''
    This is a derived class from the parent Curve class defined in
    baseClasses.py
    '''

    def _initialize(self, *arg): #,coor, barsConn, name, mergeTol=1e-2):

        '''
        This method initializes a TSurf Curve object, which uses finite element
        data to represent curves.
        This method will be called by the __init__ method of the parent Curve class

        REQUIRED INPUTS:
        name: string -> Curve name

        coor : array[nNodes,3] -> Nodal X,Y,Z coordinates.

        barsConn : array[nBars,2] -> Element connectivity matrix.

        OPTIONAL INPUTS:
        mergeTol: float -> Tolerance to merge nodes.

        John Jasa 2016-08
        Ney Secco 2016-08
        '''

        # Optional inputs
        mergeTol = 1e-4

        # Parse inputs
        for currArg in arg:

            if type(currArg) == str:
                # This is probably the name!
                name = currArg

            elif type(currArg) == np.ndarray:
                # This is an array.
                # We still need to figure out if we have coordinates or
                # bar connectivities

                # Get dimension of the current array
                dim = currArg.shape[0]

                if dim == 3:
                    # We have nodal coordinates. We need to make sure they
                    # are defined as floats, otherwise Fortran will not recognize them.
                    coor = np.array(currArg, dtype=float, order='F')

                elif dim == 2:
                    # We have bar connectivities. We need to make sure they
                    # are defined as integers, otherwise Fortran will not recognize them.
                    barsConn = np.array(currArg, dtype=np.int32, order='F')

                else:
                    print ' ERROR: Could not recognize array inputs when initializing TSurfCurve.'
                    print ' Please provide an [3 x n] array with nodal coordinates and an'
                    print ' [2 x m] array for bar connectivities.'
                    quit()

            elif type(currArg) == float:
                # This is probably mergeTol
                mergeTol = currArg

            else:
                    print ' ERROR: Could not recognize inputs when initializing TSurfCurve.'
                    quit()
 
        # Remove unused points (This will update coor and barsConn)
        coor = tst.remove_unused_points(coor, barsConn=barsConn)

        # Call a Fortran code to merge close nodes (This will update coor and barsConn)
        nUniqueNodes = utilitiesAPI.utilitiesapi.condensebarnodes(mergeTol, coor, barsConn)

        # Take bar connectivity and reorder it
        sortedConn = tst.FEsort(barsConn.T.tolist())

        # Check for failures
        if len(sortedConn)==1:
            # The sorting works just fine and it found a single curve!
            # So we just store this single curve (sortedConn[0]).
            sortedConn = np.array(sortedConn[0], dtype=type(barsConn[0,0]), order='F')
        else:
            # Just use the original connectivity
            sortedConn = barsConn
            # Print message
            print ''
            print 'Curve','"'+name+'"','could not be sorted. It might be composed by disconnect curves.'
            print ''

        # Assing coor and barsConn. Remember to crop coor to get unique nodes only
        self.coor = coor[:,:nUniqueNodes]
        self.barsConn = sortedConn
        self.name = name

    def extract_points(self):

        '''
        We will just return the nodes that define the curve.
        If the user wants another set of points, then he should use self.remesh first.
        '''

        # Get the number of elements
        numElems = self.barsConn.shape[1]

        # Initialize coordinates matrix
        pts = np.zeros((numElems+1, 3))

        # Get coordinates
        for ii in range(numElems):
            pts[ii,:] = self.coor[:, self.barsConn[0,ii]-1]

        # Get the last point
        pts[-1,:] = self.coor[:, self.barsConn[-1,-1]-1]

        # Return coordinates
        return pts

    def update_dvs(self, coor):
        self.coor = coor

    def flip(self):

        # Flip elements
        self.barsConn = self.barsConn[::-1]
        # Flip nodes
        for ii in range(len(self.barsConn)):
            self.barsConn[ii] = self.barsConn[ii][::-1]

    def translate(self, x, y, z):
        tst.translate(self, x, y, z)

    def scale(self, factor):
        tst.scale(self, factor)

    def rotate(self, angle, axis):
        tst.rotate(self, angle, axis)


    def project(self, xyz, dist2=None, xyzProj=None, tangents=None):

        '''
        This function will take the points given in xyz and project them to the surface.

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected

        dist2 -> float[nPoints] : values of best distance**2 found so far. If the distance of the projected
                                  point is less than dist2, then we will take the projected point. This allows us to find the best
                                  projection even with multiple surfaces.
                                  If you don't have previous values to dist2, just initialize all elements to
                                  a huge number (1e10).

        xyzProj -> float[nPoints,3] : coordinates of projected points found so far. These projections could be on other curves as well. This function will only replace projections whose dist2 are smaller
                                      than previous dist2. This allows us to use the same array while working with multiple curves.
                                      If you don't have previous values, just initialize all elements to zero. Also
                                      remember to set dist2 to a huge number so that all values are replaced.

        tangents -> float[nPoints,3] : tangent directions for the curve at the projected points.

        This function has no explicit outputs. It will just update dist2, xyzProj, and tangents
        '''

        # Get number of points
        nPoints = xyz.shape[0]

        # Initialize references if user provided none
        if dist2 == None:
            dist2 = np.ones(nPoints)*1e10

        if xyzProj==None:
            xyzProj = np.zeros((nPoints,3))

        if tangents==None:
            tangents = np.zeros((nPoints,3))

        # Call fortran code
        curveSearch.curveproj.mindistancecurve(xyz.T, self.coor, self.barsConn, xyzProj.T, tangents.T, dist2)

    def remesh(self, nNewNodes=None, method='linear', spacing='linear'):

        '''
        This function will redistribute the nodes along the curve to get
        better spacing among the nodes.
        This assumes that the FE data is ordered.

        The first and last nodes will remain at the same place.

        Consider using self.shift_end_nodes first so you can preserve the start and end points
        of a periodic curve.

        INPUTS:

        nNewNodes: integer -> Number of new nodes desired in the new curve definition.

        method: string -> Method used to interpolate new nodes with respect to
                existing ones. Check scipy.interpolate.interp1d for options.

        spacing: string -> Desired spacing criteria for new nodes. Current options are:
                 ['linear']

        OUTPUTS:
        This method has no explicit outputs. It will update self.coor and self.barsConn instead.

        Ney Secco 2016-08
        '''

        # Get connectivities and coordinates of the current Curve object
        coor = self.coor
        barsConn = self.barsConn

        # Get the number of nodes and elements in the curve
        nElem = barsConn.shape[1]
        nNodes = nElem + 1

        # CHECKING INPUTS

        # First we check if the FE data is ordered
        for elemID in range(2,nElem):

            # Get node indices
            prevNodeID = barsConn[1,elemID-1]
            currNodeID = barsConn[0,elemID]

            # Check if the FE data is ordered
            if prevNodeID != currNodeID:

                # Print warning
                print 'WARNING: Could not remesh curve because it has unordered FE data.'
                print '         Call FEsort first.'
                return

        # COMPUTE ARC-LENGTH

        # We can proceed if FE data is ordered

        # Initialize an array that will store the arclength of each node.
        # That is, the distance, along the curve from the current node to
        # the first node of the curve.
        arcLength = np.zeros(nNodes)

        # Also initialize an array to store nodal coordinates in the curve order
        nodeCoor = np.zeros((3,nNodes))

        # Store position of the first node (the other nodes will be covered in the loop)
        # (the -1 is due Fortran indexing)
        nodeCoor[:,0] = coor[:,barsConn[0,0]-1]

        # Loop over each element to increment arcLength
        for elemID in range(nElem):

            # Get node positions (the -1 is due Fortran indexing)
            node1 = coor[:,barsConn[0,elemID]-1]
            node2 = coor[:,barsConn[1,elemID]-1]

            # Compute distance between nodes
            dist = np.linalg.norm(node1 - node2)

            # Store nodal arc-length
            arcLength[elemID+1] = arcLength[elemID] + dist

            # Store coordinates of the next node
            nodeCoor[:,elemID+1] = node2

        # SAMPLING POSITION FOR NEW NODES

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        # The arcLength will be used as our parametric coordinate to sample the new nodes.

        # First we will initialize arrays for the new coordinates and arclengths
        newNodeCoor = np.zeros((3,nNewNodes), order='F')
        newArcLength = np.zeros(nNewNodes, order='F')

        # Now that we know the initial and final arcLength, we can redistribute the
        # parametric coordinates based on the used defined spacing criteria.
        # These statements should initially create parametric coordinates in the interval
        # [0.0, 1.0]. We will rescale it after the if statements.
        if spacing == 'linear':
            newArcLength = np.linspace(0.0, 1.0, nNewNodes)

        elif spacing == 'cosine':
            newArcLength = 0.5*(1.0 - np.cos(np.linspace(0.0, np.pi, nNewNodes)))

        # Rescale newArcLength based on the final distance
        newArcLength = arcLength[-1]*newArcLength

        # INTERPOLATE NEW NODES

        # Now we sample the new coordinates based on the interpolation method given by the user

        # Import interpolation function
        from scipy.interpolate import interp1d

        # Create interpolants for x, y, and z
        fX = interp1d(arcLength, nodeCoor[0,:], kind=method)
        fY = interp1d(arcLength, nodeCoor[1,:], kind=method)
        fZ = interp1d(arcLength, nodeCoor[2,:], kind=method)

        # Sample new points using the interpolation functions
        newNodeCoor[0,:] = fX(newArcLength)
        newNodeCoor[1,:] = fY(newArcLength)
        newNodeCoor[2,:] = fZ(newArcLength)

        # ASSIGN NEW COORDINATES AND CONNECTIVITIES

        # Set new nodes
        self.coor = newNodeCoor

        # Check if the baseline curve is periodic
        if barsConn[0,0] == barsConn[1,-1]:
            periodic = True
        else:
            periodic = False

        # Generate new connectivity (the nodes are already in order so we just
        # need to assign an ordered set to barsConn).
        barsConn = np.zeros((2,nNewNodes-1), order='F', dtype=int)
        barsConn[0,:] = range(1,nNewNodes)
        barsConn[1,:] = range(2,nNewNodes+1)

        # We still need to keep periodicity if the original curve is periodic
        if periodic:
            barsConn[1,-1] = barsConn[0,0]

        self.barsConn = barsConn

    def shift_end_nodes(self, criteria='maxX'):

        '''
        This method will shift the finite element ordering of
        a periodic curve so that the first and last elements
        satisfy a given criteria.
        This is useful when we want to set boundary conditions for
        mesh generation in specific points of this periodic curve,
        such as a trailing edge.

        This method only works for periodic curves with ordered FE.

        INPUTS:

        criteria: string -> criteria used to reorder FE data. Available options are:
                  ['maxX', 'maxY', 'maxZ']

        OUTPUTS:
        This function has no explicit outputs. It changes self.barsConn.

        Ney Secco 2016-08
        '''

        # Get current coordinates and connectivities
        coor = self.coor
        barsConn = self.barsConn

        # Get number of elements
        nElem = barsConn.shape[1]

        # CHECKING INPUTS

        # First we check if the FE data is ordered and periodic
        for elemID in range(2,nElem):

            # Get node indices
            prevNodeID = barsConn[1,elemID-1]
            currNodeID = barsConn[0,elemID]

            # Check if the FE data is ordered
            if prevNodeID != currNodeID:

                # Print warning
                print 'WARNING: Could not shift curve because it has unordered FE data.'
                print '         Call FEsort first.'
                return

        # The first and last nodes should be the same to have a periodic curve.
        if barsConn[0,0] != barsConn[1,-1]:

            # Print warning
            print 'WARNING: Could not shift curve because it is not periodic.'
            return

        # REORDERING

        # Reorder according to the given criteria. We need to find the
        # reference node (startNodeID) to reorder the FE data

        if criteria in ['maxX','maxY','maxZ']:

            if criteria == 'maxX':
                coorID = 0
            elif criteria == 'maxY':
                coorID = 1
            elif criteria == 'maxZ':
                coorID = 2

            # Get maximum value of the desired coordinate.
            # The +1 is to make it consistent with Fortran ordering.
            startNodeID = np.argmax(coor[coorID,:]) + 1

        elif criteria in ['minX','minY','minZ']:

            if criteria == 'minX':
                coorID = 0
            elif criteria == 'minY':
                coorID = 1
            elif criteria == 'minZ':
                coorID = 2

            # Get maximum value of the desired coordinate.
            # The +1 is to make it consistent with Fortran ordering.
            startNodeID = np.argmin(coor[coorID,:]) + 1

        # Now we look for an element that starts at the reference node
        startElemID = np.where(barsConn[0,:] == startNodeID)[0][0]

        # Then we shift the connectivity list so that startElemID becomes the first element
        self.barsConn[:,:] = np.hstack([barsConn[:,startElemID:], barsConn[:,:startElemID]])

    def export_plot3d(self,outputName='curve'):

        '''
        This function will export the current curve in plot3d format

        Ney Secco 2016-08
        '''

        # Create plot3d curve object and export
        p3dCurve = plot3d_interface.Curve(self.coor, self.barsConn)
        p3dCurve.export_plot3d(outputName)
