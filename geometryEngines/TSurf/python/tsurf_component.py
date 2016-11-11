from __future__ import division
import os
import numpy as np
from mpi4py import MPI
from ...baseClasses import Geometry, Curve
from ....utilities import plot3d_interface
import utilitiesAPI, curveSearchAPI
import tsurf_tools as tst
import adtAPI
import copy

fortran_flag = True

class TSurfGeometry(Geometry):

    def _initialize(self, *arg, **kwargs):

        '''
        This function initializes and TSurfGeometry.
        It is called by the __init__ method of the parent Geometry
        class defined in classes.py.

        The expected arguments for the initialization function are:
        TSurfGeometry(fileName, sectionsList, comm)

        REQUIRED INPUTS:
        fileName: string -> Name of the CGNS file that contains the
                  triangulated surface definition.

        OPTIONAL INPUTS:
        sectionsList: list of strings -> List of strings containing
                  the names of the sections in the CGNS file that
                  should be included in the current ADTGeometry.
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
                print 'Reading all CGNS sections in ADTGeometry assigment.'

            elif type(optarg) == str:
                # Get filename
                filename = optarg

        # Check if the user provided no input file
        if filename == None:
            print ' ERROR: Cannot initialize TSurf Geometry as no input file'
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
        Adds a given curve instance to the self.curves dictionary.
        '''

        self.curves[name] = copy.deepcopy(curve)

    def update(self, coor):

        '''
        This function updates the nodal coordinates used by both surface and curve objects.
        '''

        # First check if we have the same number of new coordinates
        if not self.coor.shape == coor.shape:
            print ''
            print 'WARNING: self.update in TSurfGeometry class'
            print '         The new set of coordinates does not have the'
            print '         same number of points as the original set.'
            print ''

        # Update coordinates
        self.coor = coor

        # Update surface definition
        tst.update_surface(self)


    def translate(self, x, y, z):
        tst.translate(self, x, y, z)
        tst.update_surface(self)
        for curve in self.curves.itervalues():
            curve.translate(x, y, z)

    def scale(self, factor):
        tst.scale(self, factor)
        tst.update_surface(self)
        for curve in self.curves.itervalues():
            curve.scale(factor)

    def rotate(self, angle, axis):
        tst.rotate(self, angle, axis)
        tst.update_surface(self)
        for curve in self.curves.itervalues():
            curve.rotate(angle, axis)

    def project_on_surface(self, xyz):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        OUTPUTS:
        xyzProj -> float[numPts,3] : Coordinates of the projected points

        normProj -> float[numPts,3] : Surface normal at projected points

        projDict -> dictionary : Dictionary containing intermediate values that are
                                 required by the differentiation routines.
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
        normProjNotNorm = np.zeros((numPts,3))

        # Call projection function
        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(xyz.T, self.name,
                                                                                 dist2, xyzProj.T,
                                                                                 self.nodal_normals, normProjNotNorm.T)

        # Store additional outputs in a dictionary to make outputs cleaner
        projDict = {'procID':procID,
                    'elementType':elementType,
                    'elementID':elementID,
                    'uvw':uvw,
                    'dist2':dist2,
                    'normProjNotNorm':normProjNotNorm}


        # Normalize the normals
        normProj = tst.normalize(normProjNotNorm)

        # Return projections
        return xyzProj, normProj, projDict

    def project_on_surface_d(self, xyz, xyzd, xyzProj, normProj, projDict, coord):

        '''
        This function will compute derivatives of the projection algorithm in forward mode.

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        xyzd -> float[numPts, 3] : Derivative seeds for coordinates of the points
                                   that should be projected.

        xyzProj -> float[numPts,3] : Coordinates of the projected points (can be obtained
                                     with the original projection function)

        normProj -> float[numPts,3] : Surface normal at projected points (can be obtained
                                      with the original projection function)

        projDict -> dictionary : Dictionary containing intermediate values that are
                                 required by the differentiation routines.  (can be obtained
                                 with the original projection function)

        coord -> float[numNodes,3] : Derivative seeds of the nodal coordinates of the baseline
                                     surface.

        OUTPUTS:
        xyzProjd -> float[numPts,3] : Derivative seeds of the coordinates of the projected points

        normProjd -> float[numPts,3] : Derivative seeds of the surface normal at projected points

        Ney Secco 2016-10
        '''

        # Rename variables to make code more readable
        procID = projDict['procID']
        elementType = projDict['elementType']
        elementID = projDict['elementID']
        uvw = projDict['uvw']
        dist2 = projDict['dist2']
        normProjNotNorm = projDict['normProjNotNorm']

        # Compute derivatives of the normal vectors
        nodal_normals, nodal_normalsd = adtAPI.adtapi.adtcomputenodalnormals_d(self.coor, coord,
                                                                               self.triaConn, self.quadsConn)

        # Call projection function
        # ATTENTION: The variable "xyz" here in Python corresponds to the variable "coor" in the Fortran code.
        # On the other hand, the variable "coor" here in Python corresponds to the variable "adtCoor" in Fortran.
        # I could not change this because the original ADT code already used "coor" to denote nodes that should be
        # projected.
        xyzProjd, normProjNotNormd = adtAPI.adtapi.adtmindistancesearch_d(xyz.T, xyzd.T, self.name, coord,
                                                                          procID, elementType,
                                                                          elementID, uvw,
                                                                          dist2, xyzProj.T, self.nodal_normals,
                                                                          nodal_normalsd, normProjNotNorm.T)

        # Transpose results to make them consistent
        xyzProjd = xyzProjd.T
        normProjNotNormd = normProjNotNormd.T


        # Now we need to compute derivatives of the normalization process
        normProj, normProjd = tst.normalize_d(normProjNotNorm, normProjNotNormd)

        # Return projections derivatives
        return xyzProjd, normProjd

    def project_on_surface_b(self, xyz, xyzProj, xyzProjb, normProj, normProjb, projDict):

        '''
        This function will compute derivatives of the projection algorithm in forward mode.

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        xyzProj -> float[numPts,3] : Coordinates of the projected points (can be obtained
                                     with the original projection function)

        xyzProjb -> float[numPts,3] : Derivative seeds of the coordinates of the projected points

        normProj -> float[numPts,3] : Surface normal at projected points (can be obtained
                                      with the original projection function)

        normProjb -> float[numPts,3] : Derivative seeds of the surface normal at projected points

        projDict -> dictionary : Dictionary containing intermediate values that are
                                 required by the differentiation routines. (can be obtained
                                 with the original projection function)

        OUTPUTS:

        xyzb -> float[numPts, 3] : Derivative seeds for coordinates of the points
                                   that should be projected.

        coorb -> float[numNodes,3] : Derivative seeds of the nodal coordinates of the baseline
                                     surface.

        Ney Secco 2016-10
        '''

        # Rename variables to make code more readable
        procID = projDict['procID']
        elementType = projDict['elementType']
        elementID = projDict['elementID']
        uvw = projDict['uvw']
        dist2 = projDict['dist2']
        normProjNotNorm = projDict['normProjNotNorm']

        # Compute derivatives of the normalization process
        normProjNotNormb = tst.normalize_b(normProjNotNorm, normProjb)

        # Call projection function
        # ATTENTION: The variable "xyz" here in Python corresponds to the variable "coor" in the Fortran code.
        # On the other hand, the variable "coor" here in Python corresponds to the variable "adtCoor" in Fortran.
        # I could not change this because the original ADT code already used "coor" to denote nodes that should be
        # projected.
        xyzb, coorb, nodal_normalsb = adtAPI.adtapi.adtmindistancesearch_b(xyz.T, self.name,
                                                                           procID, elementType,
                                                                           elementID, uvw,
                                                                           dist2, xyzProj.T,
                                                                           xyzProjb.T, self.nodal_normals,
                                                                           normProjNotNorm.T, normProjNotNormb.T)

        # Transpose results to make them consistent
        xyzb = xyzb.T

        # Compute derivatives of the normal vectors
        coorb = coorb + adtAPI.adtapi.adtcomputenodalnormals_b(self.coor,
                                                               self.triaConn, self.quadsConn,
                                                               self.nodal_normals, nodal_normalsb)

        # Return projections derivatives
        return xyzb, coorb

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
            curveCandidates = self.curves.keys()

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        tanProj = np.zeros((numPts,3))
        elemIDs = np.zeros((numPts),dtype='int32')

        # Check if the candidates are actually defined
        curveKeys = self.curves.keys()
        for curve in curveCandidates:
            if curve not in curveKeys:
                print 'ERROR: Curve',curve,'is not defined. Check the curve names in your CGNS file.'
                print '       Also check if you included this curve name when selecting CGNS sections.'
                quit()

        # Initialize list that will contain the names of the curves that got the best projections for
        # each point
        curveIDs = [0]*numPts

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:
            if curveName in curveCandidates:

                # Run projection code
                curveMask = self.curves[curveName].project(xyz, dist2, xyzProj, tanProj, elemIDs)

                # Use curveMask to update names of curves that got best projections
                for ii in range(numPts):
                    if curveMask[ii] == 1:
                        curveIDs[ii] = curveName

        # Create dictionary with useful outputs for differentiated codes
        curveProjDict = {}
        curveProjDict['elemIDs'] = elemIDs
        curveProjDict['curveIDs'] = curveIDs
        curveProjDict['curveCandidates'] = curveCandidates

        # Return projections
        return xyzProj, tanProj, curveProjDict

    def project_on_curve_d(self, xyz, xyzd, allCoord, xyzProj, tanProj, curveProjDict):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        xyzProj -> float[numPts,3] : Coordinates of the projected points

        tanProj -> float[numPts,3] : Curve tangent at projected points

        allCoord -> dictionary : This is a dictionary whose keys are all curve names. Each entry
                    should contain an array [nNodes x 3] with derivative seeds for all nodes of the curve

        Ney Secco 2016-11
        '''

        # Retrieve data from dictionary
        elemIDs = curveProjDict['elemIDs']
        curveIDs = curveProjDict['curveIDs']
        curveCandidates = curveProjDict['curveCandidates']

        # Get number of points
        numPts = xyz.shape[0]

        # Initialize derivatives
        xyzProjd = np.zeros((numPts,3))
        tanProjd = np.zeros((numPts,3))

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:
            if curveName in curveCandidates:

                # Identify points that are projected to the current curve to create a mask
                curveMask = [0]*numPts
                for ii in range(numPts):
                    if curveIDs[ii] == curveName:
                        curveMask[ii] = 1

                # Get derivative seeds for the nodal points of the current curve
                coord = allCoord[curveName]

                # Run projection code in forward mode
                self.curves[curveName].project_d(xyz, xyzd, coord, xyzProj, xyzProjd, tanProj, tanProjd, elemIDs, curveMask)

        # Return projections derivatives
        return xyzProjd, tanProjd

    def project_on_curve_b(self, xyz, xyzProj, xyzProjb, tanProj, tanProjb, curveProjDict):

        '''
        This function will backpropagate projections derivatives back
        to curve nodes and unprojected points.

        Ney Secco 2016-11
        '''

        # Retrieve data from dictionary
        elemIDs = curveProjDict['elemIDs']
        curveIDs = curveProjDict['curveIDs']
        curveCandidates = curveProjDict['curveCandidates']

        # Get number of points
        numPts = xyz.shape[0]

        # Initialize derivatives
        xyzb = np.zeros((numPts,3))
        allCoorb = {}

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:

            # Initialize derivatives
            coorb = np.zeros(self.curves[curveName].coor.shape, order='F')

            if curveName in curveCandidates:

                # Identify points that are projected to the current curve to create a mask
                curveMask = [0]*numPts
                for ii in range(numPts):
                    if curveIDs[ii] == curveName:
                        curveMask[ii] = 1

                # Run projection code in reverse mode
                # This will modify xyzb and coorb
                self.curves[curveName].project_b(xyz, xyzb, coorb, xyzProj, xyzProjb, tanProj, tanProjb, elemIDs, curveMask)

            # Store curve node derivatives in the corresponding dictionary
            allCoorb[curveName] = coorb

        # Return backpropagated derivatives
        return xyzb, allCoorb

    def extract_curves(self, feature='sharpness'):
        '''
        This function will define new curves in this component based on
        surface features. The new curves will be stored in self.curves
        '''
        tst.extract_curves_from_surface(self, feature)

    def intersect(self, otherGeometry, distTol=1e-7, flip=False):
        '''
        This method will intersect the current component (self) with the provided
        TSurfGeometry. The new curves will be stored in each component.

        INPUTS:
        otherGeometry: Geometry object -> Other object that we want to intersect with.

        distTol: float -> Tolerance used to merge close nodes when joining bar elements generated
                 by the intersection.

        OUTPUTS:
        This function has no explicit outputs. However, the new curves will be added to
        self.curves and otherGeometry.curves.
        '''

        # Call the intersection function defined in tsurf_tools.py
        Intersection = tst._compute_pair_intersection(self, otherGeometry, distTol)

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
            otherGeometry.add_curve(curve.name,curve)

        # Print the number of intersections
        print 'Number of intersections between',self.name,'and',otherGeometry.name,'is',len(Intersection)

#=============================================================
#=============================================================
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
        mergeTol = 1e-7

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
        sortedConn, dummy_map = tst.FEsort(barsConn.T.tolist())

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
        self.extra_data = {}

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
        # Flip the extra data
        for data in self.extra_data:
            self.extra_data[data] = self.extra_data[data][::-1]

    def translate(self, x, y, z):
        tst.translate(self, x, y, z)

    def scale(self, factor):
        tst.scale(self, factor)

    def rotate(self, angle, axis):
        tst.rotate(self, angle, axis)

    def update(self, coor):
        self.coor = coor

    def project(self, xyz, dist2=None, xyzProj=None, tangents=None, elemIDs=None):

        '''
        This function will take the points given in xyz and project them to the surface.

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected

        dist2 -> float[nPoints] : values of best distance**2 found so far. If the distance of the projected
                                  point is less than dist2, then we will take the projected point. This allows us to find the best
                                  projection even with multiple surfaces.
                                  If you don't have previous values to dist2, just initialize all elements to
                                  a huge number (1e10).

        xyzProj -> float[nPoints,3] : coordinates of projected points found so far. These projections could be on other
                                      curves as well. This function will only replace projections whose dist2 are smaller
                                      than previous dist2. This allows us to use the same array while working with multiple curves.
                                      If you don't have previous values, just initialize all elements to zero. Also
                                      remember to set dist2 to a huge number so that all values are replaced.

        tangents -> float[nPoints,3] : tangent directions for the curve at the projected points.

        elemIDs -> int[nPoints] : ID of the bar elements that received projections. This is also given by
                                  the execution of the primal routine.

        This function has no explicit outputs. It will just update dist2, xyzProj, tangents, and elemIDs
        '''

        # Get number of points
        nPoints = xyz.shape[0]

        # Initialize references if user provided none
        if dist2 is None:
            dist2 = np.ones(nPoints)*1e10

        if xyzProj is None:
            xyzProj = np.zeros((nPoints,3))

        if tangents is None:
            tangents = np.zeros((nPoints,3))

        if elemIDs is None:
            elemIDs = np.zeros((nPoints),dtype='int32')

        # Call fortran code
        curveMask = curveSearchAPI.curvesearchapi.mindistancecurve(xyz.T, self.coor, self.barsConn,
                                                                   xyzProj.T, tangents.T, dist2, elemIDs)

        # curveMask is an array of length nPoints. If point ii finds a better projection on this curve, then curveMask[ii]=1.
        # Otherwise, curveMask[ii]=0

        return curveMask


    def project_d(self, xyz, xyzd, coord, xyzProj, xyzProjd, tanProj, tanProjd, elemIDs, curveMask):

        '''
        This function will run forward mode AD to propagate derivatives from inputs (xyz and coor) to outputs (xyzProj).

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected.

        xyzd -> float[nPoints,3] : derivatives seeds of the coordinates.

        coord -> float[nNodes, 3] : derivative seeds of the nodes that constitutes the bar elements.

        xyzProj -> float[nPoints,3] : coordinates of projected points found with the primal routine.

        xyzProjd -> float[nPoints,3] : Derivative seeds of the projected points.

        elemIDs -> int[nPoints] : ID of the bar elements that received projections. This is also given by
                                  the execution of the primal routine.

        curveMask -> int[nPoints] : curveMask[ii] should be 1 if the ii-th point was actually projected onto
                                    this curve. Otherwise, curveMaks[ii] = 0, then the code will not compute
                                    derivatives for this point.

        OUTPUTS:

        This function has no explicit outputs. It will just update xyzProjd.

        Ney Secco 2016-11
        '''

        # Get number of points
        nPoints = xyz.shape[0]

        if self.coor.shape != coord.shape:
            print 'ERROR: Derivative seeds should have the same dimension of the original'
            print 'variable. The number of derivatives for the bar element nodes does not match'
            print 'with the number of nodes.'

        # Call fortran code
        curveSearchAPI.curvesearchapi.mindistancecurve_d(xyz.T, xyzd.T,
                                                         self.coor, coord,
                                                         self.barsConn,
                                                         xyzProj.T, xyzProjd.T,
                                                         tanProj.T, tanProjd.T,
                                                         elemIDs, curveMask)

    def project_b(self, xyz, xyzb, coorb, xyzProj, xyzProjb, tanProj, tanProjb, elemIDs, curveMask):

        '''
        This function will run forward mode AD to propagate derivatives from inputs (xyz and coor) to outputs (xyzProj).

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected.

        xyzb -> float[nPoints,3] : derivatives seeds of the coordinates.

        coorb -> float[nNodes, 3] : derivative seeds of the nodes that constitutes the bar elements.

        xyzProj -> float[nPoints,3] : coordinates of projected points found with the primal routine.

        xyzProjb -> float[nPoints,3] : Derivative seeds of the projected points.

        elemIDs -> int[nPoints] : ID of the bar elements that received projections. This is also given by
                                  the execution of the primal routine.

        curveMask -> int[nPoints] : curveMask[ii] should be 1 if the ii-th point was actually projected onto
                                    this curve. Otherwise, curveMaks[ii] = 0, then the code will not compute
                                    derivatives for this point.

        OUTPUTS:

        This function has no explicit outputs. It will just update xyzb and coorb.

        Ney Secco 2016-11
        '''

        # Get number of points
        nPoints = xyz.shape[0]

        if xyzProj.shape != xyzProjb.shape:
            print 'ERROR: Derivative seeds should have the same dimension of the original'
            print 'variable. The number of derivatives for the projected points does not match'
            print 'with the number of points.'

        # Call fortran code
        curveSearchAPI.curvesearchapi.mindistancecurve_b(xyz.T, xyzb.T,
                                                         self.coor, coorb,
                                                         self.barsConn,
                                                         xyzProj.T, xyzProjb.T,
                                                         tanProj.T, tanProjb.T,
                                                         elemIDs, curveMask)

    def remesh(self, nNewNodes=None, method='linear', spacing='linear',
               initialSpacing=0.1, finalSpacing=0.1):

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
                 ['linear', 'cosine', 'hypTan', 'tangent']

        initialSpacing: float -> Desired distance between the first and second nodes. Only
                                 used by 'hypTan' and 'tangent'.
        finalSpacing: float -> Desired distance between the last two nodes. Only
                               used by 'hypTan' and 'tangent'.

        OUTPUTS:
        This method has no explicit outputs. It will update self.coor and self.barsConn instead.

        Ney Secco 2016-08
        '''

        # Get connectivities and coordinates of the current Curve object
        coor = np.array(self.coor,dtype=type(self.coor[0,0]),order='F')
        barsConn = np.array(self.barsConn,dtype=type(self.barsConn[0,0]),order='F')

        # Get the number of elements in the curve
        nElem = barsConn.shape[1]
        nNodes = coor.shape[1]

        # Check if the baseline curve is periodic. If this is the case, we artificially repeat
        # the last point so that we could use the same code of the non-periodic case
        if barsConn[0,0] == barsConn[1,-1]:
            periodic = True
            coor = np.array(np.hstack([coor, coor[:,barsConn[0,0]-1].reshape((3,1))]),dtype=type(self.coor[0,0]),order='F')
            barsConn[-1,-1] = nNodes+1
        else:
            periodic = False

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        if fortran_flag:
            newCoor, newBarsConn = utilitiesAPI.utilitiesapi.remesh(nNewNodes, coor, barsConn, method, spacing)

        else:

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


            # The arcLength will be used as our parametric coordinate to sample the new nodes.

            # First we will initialize arrays for the new coordinates and arclengths
            newNodeCoor = np.zeros((3,nNewNodes), order='F')
            newArcLength = np.zeros(nNewNodes, order='F')

            # Now that we know the initial and final arcLength, we can redistribute the
            # parametric coordinates based on the used defined spacing criteria.
            # These statements should initially create parametric coordinates in the interval
            # [0.0, 1.0]. We will rescale it after the if statements.
            if spacing.lower() == 'linear':
                newArcLength = np.linspace(0.0, 1.0, nNewNodes)

            elif spacing.lower() == 'cosine':
                newArcLength = 0.5*(1.0 - np.cos(np.linspace(0.0, np.pi, nNewNodes)))

            elif spacing.lower() == 'hypTan':
                newArcLength = tst.hypTanDist(initialSpacing/arcLength[-1], finalSpacing/arcLength[-1], nNewNodes)

            elif spacing.lower() == 'tangent':
                newArcLength = tst.tanDist(initialSpacing/arcLength[-1], finalSpacing/arcLength[-1], nNewNodes)

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
            newCoor = newNodeCoor

            # Generate new connectivity (the nodes are already in order so we just
            # need to assign an ordered set to barsConn).
            barsConn = np.zeros((2,nNewNodes-1), order='F', dtype=int)
            barsConn[0,:] = range(1,nNewNodes)
            barsConn[1,:] = range(2,nNewNodes+1)

            newBarsConn = barsConn

        # Adjust connectivities if curve is periodic
        if periodic:
            newBarsConn[-1,-1] = newBarsConn[0,0]
            newCoor = newCoor[:,:-1]

        # Create a new curve object and return it
        # This way the original curve coordinates and connectivities remain the same
        newCurve = copy.deepcopy(self)
        newCurve.coor = newCoor
        newCurve.barsConn = newBarsConn

        return newCurve


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
