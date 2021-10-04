import os
import numpy as np
from mpi4py import MPI
from .baseClasses import Geometry, Curve
from . import tecplot_interface
from . import utilitiesAPI, curveSearchAPI
from . import tsurf_tools as tst
from . import adtAPI
import copy

fortran_flag = True


class TSurfGeometry(Geometry):
    def _initialize(self, *arg, **kwargs):

        """
        This function initializes and TSurfGeometry.
        It is called by the __init__ method of the parent Geometry
        class defined in classes.py.

        The expected arguments for the initialization function are:
        TSurfGeometry(fileName, sectionsList, comm, name=name)

        Parameters
        ----------
        fileName: string, optional
            Name of the CGNS file that contains the
            triangulated surface definition.

        sectionsList: list of strings, optional
            List of strings containing
            the names of the sections in the CGNS file that
            should be included in the current ADTGeometry.
            If nothing is provided, or if sectionList is None
            or an empty list, then all sections will be included.

        comm: MPI communicator, optional
            An MPI communicator, such as MPI.COMM_WORLD

        name: string, optional
            Name that can be assigned to the geometry. This name
            should match the wall boundary condition family names used by the solver.
        """

        # Set dummy value to filename, so we can check later on if
        # the user actually gave an input name.
        filename = None

        # Set default values in case we have no additional arguments
        selectedSections = None  # We will update this later if necessary

        # Check the optional arguments and do the necessary changes
        for optarg in arg:

            if type(optarg) == MPI.Intracomm:  # We have an MPI Communicator
                self.comm = optarg

            elif type(optarg) == list:
                selectedSections = optarg

            elif optarg in (None, []):
                print("Reading all CGNS sections in ADTGeometry assigment.")

            elif type(optarg) == str:
                # Get filename
                filename = optarg

        # SERIALIZATION
        # Here we create a new communicator just for the root proc because TSurf currently runs in
        # a single proc.
        comm = MPI.COMM_WORLD

        # Create colors so that only the root proc is the worker
        if comm.Get_rank() == 0:
            color = 0
        else:
            color = MPI.UNDEFINED

        newComm = comm.Split(color)

        self.comm = newComm

        # Assign component name based on its CGNS file
        self.name = os.path.splitext(os.path.basename(filename))[0]

        # We create a custom name if the user selected a subgroup of sections
        if selectedSections is not None:
            self.name = self.name + "__" + "_".join(selectedSections)

        # Rename everything if the user gives a name
        if "name" in kwargs:
            self.name = kwargs["name"]

        # Only the root proc will work here
        if self.myID == 0:

            # Check if the user provided no input file
            if filename is None:
                print(" ERROR: Cannot initialize TSurf Geometry as no input file")
                print(" was specified.")
                quit()

            # Get file extension
            fileExt = os.path.splitext(os.path.basename(filename))[1]

            if fileExt == ".cgns":

                # Read CGNS file
                self.coor, sectionDict = tst.getCGNSsections(filename, self.comm)

                # Select all section names in case the user provided none
                if selectedSections is None:
                    selectedSections = list(sectionDict.keys())

                # Now we call an auxiliary function to merge selected surface sections in a single
                # connectivity array
                self.triaConnF, self.quadsConnF = tst.merge_surface_sections(sectionDict, selectedSections)

                # Initialize curves
                tst.initialize_curves(self, sectionDict, selectedSections)

            elif fileExt == ".plt":

                # Read Tecplot file
                self.coor, triaConn, quadsConn = tecplot_interface.readTecplotFEdataSurf(filename)

                # Adjust indices to Fortran ordering
                self.triaConnF = triaConn + 1
                self.quadsConnF = quadsConn + 1

                # Assign curves
                self.curves = {}

            # Now we remove unused points
            self.coor, usedPtsMask = tst.remove_unused_points(
                self.coor, triaConnF=self.triaConnF, quadsConnF=self.quadsConnF
            )

            # Create arrays to store derivative seeds
            self.coord = np.zeros(self.coor.shape)
            self.coorb = np.zeros(self.coor.shape)

            # Create ADT for the current surface, now that coor, triaConn, and quadsConn are correct
            tst.initialize_surface(self)

        else:

            # Other procs create dummy data
            self.coor = np.zeros((1, 3))

            # Create arrays to store derivative seeds
            self.coord = np.zeros(self.coor.shape)
            self.coorb = np.zeros(self.coor.shape)

    def rename(self, name):
        """
        Renames the current surface object.
        """

        # Only the root proc will work here
        if self.myID == 0:

            # Deallocate previous tree
            adtAPI.adtapi.adtdeallocateadts(self.name)

            # Update name
            self.name = name

            # Reinitialize ADT with new name
            tst.update_surface(self)

    def add_curve(self, curve):
        """
        Adds a given curve instance to the self.curves dictionary.
        """

        # Only the root proc will work here
        # if self.myID == 0:

        self.curves[curve.name] = copy.deepcopy(curve)

    def remove_curve(self, name):
        """
        Removes a given curve instance from the self.curves dictionary.
        """

        # Only the root proc will work here
        if self.myID == 0:

            self.curves.pop(name)

    def rename_curve(self, oldName, newName):
        """
        Renames a given curve instance to the self.curves dictionary.
        """

        # Only the root proc will work here
        if self.myID == 0:

            self.curves[newName] = self.curves.pop(oldName)
            self.curves[newName].rename(newName)

    def update(self, coor=None):

        """
        This function updates the nodal coordinates used by the surface object.
        """

        # Only the root proc will work here
        if self.myID == 0:

            if coor is not None:

                # First check if we have the same number of new coordinates
                if not self.coor.shape == coor.shape:
                    print("")
                    print("WARNING: self.update in TSurfGeometry class")
                    print("         The new set of coordinates does not have the")
                    print("         same number of points as the original set.")
                    print("")

                # Update coordinates
                self.coor = coor.copy()

            # Update surface definition
            tst.update_surface(self)

    def translate(self, x, y, z):

        # Only the root proc will work here
        if self.myID == 0:

            tst.translate(self, x, y, z)
            tst.update_surface(self)
            for curve in self.curves.values():
                curve.translate(x, y, z)

    def scale(self, factor):

        # Only the root proc will work here
        if self.myID == 0:

            tst.scale(self, factor)
            tst.update_surface(self)
            for curve in self.curves.values():
                curve.scale(factor)

    def rotate(self, angle, axis, point=None):

        # Only the root proc will work here
        if self.myID == 0:

            tst.rotate(self, angle, axis, point)
            tst.update_surface(self)
            for curve in self.curves.values():
                curve.rotate(angle, axis, point)

    # ===========================================================#
    # SURFACE PROJECTION METHODS

    def project_on_surface(self, xyz):

        """
        This function will compute projections and surface Normals

        Parameters
        ----------
        xyz: float[numPts, 3]
            Coordinates of the points that should be projected.

        Returns
        -------
        xyzProj: float[numPts,3]
            Coordinates of the projected points

        normProj: float[numPts,3]
            Surface normal at projected points

        projDict: dictionary
            Dictionary containing intermediate values that are required by the differentiation routines.
        """

        """
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
        """

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts) * 1e10
        xyzProj = np.zeros((numPts, 3))
        normProjNotNorm = np.zeros((numPts, 3))

        # Call projection function
        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(
            xyz.T, self.name, dist2, xyzProj.T, self.nodal_normals.T, normProjNotNorm.T
        )

        # Adjust indices and ordering
        elementID = elementID - 1
        uvw = uvw.T

        # Store additional outputs in a dictionary to make outputs cleaner
        projDict = {
            "procID": procID,
            "elementType": elementType,
            "elementID": elementID,
            "uvw": uvw,
            "dist2": dist2,
            "normProjNotNorm": normProjNotNorm,
        }

        # Normalize the normals
        normProj = tst.normalize(normProjNotNorm)

        # Return projections
        return xyzProj, normProj, projDict

    def project_on_surface_d(self, xyz, xyzd, xyzProj, normProj, projDict):

        """
        This function will compute derivatives of the projection algorithm in forward mode.

        Parameters
        ----------
        xyz: float[numPts, 3]
            Coordinates of the points that should be projected.

        xyzd: float[numPts, 3]
            Derivative seeds for coordinates of the points
                                   that should be projected.

        xyzProj: float[numPts,3]
            Coordinates of the projected points (can be obtained
                                     with the original projection function)

        normProj: float[numPts,3]
            Surface normal at projected points (can be obtained
                                      with the original projection function)

        projDict: dictionary
            Dictionary containing intermediate values that are
            required by the differentiation routines.  (can be obtained
            with the original projection function)

        Returns
        -------
        xyzProjd: float[numPts,3]
            Derivative seeds of the coordinates of the projected points

        normProjd: float[numPts,3]
            Derivative seeds of the surface normal at projected points


        ATTENTION:

        The user should set coord and curveCoord using the self.set_forwardADSeeds method first.

        coord -> float[numNodes,3] : Derivative seeds of the nodal coordinates of the baseline
                                     surface.

        Ney Secco 2016-10
        """

        # Rename variables to make code more readable
        procID = projDict["procID"]
        elementType = projDict["elementType"]
        elementID = projDict["elementID"]
        uvw = projDict["uvw"]
        dist2 = projDict["dist2"]
        normProjNotNorm = projDict["normProjNotNorm"]

        # Compute derivatives of the normal vectors
        nodal_normals, nodal_normalsd = adtAPI.adtapi.adtcomputenodalnormals_d(
            self.coor.T, self.coord.T, self.triaConnF.T, self.quadsConnF.T
        )
        # Transpose Fortran outputs
        nodal_normals = nodal_normals.T
        nodal_normalsd = nodal_normalsd.T

        # Call projection function
        # ATTENTION: The variable "xyz" here in Python corresponds to the variable "coor" in the Fortran code.
        # On the other hand, the variable "coor" here in Python corresponds to the variable "adtCoor" in Fortran.
        # I could not change this because the original ADT code already used "coor" to denote nodes that should be
        # projected.
        xyzProjd, normProjNotNormd = adtAPI.adtapi.adtmindistancesearch_d(
            xyz.T,
            xyzd.T,
            self.name,
            self.coord.T,
            procID,
            elementType,
            elementID + 1,
            uvw.T,
            dist2,
            xyzProj.T,
            self.nodal_normals.T,
            nodal_normalsd.T,
            normProjNotNorm.T,
        )

        # Transpose results to make them consistent
        xyzProjd = xyzProjd.T
        normProjNotNormd = normProjNotNormd.T

        # Now we need to compute derivatives of the normalization process
        normProj, normProjd = tst.normalize_d(normProjNotNorm, normProjNotNormd)

        # Return projections derivatives
        return xyzProjd, normProjd

    def project_on_surface_b(self, xyz, xyzProj, xyzProjb, normProj, normProjb, projDict):

        """
        This function will compute derivatives of the projection algorithm in reverse mode.

        Parameters
        ----------
        xyz: float[numPts, 3]
            Coordinates of the points that should be projected.

        xyzProj: float[numPts,3]
            Coordinates of the projected points (can be obtained
            with the original projection function)

        xyzProjb: float[numPts,3]
            Derivative seeds of the coordinates of the projected points

        normProj: float[numPts,3]
            Surface normal at projected points (can be obtained
            with the original projection function)

        normProjb: float[numPts,3]
            Derivative seeds of the surface normal at projected points

        projDict: dictionary
            Dictionary containing intermediate values that are
            required by the differentiation routines. (can be obtained
            with the original projection function)

        Returns
        -------

        xyzb: float[numPts, 3]
            Derivative seeds for coordinates of the points that should be projected.

        coorb: float[numNodes,3]
            Derivative seeds of the nodal coordinates of the baseline surface.

        ATTENTION:
        This function just returns the derivate seeds of the points to be projected.
        The derivative seeds of the triangulated surface nodes are acumulated internally.

        Ney Secco 2016-10
        """

        # Rename variables to make code more readable
        procID = projDict["procID"]
        elementType = projDict["elementType"]
        elementID = projDict["elementID"]
        uvw = projDict["uvw"]
        dist2 = projDict["dist2"]
        normProjNotNorm = projDict["normProjNotNorm"]

        # Compute derivatives of the normalization process
        normProjNotNormb = tst.normalize_b(normProjNotNorm, normProjb)

        # Call projection function
        # ATTENTION: The variable "xyz" here in Python corresponds to the variable "coor" in the Fortran code.
        # On the other hand, the variable "coor" here in Python corresponds to the variable "adtCoor" in Fortran.
        # I could not change this because the original ADT code already used "coor" to denote nodes that should be
        # projected.
        xyzb, coorb, nodal_normalsb = adtAPI.adtapi.adtmindistancesearch_b(
            xyz.T,
            self.name,
            procID,
            elementType,
            elementID + 1,
            uvw.T,
            dist2,
            xyzProj.T,
            xyzProjb.T,
            self.nodal_normals.T,
            normProjNotNorm.T,
            normProjNotNormb.T,
        )

        # Transpose results to make them consistent
        xyzb = xyzb.T
        coorb = coorb.T
        nodal_normalsb = nodal_normalsb.T

        # Compute derivative seed contributions of the normal vectors
        deltaCoorb = adtAPI.adtapi.adtcomputenodalnormals_b(
            self.coor.T, self.triaConnF.T, self.quadsConnF.T, self.nodal_normals.T, nodal_normalsb.T
        )

        # Transpose Fortran results to make them consistent
        deltaCoorb = deltaCoorb.T

        # Accumulate normal vector contributions
        coorb = coorb + deltaCoorb

        # Accumulate seeds into the geometry object
        self.accumulate_reverseADSeeds(coorb, curveCoorb=None)

        # Return projection derivatives
        return xyzb

    # ===========================================================#
    # CURVE PROJECTION METHODS

    def project_on_curve(self, xyz, curveCandidates=None):

        """
        This function will compute projections and surface Normals

        Parameters
        ----------
        xyz: float[numPts, 3]
            Coordinates of the points that should be projected.

        Returns
        -------
        xyzProj: float[numPts,3]
            Coordinates of the projected points

        tanProj: float[numPts,3]
            Curve tangent at projected points
        """

        """
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
        """

        # Use all curves if None is provided by the user
        if curveCandidates is None:
            curveCandidates = list(self.curves.keys())

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts) * 1e10
        xyzProj = np.zeros((numPts, 3))
        tanProj = np.zeros((numPts, 3))
        elemIDs = np.zeros((numPts), dtype="int32")

        # Check if the candidates are actually defined
        curveKeys = list(self.curves.keys())
        for curve in curveCandidates:
            if curve not in curveKeys:
                print("ERROR: Curve", curve, "is not defined. Check the curve names in your CGNS file.")
                print("       Also check if you included this curve name when selecting CGNS sections.")
                quit()

        # Initialize list that will contain the names of the curves that got the best projections for
        # each point
        curveIDs = [0] * numPts

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:
            if curveName in curveCandidates:

                # Run projection code
                _, _, _, curveMask = self.curves[curveName].project(xyz, dist2, xyzProj, tanProj, elemIDs)

                # Use curveMask to update names of curves that got best projections
                for ii in range(numPts):
                    if curveMask[ii] == 1:
                        curveIDs[ii] = curveName

        # Create dictionary with useful outputs for differentiated codes
        curveProjDict = {}
        curveProjDict["elemIDs"] = elemIDs
        curveProjDict["curveIDs"] = curveIDs
        curveProjDict["curveCandidates"] = curveCandidates

        # Return projections
        return xyzProj, tanProj, curveProjDict

    def project_on_curve_d(self, xyz, xyzd, xyzProj, tanProj, curveProjDict):

        """
        This function will compute projections and surface Normals

        Parameters
        ----------
        xyz: float[numPts, 3]
            Coordinates of the points that should be projected.

        xyzProj: float[numPts,3]
            Coordinates of the projected points

        tanProj: float[numPts,3]
            Curve tangent at projected points


        ATTENTION:

        The user should set coord and curveCoord using the self.set_forwardADSeeds method first.

        allCoord -> dictionary : This is a dictionary whose keys are all curve names. Each entry
                    should contain an array [nNodes x 3] with derivative seeds for all nodes of the curve

        Ney Secco 2016-11
        """

        # Retrieve data from dictionary
        elemIDs = curveProjDict["elemIDs"]
        curveIDs = curveProjDict["curveIDs"]
        curveCandidates = curveProjDict["curveCandidates"]

        # Get number of points
        numPts = xyz.shape[0]

        # Initialize derivatives
        xyzProjd = np.zeros((numPts, 3))
        tanProjd = np.zeros((numPts, 3))

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:
            if curveName in curveCandidates:

                # Identify points that are projected to the current curve to create a mask
                curveMask = [0] * numPts
                for ii in range(numPts):
                    if curveIDs[ii] == curveName:
                        curveMask[ii] = 1

                # Run projection code in forward mode
                self.curves[curveName].project_d(xyz, xyzd, xyzProj, xyzProjd, tanProj, tanProjd, elemIDs, curveMask)

        # Return projections derivatives
        return xyzProjd, tanProjd

    def project_on_curve_b(self, xyz, xyzProj, xyzProjb, tanProj, tanProjb, curveProjDict):

        """
        This function will backpropagate projections derivatives back
        to curve nodes and unprojected points.

        Ney Secco 2016-11
        """

        # Retrieve data from dictionary
        elemIDs = curveProjDict["elemIDs"]
        curveIDs = curveProjDict["curveIDs"]
        curveCandidates = curveProjDict["curveCandidates"]

        # Get number of points
        numPts = xyz.shape[0]

        # Initialize derivatives
        xyzb = np.zeros((numPts, 3))

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.curves:

            if curveName in curveCandidates:

                # Identify points that are projected to the current curve to create a mask
                curveMask = [0] * numPts
                for ii in range(numPts):
                    if curveIDs[ii] == curveName:
                        curveMask[ii] = 1

                # Run projection code in reverse mode
                # This will modify xyzb and coorb
                self.curves[curveName].project_b(xyz, xyzb, xyzProj, xyzProjb, tanProj, tanProjb, elemIDs, curveMask)

        # Return backpropagated derivatives
        return xyzb

    def extract_curves(self, feature="sharpness"):
        """
        This function will define new curves in this component based on
        surface features. The new curves will be stored in self.curves
        """
        tst.extract_curves_from_surface(self, feature)

    # ===========================================================#
    # INTERSECTION METHODS

    def intersect(self, otherGeometry, distTol=1e-7):
        """
        This method will intersect the current component (self) with the provided
        TSurfGeometry. The new curves will be stored in each component.

        ATTENTION:
        For every geomA.intersect(geomB) in a forward code, there should be a
        corresponding geomA.intersect_d(geomB) or geomA.intersect_b(geomB) in the reverse code.

        Parameters
        ----------
        otherGeometry: Geometry object
            Other object that we want to intersect with.

        distTol: float
            Tolerance used to merge close nodes when joining bar elements generated by the intersection.

        Returns
        -------

        Intersection: list of curve objects generated by the intersection operation.

        """

        # Call the intersection function defined in tsurf_tools.py
        Intersection = tst._compute_pair_intersection(self, otherGeometry, distTol)

        # Print the number of intersections
        print("Number of intersections between", self.name, "and", otherGeometry.name, "is", len(Intersection))

        # Return the intersection curves
        return Intersection

    def intersect_d(self, otherGeometry, intCurve, distTol=1e-7):
        """
        This method will propagate forward mode AD derivatives for all intersection curves
        between the current component and the otherGeometry component.
        Remember to use the same distTol of the forward pass.

        Ney Secco 2017-01
        """

        # Check if parents are correct

        if (self.name == intCurve.extra_data["parentGeoms"][0]) and (
            otherGeometry.name == intCurve.extra_data["parentGeoms"][1]
        ):

            # Print log
            print("")
            print("Propagating forward AD seeds for intersection curve:")
            print(intCurve.name)
            print("")

        else:

            raise NameError("Trying to run derivative code with incorrect parents")

        # Then we can call the AD code defined in tsurf_tools.py
        coorIntd = tst._compute_pair_intersection_d(
            self, otherGeometry, intCurve, self.coord, otherGeometry.coord, distTol
        )

        # Place the derivative seeds into the curve object
        intCurve._set_forwardADSeeds(coorIntd)

    def intersect_b(self, otherGeometry, intCurve, distTol=1e-7, accumulateSeeds=True):
        """
        This method will propagate reverse mode AD derivatives for all intersection curves
        between the current component and the otherGeometry component.
        Remember to use the same distTol of the forward pass.

        Ney Secco 2017-01
        """

        # Check if parents are correct

        if (self.name == intCurve.extra_data["parentGeoms"][0]) and (
            otherGeometry.name == intCurve.extra_data["parentGeoms"][1]
        ):

            # Print log
            print("")
            print("Propagating reverse AD seeds for intersection curve:")
            print(intCurve.name)
            print("")

        else:

            raise NameError("Trying to run derivative code with incorrect parents")

        # Copy the intersection curve seeds to avoid any modification
        intCoorb = np.array(intCurve.coorb)

        # Then we can call the AD code defined in tsurf_tools.py
        coorAb, coorBb = tst._compute_pair_intersection_b(self, otherGeometry, intCurve, intCoorb, distTol)

        # Accumulate or set derivative seeds into the curve object on both sides
        if accumulateSeeds:
            self.accumulate_reverseADSeeds(coorb=coorAb)
            otherGeometry.accumulate_reverseADSeeds(coorb=coorBb)
        else:
            self.set_reverseADSeeds(coorb=coorAb)
            otherGeometry.set_reverseADSeeds(coorb=coorBb)

    # ===========================================================#
    # POINT METHODS

    def set_points(self, coor):

        """
        This will replace the nodal coordinates that define the surface.

        Parameters
        ----------
        coor: float[nPts,3]
            Nodal coordinates (X,Y,Z)
        """

        self.update(np.array(coor))

    def get_points(self):

        """
        This will give the nodal coordinates that define the surface.
        """

        return self.coor

    # ===========================================================#
    # DERIVATIVE SEED MANIPULATION METHODS

    def set_forwardADSeeds(self, coord=None, curveCoord=None):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        Parameters
        ----------
        coord: float[numSurfNodes, 3]
            Derivative seeds of the nodal coordinates of the reference
            surface. It should have the same shape as self.coor.

        curveCoord : [curveName]{float[numCurveNodes, 3]}
            Dictionary containing derivative seeds
            of the nodal coordinates of each curve
            present in the current object. The dictionary keys
            are the names of the curves.
        """

        if coord is not None:
            self.coord[:, :] = np.array(coord)

        if curveCoord is not None:
            for curveName in self.curves:
                if curveName in list(curveCoord.keys()):
                    if curveCoord[curveName] is not None:
                        self.curves[curveName].set_forwardADSeeds(curveCoord[curveName])

    def get_forwardADSeeds(self):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        coord = np.array(self.coord)

        curveCoord = {}
        for curveName in self.curves:
            curveCoord[curveName] = self.curves[curveName].get_forwardADSeeds()

        return coord, curveCoord

    def set_reverseADSeeds(self, coorb=None, curveCoorb=None):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        if coorb is not None:
            self.coorb = np.array(coorb)

        if curveCoorb is not None:
            for curveName in self.curves:
                if curveName in list(curveCoorb.keys()):
                    if curveCoorb[curveName] is not None:
                        self.curves[curveName].set_reverseADSeeds(curveCoorb[curveName])

    def get_reverseADSeeds(self, clean=True):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # We use np.arrays to make full copies
        coorb = np.array(self.coorb)

        # Check if we need to clean the derivative seeds
        if clean:
            self.coorb[:, :] = 0.0

        # Get curve derivatives
        curveCoorb = {}
        for curveName in self.curves:
            curveCoorb[curveName] = self.curves[curveName].get_reverseADSeeds(clean)

        # We have no mesh information to return
        return coorb, curveCoorb

    def accumulate_reverseADSeeds(self, coorb=None, curveCoorb=None):
        """
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoorb: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        if coorb is not None:
            self.coorb = self.coorb + np.array(coorb)

        if curveCoorb is not None:
            for curveName in self.curves:
                if curveCoorb[curveName] is not None:
                    self.curves[curveName].coorb = self.accumulate_reverseADSeeds(curveCoorb[curveName])

    def clean_reverseADSeeds(self):
        """
        This will erase all reverse derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        self.coorb[:, :] = 0.0

        for curveName in self.curves:
            self.curves[curveName].coorb[:, :] = 0.0

    def set_randomADSeeds(self, mode="both", fixedSeed=True):

        """
        This will set random normalized seeds to all variables.
        This can be used for testing purposes.

        mode: ['both','forward','reverse'] -> Which mode should have
        its derivatives replaced.
        """

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # Set forward AD seeds
        if mode == "forward" or mode == "both":

            coord = np.random.random_sample(self.coor.shape)
            coord = coord / np.sqrt(np.sum(coord ** 2))
            self.coord = coord

            curveCoord = {}
            for curveName in self.curves:
                self.curves[curveName].set_randomADSeeds(mode="forward", fixedSeed=fixedSeed)
                curveCoord[curveName] = self.curves[curveName].get_forwardADSeeds()

        # Set reverse AD seeds
        if mode == "reverse" or mode == "both":

            coorb = np.random.random_sample(self.coor.shape)
            coorb = coorb / np.sqrt(np.sum(coorb ** 2))
            self.coorb = coorb

            curveCoorb = {}
            for curveName in self.curves:
                self.curves[curveName].set_randomADSeeds(mode="reverse", fixedSeed=fixedSeed)
                curveCoorb[curveName] = self.curves[curveName].get_reverseADSeeds(clean=False)

        # Return the randomly generated seeds

        if mode == "forward":
            return coord, curveCoord

        elif mode == "reverse":
            return coorb, curveCoorb

        elif mode == "both":
            return coord, curveCoord, coorb, curveCoorb

    # ===========================================================#
    # VISUALIZATION METHODS

    def export_tecplot(self, fileName):

        """
        This method will export the triangulated surface data into a tecplot format.
        This will write all elements as quad data. The triangle elements will be exported
        as degenerated quads (quads with a repeated node).

        Ney Secco 2017-02
        """

        # Only the root proc will work here
        if self.myID == 0:

            # Call the function that export FE data
            # we add -1 here to convert from Fortran to Python ordering
            # The Tecplot function shouldn't assume we have Fortran ordered data.
            tecplot_interface.writeTecplotSurfaceFEData(
                self.coor, self.triaConnF - 1, self.quadsConnF - 1, self.name, fileName
            )


# =============================================================
# =============================================================
# =============================================================


class TSurfCurve(Curve):

    """
    This is a derived class from the parent Curve class defined in
    baseClasses.py
    """

    def _initialize(self, *arg):  # ,coor, barsConn, name, mergeTol=1e-2):

        """
        This method initializes a TSurf Curve object, which uses finite element
        data to represent curves.
        This method will be called by the __init__ method of the parent Curve class

        Parameters
        ----------
        name: string
            Curve name

        coor : array[nNodes,3],dtype='float'
            Nodal X,Y,Z coordinates.

        barsConn : array[nBars,2],dtype='int32'
            Element connectivity matrix.

        mergeTol: float, optional
            Tolerance to merge nodes.

        John Jasa 2016-08
        Ney Secco 2016-08
        """

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
                dim = currArg.shape[1]

                if dim == 3:
                    # We have nodal coordinates. We need to make sure they
                    # are defined as floats, otherwise Fortran will not recognize them.
                    coor = np.array(currArg, dtype=float)

                elif dim == 2:
                    # We have bar connectivities. We need to make sure they
                    # are defined as integers, otherwise Fortran will not recognize them.
                    barsConn = np.array(currArg, dtype=np.int32)

                else:
                    print(" ERROR: Could not recognize array inputs when initializing TSurfCurve.")
                    print(" Please provide an [n x 3] array with nodal coordinates and an")
                    print(" [m x 2] array for bar connectivities.")
                    quit()

            elif type(currArg) == float:
                # This is probably mergeTol
                mergeTol = currArg

            else:
                print(" ERROR: Could not recognize inputs when initializing TSurfCurve.")
                quit()

        # Remove unused points (This will update and barsConn)
        coor, usedPtsMask = tst.remove_unused_points(coor, barsConn=barsConn)

        # Call a Fortran code to merge close nodes (This will update coor and barsConn)
        barsConnF = (
            barsConn.T + 1
        )  # Adjust indices since Fortran is 1-based (we need to do this separately as Fortran changes it)
        nUniqueNodes, linkOld2New = utilitiesAPI.utilitiesapi.condensebarnodes(mergeTol, coor.T, barsConnF)

        # Readjust back to Python indexing
        linkOld2New = linkOld2New - 1
        barsConn = barsConnF.T - 1

        # Save the mapping
        self.extra_data = {}
        self.extra_data["linkOld2New"] = linkOld2New
        self.extra_data["usedPtsMask"] = usedPtsMask

        # Initialize other fields that will be used by the projection code
        self.extra_data["parentGeoms"] = []
        self.extra_data["parentTria"] = []
        self.extra_data["childMeshes"] = []

        # Take bar connectivity and reorder it
        sortedConn, dummy_map = tst.FEsort(barsConn.tolist())

        # Check for failures
        if len(sortedConn) == 1:
            # The sorting works just fine and it found a single curve!
            # So we just store this single curve (sortedConn[0]).
            sortedConn = np.array(sortedConn[0], dtype=type(barsConn[0, 0]))
        else:
            # Just use the original connectivity
            sortedConn = barsConn
            # Print message
            print("")
            print("Curve", '"' + name + '"', "could not be sorted. It might be composed by disconnected curves.")
            print("")

        # Assing coor and barsConn. Remember to crop coor to get unique nodes only
        self.coor = coor[:nUniqueNodes, :]
        self.barsConn = np.array(sortedConn, dtype="int32")
        self.name = name
        self.numNodes = self.coor.shape[0]

        # Create arrays to store derivatives of the nodes
        self.coord = np.zeros(self.coor.shape)
        self.coorb = np.zeros(self.coor.shape)

        # Check if curve is periodic by seeing if the initial and final points
        # are the same
        if self.barsConn[0, 0] == self.barsConn[-1, -1]:
            self.isPeriodic = True
        else:
            self.isPeriodic = False

    def update(self, coor):
        self.coor = np.array(coor)

    def flip(self):

        self.barsConn = self.barsConn[::-1, ::-1]

        # Flip the extra data associated with element ordering
        for data in self.extra_data:
            if data in ["parentTria"]:
                self.extra_data[data] = self.extra_data[data][::-1]

    def translate(self, x, y, z):
        tst.translate(self, x, y, z)

    def scale(self, factor):
        tst.scale(self, factor)

    def rotate(self, angle, axis, point=None):
        tst.rotate(self, angle, axis, point)

    def rename(self, name):
        self.name = name

    # =================================================================
    # PROJECTION FUNCTIONS
    # =================================================================

    def project(self, xyz, dist2=None, xyzProj=None, tanProj=None, elemIDs=None):

        """
        This function will take the points given in xyz and project them to the curve.

        Parameters
        ----------
        xyz: float[nPoints,3]
            coordinates of the points that should be projected

        dist2: float[nPoints]
            values of best distance**2 found so far. If the distance of the projected
            point is less than dist2, then we will take the projected point. This allows us to find the best
            projection even with multiple surfaces.
            If you don't have previous values to dist2, just initialize all elements to
            a huge number (1e10).

        xyzProj: float[nPoints,3]
            coordinates of projected points found so far. These projections could be on other
            curves as well. This function will only replace projections whose dist2 are smaller
            than previous dist2. This allows us to use the same array while working with multiple curves.
            If you don't have previous values, just initialize all elements to zero. Also
            remember to set dist2 to a huge number so that all values are replaced.

        tanProj: float[nPoints,3]
            tangent directions for the curve at the projected points.

        elemIDs: int[nPoints]
            ID of the bar elements that received projections. This is also given by
            the execution of the primal routine.

        This function has no explicit outputs. It will just update dist2, xyzProj, tangents, and elemIDs
        """

        # Get number of points
        nPoints = xyz.shape[0]

        # Initialize references if user provided none
        if dist2 is None:
            dist2 = np.ones(nPoints) * 1e10
        if xyzProj is None:
            xyzProj = np.zeros((nPoints, 3))
        if tanProj is None:
            tanProj = np.zeros((nPoints, 3))
        if elemIDs is None:
            elemIDs = np.zeros((nPoints), dtype="int32")

        # Call fortran code
        # This will modify xyzProj, tanProj, dist2, and elemIDs if we find better projections than dist2.
        # Remember that we should adjust some indices before calling the Fortran code
        # Remember to use [:] to don't lose the pointer (elemIDs is an input/output variable)
        elemIDs[:] = (
            elemIDs + 1
        )  # (we need to do this separetely because Fortran will actively change elemIDs contents.
        curveMask = curveSearchAPI.curvesearchapi.mindistancecurve(
            xyz.T, self.coor.T, self.barsConn.T + 1, xyzProj.T, tanProj.T, dist2, elemIDs
        )

        # Adjust indices back to Python standards
        elemIDs[:] = elemIDs - 1

        return xyzProj, tanProj, elemIDs, curveMask

    def project_d(self, xyz, xyzd, xyzProj, xyzProjd, tanProj, tanProjd, elemIDs, curveMask):

        """
        This function will run forward mode AD to propagate derivatives from inputs (xyz and coor) to outputs (xyzProj).

        Parameters
        ----------
        xyz: float[nPoints,3]
            coordinates of the points that should be projected.

        xyzd: float[nPoints,3]
            derivatives seeds of the coordinates.

        coord: float[nNodes, 3]
            derivative seeds of the nodes that constitutes the bar elements.

        xyzProj: float[nPoints,3]
            coordinates of projected points found with the primal routine.

        xyzProjd: float[nPoints,3]
            Derivative seeds of the projected points.

        elemIDs: int[nPoints]
            ID of the bar elements that received projections. This is also given by
            the execution of the primal routine.

        curveMask: int[nPoints]
            curveMask[ii] should be 1 if the ii-th point was actually projected onto
            this curve. Otherwise, curveMaks[ii] = 0, then the code will not compute
            derivatives for this point.

        Returns
        -------
        This function has no explicit outputs. It will just update xyzProjd.

        Ney Secco 2016-11
        """

        if self.coor.shape != self.coord.shape:
            print("ERROR: Derivative seeds should have the same dimension of the original")
            print("variable. The number of derivatives for the bar element nodes does not match")
            print("with the number of nodes.")

        # Call fortran code
        curveSearchAPI.curvesearchapi.mindistancecurve_d(
            xyz.T,
            xyzd.T,
            self.coor.T,
            self.coord.T,
            self.barsConn.T + 1,
            xyzProj.T,
            xyzProjd.T,
            tanProj.T,
            tanProjd.T,
            elemIDs + 1,
            curveMask,
        )

    def project_b(self, xyz, xyzb, xyzProj, xyzProjb, tanProj, tanProjb, elemIDs, curveMask):

        """
        This function will run forward mode AD to propagate derivatives from inputs (xyz and coor) to outputs (xyzProj).

        Parameters
        ----------
        xyz: float[nPoints,3]
            coordinates of the points that should be projected.

        xyzb: float[nPoints,3]
            derivatives seeds of the coordinates.

        coorb: float[nNodes, 3]
            derivative seeds of the nodes that constitutes the bar elements.

        xyzProj: float[nPoints,3]
            coordinates of projected points found with the primal routine.

        xyzProjb: float[nPoints,3]
            Derivative seeds of the projected points.

        elemIDs: int[nPoints]
            ID of the bar elements that received projections. This is also given by
            the execution of the primal routine.

        curveMask: int[nPoints]
            curveMask[ii] should be 1 if the ii-th point was actually projected onto
            this curve. Otherwise, curveMaks[ii] = 0, then the code will not compute
            derivatives for this point.

        Returns
        -------
        This function has no explicit outputs. It will just update xyzb and coorb.

        Ney Secco 2016-11
        """

        if xyzProj.shape != xyzProjb.shape:
            print("ERROR: Derivative seeds should have the same dimension of the original")
            print("variable. The number of derivatives for the projected points does not match")
            print("with the number of points.")

        # Call fortran code (This will accumulate seeds in xyzb and self.coorb)
        xyzb_new, coorb_new = curveSearchAPI.curvesearchapi.mindistancecurve_b(
            xyz.T,
            self.coor.T,
            self.barsConn.T + 1,
            xyzProj.T,
            xyzProjb.T,
            tanProj.T,
            tanProjb.T,
            elemIDs + 1,
            curveMask,
        )

        # Accumulate derivatives
        xyzb[:, :] = xyzb + xyzb_new.T
        self.coorb = self.coorb + coorb_new.T

    # =================================================================
    # REMESHING FUNCTIONS
    # =================================================================

    def remesh(self, nNewNodes=None, method="linear", spacing="linear", initialSpacing=0.1, finalSpacing=0.1):

        """
        This function will redistribute the nodes along the curve to get
        better spacing among the nodes.
        This assumes that the FE data is ordered.

        The first and last nodes will remain at the same place.

        Consider using self.shift_end_nodes first so you can preserve the start and end points
        of a periodic curve.

        Parameters
        ----------

        nNewNodes: integer
            Number of new nodes desired in the new curve definition.

        method: string
            Method used to interpolate new nodes with respect to
            existing ones. Check scipy.interpolate.interp1d for options.

        spacing: string
            Desired spacing criteria for new nodes.
            Current options are:['linear', 'cosine', 'hypTan', 'tangent']

        initialSpacing: float
            Desired distance between the first and second nodes. Only
            used by 'hypTan' and 'tangent'.
        finalSpacing: float
            Desired distance between the last two nodes. Only
            used by 'hypTan' and 'tangent'.
        guideCurves: list of strings
            Curves to snap nodes to.
            Especially useful for blunt trailing edges. (removed since it was not differentiated)
        ref_geom: geometry object
            Container with the data for each curve used in guideCurves. (removed since it was not differentiated)

        Returns
        -------
        newCurve: curve object
            This is the remeshed curve object.

        Ney Secco 2016-08
        """

        spacing = spacing.lower()

        # Make copies of the connectivities and coordinates of the current Curve object
        coor = np.array(self.coor)
        barsConn = np.array(self.barsConn)

        # Get the number of elements in the curve
        nElem = barsConn.shape[0]
        nNodes = nElem + 1

        # Check if the baseline curve is periodic. If this is the case, we artificially repeat
        # the last point so that we could use the same code of the non-periodic case
        if barsConn[0, 0] == barsConn[-1, 1]:
            periodic = True
            coor = np.vstack([coor, coor[barsConn[0, 0], :].reshape((1, 3))])
            barsConn[-1, -1] = nNodes - 1
        else:
            periodic = False

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        if fortran_flag:

            # Call Fortran code. Remember to adjust transposes and indices
            newCoor, newBarsConn = utilitiesAPI.utilitiesapi.remesh(
                nNewNodes, coor.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
            )
            newCoor = newCoor.T
            newBarsConn = newBarsConn.T - 1

        else:

            # CHECKING INPUTS

            # First we check if the FE data is ordered
            for elemID in range(1, nElem):

                # Get node indices
                prevNodeID = barsConn[elemID - 1, 1]
                currNodeID = barsConn[elemID, 0]

                # Check if the FE data is ordered
                if prevNodeID != currNodeID:

                    # Print warning
                    print("WARNING: Could not remesh curve because it has unordered FE data.")
                    print("         Call FEsort first.")
                    return

            # COMPUTE ARC-LENGTH

            # We can proceed if FE data is ordered

            # Initialize an array that will store the arclength of each node.
            # That is, the distance, along the curve from the current node to
            # the first node of the curve.
            arcLength = np.zeros(nNodes)

            # Also initialize an array to store nodal coordinates in the curve order
            nodeCoor = np.zeros((nNodes, 3))

            # Store position of the first node (the other nodes will be covered in the loop)
            # (the -1 is due Fortran indexing)
            nodeCoor[0, :] = coor[barsConn[0, 0], :]

            # Loop over each element to increment arcLength
            for elemID in range(nElem):

                # Get node positions
                node1 = coor[barsConn[elemID, 0], :]
                node2 = coor[barsConn[elemID, 1], :]

                # Compute distance between nodes
                dist = np.linalg.norm(node1 - node2)

                # Store nodal arc-length
                arcLength[elemID + 1] = arcLength[elemID] + dist

                # Store coordinates of the next node
                nodeCoor[elemID + 1, :] = node2

            # SAMPLING POSITION FOR NEW NODES

            # The arcLength will be used as our parametric coordinate to sample the new nodes.

            # First we will initialize arrays for the new coordinates and arclengths
            # newNodeCoor = np.zeros((3, nNewNodes), order='F') #3change here
            newNodeCoor = np.zeros((nNewNodes, 3))
            newArcLength = np.zeros(nNewNodes)

            # Now that we know the initial and final arcLength, we can redistribute the
            # parametric coordinates based on the used defined spacing criteria.
            # These statements should initially create parametric coordinates in the interval
            # [0.0, 1.0]. We will rescale it after the if statements.
            if spacing.lower() == "linear":
                newArcLength = np.linspace(0.0, 1.0, nNewNodes)

            elif spacing.lower() == "cosine":
                newArcLength = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, nNewNodes)))

            elif spacing.lower() == "hyptan":
                newArcLength = tst.hypTanDist(initialSpacing / arcLength[-1], finalSpacing / arcLength[-1], nNewNodes)

            elif spacing.lower() == "tangent":
                newArcLength = tst.tanDist(initialSpacing / arcLength[-1], finalSpacing / arcLength[-1], nNewNodes)

            # Rescale newArcLength based on the final distance
            newArcLength = arcLength[-1] * newArcLength

            # INTERPOLATE NEW NODES

            # Now we sample the new coordinates based on the interpolation method given by the user

            # Import interpolation function
            from scipy.interpolate import interp1d

            # Create interpolants for x, y, and z
            fX = interp1d(arcLength, nodeCoor[:, 0], kind=method)
            fY = interp1d(arcLength, nodeCoor[:, 1], kind=method)
            fZ = interp1d(arcLength, nodeCoor[:, 2], kind=method)

            # Sample new points using the interpolation functions
            newNodeCoor[:, 0] = fX(newArcLength)
            newNodeCoor[:, 1] = fY(newArcLength)
            newNodeCoor[:, 2] = fZ(newArcLength)

            """
            # The next commands are to snap a point to guide curves. Since this was not differentiated
            # and also not added to the Fortran code, I'll leave this commented until future use.

            # If you want to use this, you should include guideCurves=[] and ref_geom=None
            # as arguments to this function.

            guideIndices = []
            for curve in guideCurves:
                guideCurve = ref_geom.curves[curve].extract_points()
                guideIndices.append(closest_node(guideCurve, newNodeCoor.T))

            if guideIndices:
                for i, index in enumerate(guideIndices):
                    curve = guideCurves[i]
                    node = newNodeCoor[:, index]
                    newNodeCoor[:, index], _, __ = ref_geom.project_on_curve((node.reshape(1, 3)), curveCandidates=[curve])
            """

            # ASSIGN NEW COORDINATES AND CONNECTIVITIES

            # Set new nodes
            newCoor = newNodeCoor

            # Generate new connectivity (the nodes are already in order so we just
            # need to assign an ordered set to barsConn).
            barsConn = np.zeros((nNewNodes - 1, 2), dtype="int32")
            barsConn[:, 0] = list(range(0, nNewNodes - 1))
            barsConn[:, 1] = list(range(1, nNewNodes))

            newBarsConn = barsConn

        # Adjust connectivities if curve is periodic
        if periodic:
            newBarsConn[-1, -1] = newBarsConn[0, 0]
            newCoor = newCoor[:-1, :]

        # Generate name for the new curve
        newCurveName = self.name + "_remeshed"

        # Create a new curve object and return it
        # This way the original curve coordinates and connectivities remain the same
        newCurve = TSurfCurve(newCurveName, newCoor, newBarsConn)

        return newCurve

    # =================================================================

    def remesh_d(
        self, newCurve, nNewNodes=None, method="linear", spacing="linear", initialSpacing=0.1, finalSpacing=0.1
    ):

        """
        The original curve (before the remesh) should be the one that calls this method.
        """

        spacing = spacing.lower()

        # Get derivative seeds
        coord = self._get_forwardADSeeds()

        # Get connectivities and coordinates of the current Curve object
        coor = np.array(self.coor)
        barsConn = np.array(self.barsConn)

        # Get the number of elements in the curve
        nElem = barsConn.shape[0]
        nNodes = nElem + 1

        # Check if the baseline curve is periodic. If this is the case, we artificially repeat
        # the last point so that we could use the same code of the non-periodic case
        if barsConn[0, 0] == barsConn[-1, 1]:
            periodic = True
            coor = np.vstack([coor, coor[barsConn[0, 0], :].reshape((1, 3))])
            coord = np.vstack([coord, coord[barsConn[0, 0], :].reshape((1, 3))])
            barsConn[-1, -1] = nNodes - 1
        else:
            periodic = False

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        nNewElems = nNewNodes - 1

        # Call Fortran code. Remember to adjust transposes and indices
        newCoor, newCoord, newBarsConn = utilitiesAPI.utilitiesapi.remesh_d(
            nNewNodes, nNewElems, coor.T, coord.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
        )

        # Transpose results coming out of Fortran
        newCoor = newCoor.T
        newCoord = newCoord.T
        newBarsConn = newBarsConn.T - 1

        # DO A LOCAL FD TEST
        stepSize = 1e-7
        coord = coord / np.sqrt(np.sum(coor ** 2))
        coor = coor + coord * stepSize

        # Call Fortran code. Remember to adjust transposes and indices
        newCoorf, newBarsConn = utilitiesAPI.utilitiesapi.remesh(
            nNewNodes, coor.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
        )

        # Transpose results coming out of Fortran
        newCoorf = newCoorf.T

        newCoord_FD = (newCoorf - newCoor) / stepSize
        print("FD test @ remesh_d")
        print(np.max(np.abs(newCoord_FD - newCoord)))

        # Adjust seeds if curve is periodic
        if periodic:
            newCoord[0, :] = 0.5 * (newCoord[0, :] + newCoord[-1, :])
            newCoord = newCoord[:-1, :]

        # Store propagated seeds
        newCurve._set_forwardADSeeds(newCoord)

    # =================================================================

    def remesh_b(
        self,
        newCurve,
        clean=True,
        accumulateSeeds=True,
        nNewNodes=None,
        method="linear",
        spacing="linear",
        initialSpacing=0.1,
        finalSpacing=0.1,
    ):

        """
        The original curve (before the remesh) should be the one that calls this method.
        """

        spacing = spacing.lower()

        # Get derivative seeds from the new curve
        newCoorb = newCurve._get_reverseADSeeds(clean=clean)

        # Get connectivities and coordinates of the current Curve object
        coor = np.array(self.coor)
        barsConn = np.array(self.barsConn)

        # Get the number of elements in the curve
        nElem = barsConn.shape[0]
        nNodes = nElem + 1

        # Check if the baseline curve is periodic. If this is the case, we artificially repeat
        # the last point so that we could use the same code of the non-periodic case
        if barsConn[0, 0] == barsConn[-1, 1]:
            periodic = True
            coor = np.vstack([coor, coor[barsConn[0, 0], :].reshape((1, 3))])
            newCoorb = np.vstack([newCoorb, newCoorb[0, :].reshape((1, 3))])
            newCoorb[0, :] = 0.5 * newCoorb[0, :]
            newCoorb[-1, :] = 0.5 * newCoorb[-1, :]
            barsConn[-1, -1] = nNodes - 1

        else:
            periodic = False

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        nNewElems = nNewNodes - 1

        # Call Fortran code. Remember to adjust transposes and indices
        _, __, coorb = utilitiesAPI.utilitiesapi.remesh_b(
            nNewElems, coor.T, newCoorb.T, barsConn.T + 1, method, spacing, initialSpacing, finalSpacing
        )

        coorb = coorb.T

        # Adjust seeds if curve is periodic
        if periodic:
            coorb[barsConn[0, 0], :] += coorb[-1, :]
            coorb = coorb[:-1, :]

        # Store propagated seeds
        if accumulateSeeds:
            self._accumulate_reverseADSeeds(coorb)
        else:
            self._set_reverseADSeeds(coorb)

    # =================================================================
    # SPLITTING METHODS

    def split(self, optionsDict={}, criteria="sharpness"):
        """
        This method will split the current curve based on a given criteria.

        ATTENTION: This function assumes that the FE data is sorted.

        Parameters
        ----------

        criteria: string -> Criteria that will be used to split curves. The options
        available for now are: ['sharpness', 'curve']

        Returns
        -------
        splitCurveDict: dictionary[curve objects] -> Dictionary containing split curves.

        Ney Secco 2017-01
        """

        # Call splitting function defined in tsurf_tools.py
        splitCurvesDict = tst.split_curve_single(self, self.name, optionsDict, criteria)

        # Return dictionary of new curves
        return splitCurvesDict

    def split_d(self, childCurve):
        """
        This method will propagate the forward AD Seeds from the parent curve
        to its child (which was generated from the splitting process).

        Ney Secco 2017-01
        """

        # Check if the given curve is actually a child
        if childCurve.name in list(self.extra_data["splitCurves"].keys()):

            # The usedPtsMask of the child curve will help us inherit derivative seeds from
            # the parent curve

            usedPtsMask = childCurve.extra_data["usedPtsMask"]

            for parentNodeID in range(len(usedPtsMask)):

                # Get child node that corresponds to the parent node
                childNodeID = usedPtsMask[parentNodeID]

                # Update derivatives only if the parent node is present in the child curve
                if childNodeID >= 0:
                    childCurve.coord[childNodeID, :] = self.coord[parentNodeID, :]

    def split_b(self, childCurve):
        """
        This method will propagate the forward AD Seeds from the child curve
        to its parent (which was generated from the splitting process).

        Ney Secco 2017-01
        """

        # Check if the given curve is actually a child
        if childCurve.name in list(self.extra_data["splitCurves"].keys()):

            # The usedPtsMask of the child curve will help us inherit derivative seeds from
            # the parent curve

            usedPtsMask = childCurve.extra_data["usedPtsMask"]

            for parentNodeID in range(len(usedPtsMask)):

                # Get child node that corresponds to the parent node
                childNodeID = usedPtsMask[parentNodeID]

                # Update derivatives only if the parent node is present in the child curve
                if childNodeID >= 0:
                    self.coorb[parentNodeID, :] = self.coorb[parentNodeID, :] + childCurve.coorb[childNodeID, :]

    # =================================================================
    # MERGING METHODS

    def merge(self, curveDict, mergedCurveName, curvesToMerge=None):
        """
        This method will merge the current curve with all other curves in
        the curveDict whose name is in curvesToMerge.
        """

        # Add the current curve to the dictionary
        if self.name not in list(curveDict.keys()):
            curveDict[self.name] = self

        # Merge all curves if user provided none
        if curvesToMerge is None:
            curvesToMerge = list(curveDict.keys())

        # Check if the current curve is in the curvesToMerge list
        if self.name not in curvesToMerge:
            curvesToMerge.append(self.name)

        # Call the merge_curves function defined in tsurf_tools.py
        mergedCurve = tst.merge_curves(curveDict, mergedCurveName, curvesToMerge)

        # This mergeCurve has mergedCurve.extra_data['linkOld2New'], which will help
        # propagate derivatives.

        # Add list of parents to the merged curve
        mergedCurve.extra_data["parentCurves"] = curvesToMerge

        # Return the new curve
        return mergedCurve

    def merge_d(self, curveDict):
        """
        This method will propagate forward AD derivatives from the parent curves to
        the merged curve.

        ATTENTION: This method should be called from the MERGED CURVE since it knows
        its parents.
        """

        # Initialize index offset. The new curves cannot reference
        # nodes of the old ones.
        indexOffset = 0

        # Get mapping between parent and child nodes
        linkOld2New = self.extra_data["linkOld2New"]

        # Initialize list that states how many nodes where merged to generate a single node
        # of the merged curve.
        numMergedNodes = np.zeros(self.coor.shape[0], dtype="int")

        # Loop over every curve we want to merge
        for curveName in self.extra_data["parentCurves"]:

            # Check if we have access to this curve
            if curveName in curveDict:

                # Get derivative seeds of the parent curve
                coord = curveDict[curveName]._get_forwardADSeeds()

                # Loop over the parent curve nodes to propagate seeds
                for parentNodeID in range(coord.shape[0]):

                    # Get corresponding ID in the child nodes
                    childNodeID = linkOld2New[indexOffset + parentNodeID]

                    # Acumulate derivative seed into the child node
                    # (we will take the average later on)
                    self.coord[childNodeID, :] = (
                        self.coord[childNodeID, :] + curveDict[curveName].coord[parentNodeID, :]
                    )

                    # Increment number of parent nodes added to this child node
                    numMergedNodes[childNodeID] = numMergedNodes[childNodeID] + 1

                # Update offset for next curve
                indexOffset = indexOffset + curveDict[curveName].coor.shape[0]

        # Compute averages
        for childNodeID in range(self.coor.shape[0]):
            self.coord[childNodeID, :] = self.coord[childNodeID, :] / numMergedNodes[childNodeID]

    def merge_b(self, curveDict):
        """
        This method will propagate reverse AD derivatives from the merged curve to
        the parent curves.

        ATTENTION: This method should be called from the MERGED CURVE since it knows
        its parents.
        """

        # Initialize index offset. The new curves cannot reference
        # nodes of the old ones.
        indexOffset = 0

        # Get mapping between parent and child nodes
        linkOld2New = self.extra_data["linkOld2New"]

        # Initialize list that states how many nodes where merged to generate a single node
        # of the merged curve.
        numMergedNodes = np.zeros(self.coor.shape[0], dtype="int")

        # We need one loop just to identify how many nodes will be seeded by a merged curve node

        # Loop over every curve we want to merge
        for curveName in self.extra_data["parentCurves"]:

            # Check if we have access to this curve
            if curveName in curveDict:

                # Loop over the parent curve nodes to propagate seeds
                for parentNodeID in range(curveDict[curveName].coor.shape[0]):

                    # Get corresponding ID in the child nodes
                    childNodeID = linkOld2New[indexOffset + parentNodeID]

                    # Increment number of parent nodes added to this child node
                    numMergedNodes[childNodeID] = numMergedNodes[childNodeID] + 1

                # Update offset for next curve
                indexOffset = indexOffset + curveDict[curveName].coor.shape[0]

        # Compute averages
        for childNodeID in range(self.coor.shape[0]):
            self.coorb[childNodeID, :] = self.coorb[childNodeID, :] / numMergedNodes[childNodeID]

        # Now we have another loop to distribute derivative seeds

        # Initialize index offset. The new curves cannot reference
        # nodes of the old ones.
        indexOffset = 0

        # Loop over every curve we want to merge
        for curveName in self.extra_data["parentCurves"]:

            # Check if we have access to this curve
            if curveName in curveDict:

                # Loop over the parent curve nodes to propagate seeds
                for parentNodeID in range(curveDict[curveName].coor.shape[0]):

                    # Get corresponding ID in the child nodes
                    childNodeID = linkOld2New[indexOffset + parentNodeID]

                    # Accumulate derivative seed into the parent node
                    # (we will take the average later on)
                    curveDict[curveName].coorb[parentNodeID, :] = (
                        curveDict[curveName].coorb[parentNodeID, :] + self.coorb[childNodeID, :]
                    )

                # Update offset for next curve
                indexOffset = indexOffset + curveDict[curveName].coor.shape[0]

    # =================================================================

    def condense_disconnect_curves(self, guessNode=0):

        """
        This function takes a curve defined by multiple discontinuous
        curves (in the FE connectivity sense) and tries to merge them in
        a single continuous curve. This is useful to heal curves exported by
        ICEM where it defines multiple overlapping bar elements.
        """

        # First we try to sort the curve to find how many disconnect FE sets define it
        sortedConn, dummy_map = tst.FEsort(self.barsConn.tolist())

        # Get number of disconnect curves
        numCurves = len(sortedConn)

        # Check the number of connected FE sets
        if numCurves == 1:
            # The sorting works just fine and it found a single curve!
            # There is no need to condense nodes
            print("")
            print("tsurf_component:condense_nodes")
            print("Curve is already composed by a single FE set.")
            print("There is no need to condense curves.")
            print("")
        else:

            # Initialize
            FEcurves = []
            # Create one curve object for each disconnect curve
            for (curveID, currConn) in enumerate(sortedConn):
                curveName = "FEcurve_%03d" % curveID
                FEcurves.append(TSurfCurve(curveName, self.coor, currConn))

            # Now we will chain the closer nodes to make a single curve.
            # We will always check both ends of the chain and append the node that is
            # closer to one of the ends.

            # Make a list out of all coordinates so we can pop nodes
            coorList = self.coor.tolist()

            # Start chain of new coordinates with an arbitrary node
            newCoor = [coorList.pop(guessNode)]

            # Loop to chain nodes
            while len(coorList) > 0:

                # Get coordinates of the first and last nodes of the chain
                firstNode = newCoor[0]
                lastNode = newCoor[-1]

                # Compute distance of all remaining nodes to the first node of the chain
                dist2first = np.linalg.norm(np.array(coorList) - np.array(firstNode), axis=1)

                # Compute distance of all remaining nodes to the first node of the chain
                dist2last = np.linalg.norm(np.array(coorList) - np.array(lastNode), axis=1)

                # Find if which end of the chain has the closest new node
                if np.min(dist2first) < np.min(dist2last):

                    # The first node of the chain has the closest point.
                    # So we will add a node before the first point

                    # Get id of the closest node
                    nextNodeID = np.argmin(dist2first)

                    # Pop the element from the remaining nodes list
                    nextNode = coorList.pop(nextNodeID)

                    # Insert the next node at the beginning of the list
                    newCoor.insert(0, nextNode)

                else:

                    # The last node of the chain has the closest point.
                    # So we will add a node after the last point

                    # Get id of the closest node
                    nextNodeID = np.argmin(dist2last)

                    # Pop the element from the remaining nodes list
                    nextNode = coorList.pop(nextNodeID)

                    # Insert the next node at the beginning of the list
                    newCoor.append(nextNode)

            # Get number of nodes
            numNodes = len(newCoor)

            # Convert nodal coordinates to array
            newCoor = np.array(newCoor)

            # Now we create bar element connectivity for these nodes
            newBarsConn = np.zeros((numNodes - 1, 2))
            newBarsConn[:, 0] = np.arange(numNodes - 1)
            newBarsConn[:, 1] = np.arange(numNodes - 1) + 1

            # Replace information of the current curve
            self.coor = newCoor
            self.barsConn = newBarsConn

    def shift_end_nodes(self, criteria="maxX", startPoint=None, curveObject=None):

        """
        This method will shift the finite element ordering of
        a periodic curve so that the first and last elements
        satisfy a given criteria.
        This is useful when we want to set boundary conditions for
        mesh generation in specific points of this periodic curve,
        such as a trailing edge.

        This method only works for periodic curves with ordered FE.

        Parameters
        ----------
        criteria: string
            criteria used to reorder FE data. Available options are:
            ['minX', 'minY', 'minZ', 'maxX', 'maxY', 'maxZ', 'startPoint', 'curve']

        Returns
        -------
        This function has no explicit outputs. It changes self.barsConn.

        Ney Secco 2016-08
        """

        # Initialize startPoint if necessary
        if startPoint is None:
            startPoint = np.zeros((3))

        # Get current coordinates and connectivities
        coor = self.coor
        barsConn = self.barsConn

        # Get number of elements
        nElem = barsConn.shape[0]

        # CHECKING INPUTS

        # First we check if the FE data is ordered and periodic
        for elemID in range(1, nElem):

            # Get node indices
            prevNodeID = barsConn[elemID - 1, 1]
            currNodeID = barsConn[elemID, 0]

            # Check if the FE data is ordered
            if prevNodeID != currNodeID:

                # Print warning
                print("WARNING: Could not shift curve because it has unordered FE data.")
                print("         Call FEsort first.")
                return

        # The first and last nodes should be the same to have a periodic curve.
        if barsConn[0, 0] != barsConn[-1, 1]:

            # Print warning
            print("WARNING: Could not shift curve because it is not periodic.")
            return

        # REORDERING

        # Reorder according to the given criteria. We need to find the
        # reference node (startNodeID) to reorder the FE data

        if criteria in ["maxX", "maxY", "maxZ"]:

            if criteria == "maxX":
                coorID = 0
            elif criteria == "maxY":
                coorID = 1
            elif criteria == "maxZ":
                coorID = 2

            # Get maximum value of the desired coordinate.
            startNodeID = np.argmax(coor[:, coorID])

        elif criteria in ["minX", "minY", "minZ"]:

            if criteria == "minX":
                coorID = 0
            elif criteria == "minY":
                coorID = 1
            elif criteria == "minZ":
                coorID = 2

            # Get maximum value of the desired coordinate.
            startNodeID = np.argmin(coor[:, coorID])

        elif criteria == "startPoint":

            # Compute the squared distances between the selected startPoint
            # the coordinates of the curve.
            deltas = coor - startPoint
            dist_2 = np.sum(deltas ** 2, axis=1)

            # Set the start node as the curve coordinate point closest to the
            # selected startPoint
            startNodeID = np.argmin(dist_2)

        elif criteria == "curve":

            # Find which node is the closest one to the reference curve
            startNodeID = closest_node(curveObject.coor, coor)

        # Now we look for an element that starts at the reference node
        startElemID = np.where(barsConn[:, 0] == startNodeID)[0][0]

        # Then we shift the connectivity list so that startElemID becomes the first element
        self.barsConn[:, :] = np.vstack([barsConn[startElemID:, :], barsConn[:startElemID, :]])

        # Do the same to the list of parent triangles if necessary
        if np.array(self.extra_data["parentTria"]).any():
            self.extra_data["parentTria"] = np.vstack(
                [self.extra_data["parentTria"][startElemID:, :], self.extra_data["parentTria"][:startElemID, :]]
            )

    def export_tecplot(self, outputName="curve"):

        """
        This function will export the current curve in tecplot FE format

        Ney Secco 2016-08
        """

        # Create tecplot curve object and export
        tecCurve = tecplot_interface.Curve(self.coor, self.barsConn, self.name)
        tecCurve.export_tecplot(outputName)

    # ===========================================================#
    # INTERFACE AND DERIVATIVE SEED MANIPULATION METHODS
    #
    # All these methods below assume that the FE data is ordered. That is
    # The last one of the previous bar in barsConn is the first node of the next bar.

    def get_points(self):

        """
        We will just return the nodes that define the curve.
        If the user wants another set of points, then he should use self.remesh first.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Initialize coordinates matrix
        pts = np.zeros((numElems + 1, 3))

        # Get the last point
        pts[-1, :] = self.coor[self.barsConn[-1, -1], :]

        # Get coordinates
        for ii in range(numElems):
            pts[ii, :] = self.coor[self.barsConn[ii, 0], :]

        # Return coordinates
        return pts

    def set_points(self, pts):

        """
        We will just return the nodes that define the curve.
        If the user wants another set of points, then he should use self.remesh first.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Set the last point
        self.coor[self.barsConn[-1, -1], :] = pts[-1, :]

        # Set other points
        for ii in range(numElems):
            self.coor[self.barsConn[ii, 0], :] = pts[ii, :]

    def set_forwardADSeeds(self, coord):
        """
        This will apply derivative seeds to the curve nodes.

        coord: array[numNodes, 3] -> Array with derivative seeds. They
        should follow the order of the curve nodes, which is not necessarily
        the self.coor ordering.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Get the last point
        self.coord[self.barsConn[-1, -1], :] = coord[-1, :]

        # Get coordinates
        for ii in range(numElems):
            self.coord[self.barsConn[ii, 0], :] = coord[ii, :]

    def get_forwardADSeeds(self):

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Initialize seed matrix
        coord = np.zeros((numElems + 1, 3))

        # Get the last point
        coord[-1, :] = self.coord[self.barsConn[-1, -1], :]

        # Get coordinates
        for ii in range(numElems):
            coord[ii, :] = self.coord[self.barsConn[ii, 0], :]

        # Return derivatives
        return coord

    def set_reverseADSeeds(self, coorb):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Reinitialize all seeds
        self.coorb[:, :] = 0.0

        # Get coordinates
        for ii in range(numElems):
            self.coorb[self.barsConn[ii, 0], :] = coorb[ii, :]

        # Add the last point
        self.coorb[self.barsConn[-1, -1], :] = self.coorb[self.barsConn[-1, -1], :] + coorb[-1, :]

    def get_reverseADSeeds(self, clean=True):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Initialize seed matrix
        coorb = np.zeros((numElems + 1, 3))

        # Get coordinates
        for ii in range(numElems):
            coorb[ii, :] = self.coorb[self.barsConn[ii, 0], :]

        # Add the last point
        coorb[-1, :] = self.coorb[self.barsConn[-1, -1], :]

        # Adjust seeds for periodic curves
        if self.isPeriodic:
            coorb[0, :] = coorb[0, :] / 2.0
            coorb[-1, :] = coorb[-1, :] / 2.0

        # Clean derivatives if necessary
        if clean:
            self.coorb[:, :] = 0.0

        # Return derivatives
        return coorb

    def accumulate_reverseADSeeds(self, coorb):
        """
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoorb: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # Get the number of elements
        numElems = self.barsConn.shape[0]

        # Get coordinates
        for ii in range(numElems):
            self.coorb[self.barsConn[ii, 0], :] = self.coorb[self.barsConn[ii, 0], :] + coorb[ii, :]

        # Add the last point (if the curve is periodic, this will accumulate the last point into the first one)
        self.coorb[self.barsConn[-1, -1], :] = self.coorb[self.barsConn[-1, -1], :] + coorb[-1, :]

    def clean_reverseADSeeds(self):
        """
        This will erase all reverse derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        self.coorb[:, :] = 0.0

    def set_randomADSeeds(self, mode="both", fixedSeed=True):

        """
        This will set random normalized seeds to all variables.
        This can be used for testing purposes.

        mode: ['both','forward','reverse'] -> Which mode should have
        its derivatives replaced.
        """

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # Set forward AD seeds
        if mode == "forward" or mode == "both":

            coord = np.random.random_sample(self.coor.shape)
            coord = coord / np.sqrt(np.sum(coord ** 2))
            self.coord = coord

        # Set reverse AD seeds
        if mode == "reverse" or mode == "both":

            coorb = np.random.random_sample(self.coor.shape)
            coorb = coorb / np.sqrt(np.sum(coorb ** 2))
            self.coorb = coorb

        # Return the randomly generated seeds

        if mode == "forward":
            return self.get_forwardADSeeds()

        elif mode == "reverse":
            return self.get_reverseADSeeds(clean=False)

        elif mode == "both":
            return self.get_forwardADSeeds(), self.get_reverseADSeeds(clean=False)

    # ===========================================================#
    # DERIVATIVE SEED MANIPULATION METHODS (unordered versions)
    #
    # These are variations of the original methods above.
    # The difference is that they operatly directly into coord and coorb, bypassing
    # the FE ordering specified in barsConn.
    # These methods are not necessary for other geometry engines.

    def _set_forwardADSeeds(self, coord=None):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        if coord is not None:
            self.coord = np.array(coord)

    def _get_forwardADSeeds(self):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        coord = np.array(self.coord)

        return coord

    def _set_reverseADSeeds(self, coorb=None):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        if coorb is not None:
            self.coorb = np.array(coorb)

    def _get_reverseADSeeds(self, clean=True):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # We use np.arrays to make hard copies, and also to enforce Fortran ordering.
        coorb = np.array(self.coorb)

        # Check if we need to clean the derivative seeds
        if clean:
            self.clean_reverseADSeeds()

        return coorb

    def _accumulate_reverseADSeeds(self, coorb=None):
        """
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoorb: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        if coorb is not None:
            self.coorb = self.coorb + coorb


### HELPER FUNCTIONS ###


def closest_node(guideCurveNodes, points):
    """Find which point in points is the closest to the set of nodes in guideCurveNodes."""
    minDist = 1e9
    points = np.asarray(points)
    for refNode in guideCurveNodes:
        deltas = points - refNode
        dist_2 = np.einsum("ij,ij->i", deltas, deltas)
        if minDist > np.min(dist_2):
            minDist = np.min(dist_2)
            ind = np.argmin(dist_2)
    return ind
