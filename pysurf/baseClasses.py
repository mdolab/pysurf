import numpy as np
from mpi4py import MPI


class Geometry(object):

    """
    This is the base Geometry class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.

    The required attributes for a Geometry object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.

    self.curves : dictionary{curveName:curveObject} -> Dictionary with component curves.
    """

    def __init__(self, *arg, **kwargs):
        """
        Call the initialization method defined in each
        derived class.

        Expected arguments for the base class are:
        comm : MPI Communicator
        """

        self.name = ""
        self.comm = None
        self.curves = {}

        if "comm" in kwargs:
            self.comm = kwargs["comm"]
        else:
            self.comm = MPI.COMM_WORLD

        # Assign proc ID
        self.myID = self.comm.Get_rank()

        # Define attribute to hold the geometry manipulation object
        self.manipulator = None

        # Define attribute to hold point sets embedded into the manipulator
        self.manipulatorPts = {}
        self.manipulatorPtsd = {}
        self.manipulatorPtsb = {}

        # Define point set name to store geometry nodes into the manipulator
        self.ptSetName = None

        self._initialize(*arg, **kwargs)

    def _initialize(self, *arg):
        """
        Virtual method
        """
        pass

    def translate(self, x, y, z):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def scale(self, factor):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def rotate(self, angle, axis, point=None):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        Parameters
        ----------
        angle: float
            Rotation angle, in degrees.

        axis: integer
            Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3]
            Coordinates of the rotation center. If point=None, the origin should be used.
        """
        pass

    # ===========================================================#
    # SURFACE PROJECTION METHODS

    def project_on_surface(self, xyz):
        """
        This function projects a point in space to the component surface.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def project_on_surface_d(self, xyz, xyzd):
        """
        This is the forward AD mode of project_on_surface.
        """
        pass

    def project_on_surface_b(self, xyz, xyzProjb, normProjb):
        """
        This is the reverse AD mode of project_on_surface.
        """
        pass

    # ===========================================================#
    # CURVE PROJECTION METHODS

    def project_on_curve(self, xyz, curveCandidates=None):
        """
        This function projects a point in space to the component curves.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def project_on_curve_d(self, xyz, xyzd, curveCandidates=None):
        """
        This is the forward AD mode of project_on_curve.
        """
        pass

    def project_on_curve_b(self, xyz, xyzProjb, tanProjb, curveCandidates=None):
        """
        This is the reverse AD mode of project_on_curve.
        """
        pass

    # ===========================================================#
    # DERIVATIVE SEED MANIPULATION METHODS

    def set_forwardADSeeds(self, coord, curveCoord):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coord = coord
        for curveName in self.curves:
            self.curves[curveName].coord = curveCoord[curveName]

    def get_forwardADSeeds(self):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        coord = np.array(self.coord, order="F")

        curveCoord = {}
        for curveName in self.curves:
            curveCoord[curveName] = self.curves[curveName].get_forwardADSeeds()

        return coord, curveCoord

    def set_reverseADSeeds(self, coorb, curveCoorb):
        """
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coorb = coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = curveCoorb[curveName]

    def get_reverseADSeeds(self, clean=True):
        """
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition

        # We use np.arrays to make hard copies, and also to enforce Fortran ordering.
        coorb = np.array(self.coorb, order="F")

        curveCoorb = {}
        for curveName in self.curves:
            curveCoorb[curveName] = self.curves[curveName].get_reverseADSeeds(clean)

        # Check if we need to clean the derivative seeds
        if clean:
            self.coorb[:, :] = 0.0

        return coorb, curveCoorb

    def accumulate_reverseADSeeds(self, coorb=None, curveCoorb=None):
        """
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        """

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coorb = self.coorb + coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = self.curves[curveName].coorb + curveCoorb[curveName]

    # ===========================================================#
    # MANIPULATOR INTERFACE METHODS

    def assign_manipulator(self, GMObj):

        """
        This function assign a geometry manipulation object (such ad DVGeo) to the
        current Geometry object.

        This function receives a geometry manipulator object and then embed the points of
        the current geometry object into this manipulator

        Parameters
        ----------

        GMObj: Geometry Manipulation Object

        Ney Secco 2017-02
        """

        # Only the root proc will embed the pySurf nodes into the manipulator
        if self.myID == 0:

            print("")
            print("Assigning ", self.name, " to manipulator object")

            # Store the geometry manipulator object
            self.manipulator = GMObj

            # Generate name for the surface point set
            ptSetName = self.name + ":surfNodes"

            # Assing the triangulated surface nodes to the geometry manipulator
            print("Assigning surface nodes from ", self.name, " to the manipulator")
            coor = self.get_points()
            self.manipulator.addPointSet(coor, ptSetName)
            print("Done")

            # Store the set name for future uses
            self.ptSetName = ptSetName

            # Now we need to assign every curve to the manipulator as well
            for curveName in self.curves:

                # Generate name for the curve point set
                ptSetName = self.name + ":curveNodes:" + curveName

                # Assign curve nodes to the geometry manipulator
                print("Assigning nodes from curve ", curveName, " to the manipulator")
                coor = self.curves[curveName].get_points()
                self.manipulator.addPointSet(coor, ptSetName)
                print("Done")

                # Store the set name for future uses
                self.curves[curveName].ptSetName = ptSetName

            print("Manipulator assignment finished")
            print("")

        else:

            # The other processors will just save the reference to the manipulator
            self.manipulator = GMObj

    def manipulator_addPointSet(self, coor, ptSetName):

        """
        This method adds an extra point set to the manipulator object.

        Ney Secco 2017-11
        """

        # First store the point set array
        self.manipulatorPts[ptSetName] = coor

        # Initialize arrays for derivative seeds
        self.manipulatorPtsd[ptSetName] = np.zeros(coor.shape)
        self.manipulatorPtsb[ptSetName] = np.zeros(coor.shape)

        # Now embed the points into the manipulator, under the same ptSetName
        self.manipulator.addPointSet(coor, ptSetName)

    def manipulator_update(self, ptSetName=None):

        """
        This method will call the update functions from the manipulator object accordingly
        to update the geometry based on the new set of design variables.

        Some geometry manipulators (such as DVGeo) do not have a general update function to
        update all embedded points at once. So we need this separate class to call the update
        function as many times as needed, based on the internal structured of the derived
        geometry class.

        If ptSetName==None, we only update the triangulated surface description. Otherwise,
        we update the triangulated surface and also the points in under the ptSetName.

        Ney Secco 2017-02
        """

        # Only the root proc should update the triangulated surfaces and curves
        if self.myID == 0:

            # Update surface nodes
            print("Updating triangulated surface nodes from ", self.name)
            coor = self.manipulator.update(self.ptSetName)
            self.set_points(coor)
            print("Done")

            # Now we need to update every curve as well
            for curveName in self.curves:

                # Update curve nodes
                print("Updating nodes from curve ", curveName)
                coor = self.manipulator.update(self.curves[curveName].ptSetName)
                self.curves[curveName].set_points(coor)
                print("Done")

        # Finally, all procs update the embeded nodes, if they have one
        if ptSetName is not None:

            if self.myID == 0:
                print("Updating embedded nodes from ", self.name, " under ptSet ", ptSetName)

            coor = self.manipulator.update(ptSetName)
            self.manipulatorPts[ptSetName] = coor

            if self.myID == 0:
                print("Done")

    def manipulator_forwardAD(self, xDVd, ptSetName=None):

        """
        This method uses forward AD to propagate derivative seeds from the design
        variables to the surface mesh coordinates.

        xDVd : dictionary whose keys are the design variable names, and whose
               values are the derivative seeds of the corresponding design variable.
        """

        # Print log
        if self.myID == 0:
            print("")
            print("Forward AD call to manipulator of ", self.name, " object")

            # Update surface nodes
            print("Updating triangulated surface nodes from ", self.name)
            coord = self.manipulator.totalSensitivityProd(xDVd, self.ptSetName)
            self.set_forwardADSeeds(coord=coord.reshape((-1, 3)))
            # self.set_forwardADSeeds(coord=coord) # Use this if we fix DVGeo shape
            print("Done")

            # Now we need to update every curve as well
            for curveName in self.curves:

                # Update curve nodes
                print("Updating nodes from curve ", curveName)
                coord = self.manipulator.totalSensitivityProd(xDVd, self.curves[curveName].ptSetName)
                self.curves[curveName].set_forwardADSeeds(coord.reshape((-1, 3)))
                # self.curves[curveName].set_forwardADSeeds(coord)
                print("Done")

        # Finally, all procs update the embeded nodes, if they have one
        if ptSetName is not None:

            if self.myID == 0:
                print("Updating embedded nodes from ", self.name, " under ptSet ", ptSetName)

            coord = self.manipulator.totalSensitivityProd(xDVd, ptSetName)
            self.manipulatorPtsd[ptSetName] = coord.reshape((-1, 3))

            if self.myID == 0:
                print("Done")

    def manipulator_reverseAD(self, xDVb, ptSetName=None, comm=None, clean=True):

        """
        This method uses reverse AD to propagate derivative seeds from the surface mesh
        coordinates to the design variables.

        xDVb : dictionary whose keys are the design variable names, and whose
               values are the derivative seeds of the corresponding design variable.
               This function will accumulate the derivative seeds contained in xDVb.
        """

        if self.myID == 0:

            # Print log
            print("")
            print("Reverse AD call to manipulator of ", self.name, " object")

            # Get derivative seeds
            coorb, curveCoorb = self.get_reverseADSeeds(clean=clean)

            # Update surface nodes
            print("Updating triangulated surface nodes from ", self.name)
            xDVb_curr = self.manipulator.totalSensitivity(coorb, self.ptSetName)
            accumulate_dict(xDVb, xDVb_curr)
            print("Done")

            # Now we need to update every curve as well
            for curveName in self.curves:

                # Update curve nodes
                print("Updating nodes from curve ", curveName)
                xDVb_curr = self.manipulator.totalSensitivity(curveCoorb[curveName], self.curves[curveName].ptSetName)
                accumulate_dict(xDVb, xDVb_curr)
                print("Done")

        # Finally, all procs update the embeded nodes, if they have one
        if ptSetName is not None:

            if self.myID == 0:
                print("Updating embedded nodes from ", self.name, " under ptSet ", ptSetName)

            xDVb_curr = self.manipulator.totalSensitivity(self.manipulatorPtsb[ptSetName], ptSetName)
            accumulate_dict(xDVb, xDVb_curr)

            if self.myID == 0:
                print("Done")

    def manipulator_getDVs(self):

        """
        This returns the design variables of the current manipulator.
        """

        dvDict = self.manipulator.getValues()

        # Make sure all entries are real (sometimes DVGeo returns complex values)
        for key in dvDict:
            dvDict[key] = np.real(dvDict[key])

        return dvDict


class Curve(object):

    """
    This is the base Curve class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.
    Remember that a Geometry object may hold Curve objects.

    The required attributes for a Curve object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.
    """

    def __init__(self):
        """
        Call the initialization method defined in each
        child class
        """
        # Initialize dictionary to hold extra information
        self.extra_data = {}

        # Initialize important extra data fields
        self.extra_data["parentGeoms"] = None  # This will be used by the manager to identify intersections meshes
        self.extra_data["childMeshes"] = None  # This will be used by the manager to identify collar meshes

    def get_points(self):
        """
        This method should return ordered coordinates along the curve in
        a 3 x numNodes array
        If a curve is periodic, a node should be repeated at the beginning and
        the end of this array.
        """
        pass

    def update_dvs(self, coor):
        self.coor = coor

    def flip(self):
        """
        This flip should flit the curve direction.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def translate(self, xyz):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def scale(self, factor):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        """
        pass

    def rotate(self, angle, axis, point=None):
        """
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        Parameters
        ----------
        angle: float
            Rotation angle, in degrees.

        axis: integer
            Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3]
            Coordinates of the rotation center. If point=None, the origin should be used.
        """
        pass

    def project(self, xyz):
        pass


# ==================================================
# AUXILIARY FUNCTIONS


def accumulate_dict(majorDict, minorDict):

    """
    This function will loop over all keys of minorDict. If majorDict
    shares the same key, then we will accumulate (add) the values into
    majorDict.
    """

    for key in minorDict:
        if key in list(majorDict.keys()):
            majorDict[key] = majorDict[key] + minorDict[key]
