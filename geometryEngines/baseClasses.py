from __future__ import division
import numpy as np
from mpi4py import MPI
import pysurf

class Geometry(object):

    '''
    This is the base Geometry class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.

    The required attributes for a Geometry object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.

    self.curves : dictionary{curveName:curveObject} -> Dictionary with component curves.
    '''


    def __init__(self, *arg, **kwargs):
        '''
        Call the initialization method defined in each
        derived class
        '''

        self.name = ''
        self.comm = None
        self.curves = {}

        if 'comm' in kwargs:
            self.comm = kwargs['comm']
        else:
            self.comm = MPI.COMM_WORLD

        # Define attribute to hold the geometry manipulation object
        self.manipulator = None

        # Define attribute to hold mesh object
        self.meshObj = None

        # Define attribute to hold structured surface mesh file name (we will use this to call pyHyp)
        self.meshFileName = None

        # Define attribute to hold structured surface mesh point set name (used by the geometry manipulator)
        self.meshPtSetName = None

        self._initialize(*arg, **kwargs)

    def _initialize(self, *arg):
        '''
        Virtual method
        '''
        pass

    def translate(self, x, y, z):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def scale(self, factor):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def rotate(self, angle, axis, point=None):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        INPUTS:
        angle: float -> Rotation angle, in degrees.

        axis: integer -> Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3] -> Coordinates of the rotation center. If point=None,
                           the origin should be used.
        '''
        pass

    #===========================================================#
    # SURFACE PROJECTION METHODS

    def project_on_surface(self, xyz):
        '''
        This function projects a point in space to the component surface.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def project_on_surface_d(self, xyz, xyzd):
        '''
        This is the forward AD mode of project_on_surface.
        '''
        pass

    def project_on_surface_b(self, xyz, xyzProjb, normProjb):
        '''
        This is the reverse AD mode of project_on_surface.
        '''
        pass

    #===========================================================#
    # CURVE PROJECTION METHODS

    def project_on_curve(self, xyz, curveCandidates=None):
        '''
        This function projects a point in space to the component curves.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def project_on_curve_d(self, xyz, xyzd, curveCandidates=None):
        '''
        This is the forward AD mode of project_on_curve.
        '''
        pass

    def project_on_curve_b(self, xyz, xyzProjb, tanProjb, curveCandidates=None):
        '''
        This is the reverse AD mode of project_on_curve.
        '''
        pass

    #===========================================================#
    # DERIVATIVE SEED MANIPULATION METHODS

    def set_forwardADSeeds(self, coord, curveCoord):
        '''
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coord = coord
        for curveName in self.curves:
            self.curves[curveName].coord = curveCoord[curveName]


    def get_forwardADSeeds(self):
        '''
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        coord = np.array(self.coord, order='F')

        curveCoord = {}
        for curveName in self.curves:
            curveCoord[curveName] = self.curves[curveName].get_forwardADSeeds()

        return coord, curveCoord


    def set_reverseADSeeds(self, coorb, curveCoorb):
        '''
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coorb = coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = curveCoord[curveName]


    def get_reverseADSeeds(self,clean=True):
        '''
        This will return the current derivative seeds of the surface design variables (coor)
        and curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition

        # We use np.arrays to make hard copies, and also to enforce Fortran ordering.
        coorb = np.array(self.coorb, order='F')

        curveCoorb = {}
        for curveName in self.curves:
            curveCoorb[curveName] = self.curves[curveName].get_reverseADSeeds(clean)

        # Check if we need to clean the derivative seeds
        if clean:
            self.coorb[:,:] = 0.0

        return coorb, curveCoorb


    def accumulate_reverseADSeeds(self, coorb=None, curveCoorb=None):
        '''
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        # This is just an example. The actual function should be overwritten
        # in the new geometry class definition
        self.coorb = self.coorb + coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = self.curves[curveName].coorb + curveCoord[curveName]

    #===========================================================#
    # MANIPULATOR INTERFACE METHODS

    def assign_manipulator(self, GMObj):

        '''
        This function assign a geometry manipulation object (such ad DVGeo) to the
        current Geometry object.

        INPUTS:

        GMObj: Geometry Manipulation Object

        Ney Secco 2017-02
        '''

        print ''
        print 'Assigning ',self.name,' to manipulator object'

        # Store the geometry manipulator object
        self.manipulator = GMObj

        # Call method to do specialized assignments according to the derived object.
        self._assign_manipulator()



        print 'Manipulator assignment finished'
        print ''

    def manipulator_update(self):

        '''
        This method will call the update functions from the manipulator object accordingly
        to update the geometry based on the new set of design variables.

        Some geometry manipulators (such as DVGeo) do not have a general update function to
        update all embedded points at once. So we need this separate class to call the update
        function as many times as needed, based on the internal structured of the derived
        geometry class.

        Ney Secco 2017-02
        '''

        pass

    def manipulator_forwardAD(self, xDVd, ptSetName):

        '''
        This function will use forward AD to propagate derivative seeds from
        the manipulator design variables to the triangulated surface (and associated
        curves) nodes.

        Ney Secco 2017-02
        '''

        pass

        #xSd = self.manipulator.totalSensitivityProd(xDVd, ptSetName)

    def manipulator_reverseAD(self, xDVb, ptSetName):

        '''
        This function will use reverse AD to propagate derivative seeds from
        the triangulated surface (and associated
        curves) nodes to the manipulator design variables.

        Ney Secco 2017-02
        '''

        pass

        #xDVb = self.manipulator.totalSensitivity(ptSetName)

    def manipulator_getDVs(self):

        '''
        This returns the design variables of the current manipulator.
        '''

        dvDict = self.manipulator.getValues()

        return dvDict

    #===========================================================#
    # MESH METHODS

    def assign_structured_surface_mesh(self, fileName):

        '''
        This method reads a CGNS or Plot3D file that contains just a structured
        SURFACE mesh and assigns its nodes to the current geometry.

        INPUTS:
        fileName: string -> It should specify the file name that
        contains the surface coordinates. If the file extension is
        ".cgns" we will use the CGNS reader. Otherwise, we will treat it
        as a formatted Plot3D file.

        Ney Secco 2017-02
        '''

        # Add spaces for printing statements
        print ''

        # Save filename
        self.meshFileName = fileName

        # Initialize mesh object
        print 'Assigning surface mesh to ',self.name
        print 'mesh file:',fileName
        self.meshObj = pysurf.SurfaceMesh('mesh', fileName)

        # Print statements
        print 'Done'
        print ''

    def embed_mesh_points(self):

        '''
        This method will assign the structured surface mesh points to the
        geometry manipulator. The manipulator and the mesh should be assigned
        BEFORE this method is called.
        '''

        # Add spaces for printing statements
        print ''

        # Generate name for the mesh point set
        ptSetName = self.name + ':surfMesh'
        self.meshPtSetName = ptSetName

        # Get surface mesh coordinates
        coor = self.meshObj.get_points()

        # Embed coordinates into the manipulator
        print 'Assigning surface mesh of ',self.name,' to the manipulator'
        self.manipulator.addPointSet(coor, self.meshPtSetName)
        print 'Done'

        # Add spaces for printing statements
        print ''

class Curve(object):


    '''
    This is the base Curve class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.
    Remember that a Geometry object may hold Curve objects.

    The required attributes for a Curve object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.
    '''

    def __init__(self, *arg, **kwargs):
        '''
        Call the initialization method defined in each
        child class
        '''
        self._initialize(*arg)

    def get_points(self):
        '''
        This method should return ordered coordinates along the curve in
        a 3 x numNodes array
        If a curve is periodic, a node should be repeated at the beginning and
        the end of this array.
        '''
        pass

    def update_dvs(self, coor):
        self.coor = coor

    def flip(self):
        '''
        This flip should flit the curve direction.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def translate(self, xyz):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def scale(self, factor):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def rotate(self, angle, axis, point=None):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        INPUTS:
        angle: float -> Rotation angle, in degrees.

        axis: integer -> Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3] -> Coordinates of the rotation center. If point=None,
                           the origin should be used.
        '''
        pass

    def project(self,xyz):
        pass

