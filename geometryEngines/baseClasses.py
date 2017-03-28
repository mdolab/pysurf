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

        # Only the root proc will work here
        if self.myID == 0:

            print ''
            print 'Assigning ',self.name,' to manipulator object'

            # Store the geometry manipulator object
            self.manipulator = GMObj

            # Generate name for the surface point set
            ptSetName = self.name + ':surfNodes'

            # Assing the triangulated surface nodes to the geometry manipulator
            print 'Assigning surface nodes from ',self.name,' to the manipulator'
            coor = self.get_points()
            self.manipulator.addPointSet(coor, ptSetName)
            print 'Done'

            # Store the set name for future uses
            self.ptSetName = ptSetName

            # Now we need to assign every curve to the manipulator as well
            for curveName in self.curves:

                # Generate name for the curve point set
                ptSetName = self.name + ':curveNodes:' + curveName

                # Assign curve nodes to the geometry manipulator
                print 'Assigning nodes from curve ',curveName,' to the manipulator'
                coor = self.curves[curveName].get_points()
                self.manipulator.addPointSet(coor, ptSetName)
                print 'Done'

                # Store the set name for future uses
                self.curves[curveName].ptSetName = ptSetName

            # Assign surface mesh node if we already have one
            if self.meshFileName is not None:
                self.embed_mesh_points()



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

        # Update surface nodes
        print 'Updating surface nodes from ',self.name
        coor = self.manipulator.update(self.ptSetName)
        self.set_points(coor)
        print 'Done'

        # Now we need to update every curve as well
        for curveName in self.curves:

            # Update curve nodes
            print 'Updating nodes from curve ',curveName
            coor = self.manipulator.update(self.curves[curveName].ptSetName)
            self.curves[curveName].set_points(coor)
            print 'Done'

        # Finally update the structured surface mesh, if we have one
        if self.meshPtSetName is not None:
            print 'Updating surface mesh from ',self.name
            self.meshObj.set_points(self.manipulator.update(self.meshPtSetName))
            print 'Done'

    def manipulator_forwardAD(self, xDVd, ptSetName=None):

        '''
        This method uses forward AD to propagate derivative seeds from the design
        variables to the surface mesh coordinates.

        xDVd : dictionary whose keys are the design variable names, and whose
               values are the derivative seeds of the corresponding design variable.
        '''

        # Print log
        print ''
        print 'Forward AD call to manipulator of ',self.name,' object'

        # Update surface nodes
        print 'Updating surface nodes from ',self.name
        coord = self.manipulator.totalSensitivityProd(xDVd, self.ptSetName)
        self.set_forwardADSeeds(coord=coord.reshape((-1,3)))
        #self.set_forwardADSeeds(coord=coord) # Use this if we fix DVGeo shape
        print 'Done'

        # Now we need to update every curve as well
        for curveName in self.curves:

            # Update curve nodes
            print 'Updating nodes from curve ',curveName
            coord = self.manipulator.totalSensitivityProd(xDVd, self.curves[curveName].ptSetName)
            self.curves[curveName].set_forwardADSeeds(coord.reshape((-1,3)))
            #self.curves[curveName].set_forwardADSeeds(coord)
            print 'Done'

        # Finally update the structured surface mesh, if we have one
        if self.meshPtSetName is not None:
            print 'Updating surface mesh from ',self.name
            coord = self.manipulator.totalSensitivityProd(xDVd, self.meshPtSetName)
            self.meshObj.set_forwardADSeeds(coord.reshape((-1,3)))
            #self.meshObj.set_forwardADSeeds(coord)
            print 'Done'

    def manipulator_reverseAD(self, xDVb, ptSetName=None, clean=True):

        '''
        This method uses reverse AD to propagate derivative seeds from the surface mesh
        coordinates to the design variables.

        xDVb : dictionary whose keys are the design variable names, and whose
               values are the derivative seeds of the corresponding design variable.
               This function will acuumulate the derivative seeds contained in xDVb.
        '''

        # Print log
        print ''
        print 'Reverse AD call to manipulator of ',self.name,' object'

        # Get derivative seeds
        if self.meshPtSetName is not None:
            coorb, curveCoorb, meshCoorb = self.get_reverseADSeeds(clean=clean)
        else:
            coorb, curveCoorb = self.get_reverseADSeeds(clean=clean)

        # Update surface nodes
        print 'Updating triangulated surface nodes from ',self.name
        xDVb_curr = self.manipulator.totalSensitivity(coorb, self.ptSetName)
        accumulate_dict(xDVb, xDVb_curr)
        print 'Done'

        # Now we need to update every curve as well
        for curveName in self.curves:

            # Update curve nodes
            print 'Updating nodes from curve ',curveName
            xDVb_curr = self.manipulator.totalSensitivity(curveCoorb[curveName], self.curves[curveName].ptSetName)
            accumulate_dict(xDVb, xDVb_curr)
            print 'Done'

        # Finally update the structured surface mesh, if we have one
        if self.meshPtSetName is not None:
            print 'Updating surface mesh from ',self.name
            xDVb_curr = self.manipulator.totalSensitivity(meshCoorb, self.meshPtSetName)
            accumulate_dict(xDVb, xDVb_curr)
            print 'Done'

    def manipulator_getDVs(self):

        '''
        This returns the design variables of the current manipulator.
        '''

        dvDict = self.manipulator.getValues()

        # Make sure all entries are real (sometimes DVGeo returns complex values)
        for key in dvDict:
            dvDict[key] = np.real(dvDict[key])

        return dvDict

    #===========================================================#
    # MESH METHODS

    def assign_structured_surface_mesh(self, fileName, extrusionOptions={}):

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

        # Only the root proc will work here
        if self.myID == 0:

            # Add spaces for printing statements
            print ''

            # Save filename
            self.meshFileName = fileName

            # Initialize mesh object
            print 'Assigning surface mesh to ',self.name
            print 'mesh file:',fileName
            self.meshObj = pysurf.SurfaceMesh('mesh', fileName, extrusionOptions)

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
        # Initialize dictionary to hold extra information
        self.extra_data = {}

        self._initialize(*arg)

        # Initialize important extra data fields
        self.extra_data['parentGeoms'] = None # This will be used by the manager to indentify intersections meshes
        self.extra_data['childMeshes'] = None # This will be used by the manager to indentify collar meshes

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

#==================================================
# AUXILIARY FUNCTIONS

def accumulate_dict(majorDict, minorDict):

    '''
    This function will loop over all keys of minorDict. If majorDict
    shares the same key, then we will accumulate (add) the values into
    majorDict.
    '''

    for key in minorDict:
        if key in majorDict.keys():
            majorDict[key] = majorDict[key] + minorDict[key]
