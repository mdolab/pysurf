from __future__ import division
import numpy as np
from mpi4py import MPI

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

    def set_forwardADSeeds(self, coord, curveCoord):
        '''
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        self.coord = coord
        for curveName in self.curves:
            self.curves[curveName].coord = curveCoord[curveName]

    def set_reverseADSeeds(self, coorb, curveCoorb):
        '''
        This will apply derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        self.coorb = coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = curveCoord[curveName]

    def accumulate_reverseADSeeds(self, coorb=None, curveCoorb=None):
        '''
        This will accumulate derivative seeds to the surface design variables (coor)
        and the curve design variables (curve.coor).

        curveCoord: Is a dictionary whose keys are curve names, and whose fields
        are derivative seeds to be applied on the control points of each curve.
        '''

        self.coorb = self.coorb + coorb
        for curveName in self.curves:
            self.curves[curveName].coorb = self.curves[curveName].coorb + curveCoord[curveName]

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

    def extract_points(self):
        '''
        This method should return coordinates along the curve.
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

class Manager(object):

    '''
    This is a geometry manager object. The geometry engines DO NOT NEED a
    derived Manager class.
    This Manager class will manage information passing among the several
    geometries and intersections curves defined in the problem.
    '''

    def __init__(self):

        self.geoms = {}
        self.intCurves = {}

        # Define a task list.
        # This list will store all tasks done during the forward pass so
        # that we could repeat the same steps when propagating derivatives.
        # The even indices of self.tasks have the task type, while the even
        # indices contain arguments used by the task that precedes it.
        self.tasks = []

        pass

    def add_geometry(self, geom):

        '''
        This method adds a Geometry object to the current Manager's dictionary.
        '''

        # Check if we already have a curve of same name
        if geom.name in self.geoms.keys():
            raise NameError('Trying to add geometry of same name.')

        self.geoms[geom.name] = geom

    def remove_geometry(self, geomName):

        '''
        This method removes a Geometry object from the current Manager's dictionary.
        '''

        del self.geoms[geomName]

    def add_curve(self, curve):

        '''
        This method adds a Curve object to the current Manager's dictionary.
        In general, this should be an intersection curve.
        '''

        # Check if we already have a curve of same name
        if curve.name in self.intCurves.keys():
            raise NameError('Trying to add curves of same name.')

        self.intCurves[curve.name] = curve

    def remove_curve(self, curveName):

        '''
        This method removes a Curve object from the current Manager's dictionary.
        '''

        del self.intCurves[curveName]

    #=====================================================
    # AD METHODS

    def forwardAD(self):
        '''
        This step will execute forward AD for all steps stored in self.tasks.
        '''

        # Get the number of tasks
        numTasks = int(len(self.tasks)/2)

        # Execute forward AD for every task
        for taskID in range(numTasks):

            # Get the name and the arguments of the task
            taskName = self.tasks[2*taskID]
            taskArg = self.tasks[2*taskID+1]

            # Run the corresponding AD code
            if taskName == 'intersect':

                # Get arguments
                distTol = taskArg

                # Run the intersection code
                self.intersect_d(distTol)

    def reverseAD(self):
        '''
        This step will execute reverse AD for all steps stored in self.tasks.
        '''

        # Get the number of tasks
        numTasks = int(len(self.tasks)/2)

        # Execute reverse AD for every task (in reverse order)
        for taskID in reversed(range(numTasks)):

            # Get the name and the arguments of the task
            taskName = self.tasks[2*taskID]
            taskArg = self.tasks[2*taskID+1]

            # Run the corresponding AD code
            if taskName == 'intersect':

                # Get arguments
                distTol = taskArg

                # Run the intersection code
                self.intersect_b(distTol)

    #=====================================================
    # INTERSECTION METHODS

    def intersect(self, geomList=None, distTol=1e-7):

        '''
        This method intersects all geometries contained in the current Manager,
        provided that their names are in geomList.
        All geometry objects should be of same type.

        if geomList==None, all geometries will be intersected.

        distTol is a distance tolerance used to merge nearby nodes when
        generating the intersection finite element data.
        '''

        # Generate list of geometry names if user provided None
        if geomList == None:
            geomList = self.geoms.keys()

        # Make list of geometry objects
        geomObjList = []

        for geomName in self.geoms:
            
            # Detect if the user want to use the current geometry
            if geomName in geomList:

                # Add the corresponding geometry object to the list
                geomObjList.append(self.geoms[geomName])

        # Get number of components
        numGeometries = len(geomObjList)

        # Stop if user gives only one component
        if numGeometries < 2:
            print 'ERROR: Cannot compute intersections with just one component'
            quit()

        # Initialize number of curves computed so far
        numCurves = 0

        # Call intersection function for each pair
        for ii in range(numGeometries):
            for jj in range(ii+1,numGeometries):

                # Gather names of parent components
                name1 = geomObjList[ii].name
                name2 = geomObjList[jj].name

                # Compute new intersections for the current pair
                newIntersections = geomObjList[ii].intersect(geomObjList[jj],distTol)

                # Append new curve objects to the dictionary
                for curve in newIntersections:

                    # Increment curve counter
                    numCurves = numCurves+1

                    # Store name of parent components
                    curve.extra_data['parentGeoms'] = [name1, name2]

                    # Add curve to the manager object
                    self.add_curve(curve)

        # Print log
        print 'Computed',numCurves,'intersection curves.'

        # Save the current task and its argument
        if 'intersect' not in self.tasks:
            self.tasks.append('intersect')
            self.tasks.append(distTol)

    def intersect_d(self, distTol):

        '''
        This method will execute the forward AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curve in self.intCurves.itervalues():

            # Get pointers to the parent objects
            geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
            geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

            # Run the AD intersection code
            geom1.intersect_d(geom2, curve, distTol)

    def intersect_b(self, distTol, accumulateSeeds=True):

        '''
        This method will execute the reverse AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curve in self.intCurves.itervalues():

            # Get pointers to the parent objects
            geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
            geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

            # Run the AD intersection code
            geom1.intersect_b(geom2, curve, distTol, accumulateSeeds)
