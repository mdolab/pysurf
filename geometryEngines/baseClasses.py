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

        # Define dictionary that will hold all geometries
        self.geoms = {}

        # Define dictionary to hold intersection curves
        self.intCurves = {}

        # Define a task list.
        # This list will store all tasks done during the forward pass so
        # that we could repeat the same steps when propagating derivatives.
        # Every element of this list contains a sub-list used to define the task.
        # The first element of the sub-list should be the task type, and the
        # remaining elements could be the arguments used by the task.
        # For instance, one sub-list could be: ['remesh',optionsDict], or
        # ['intersect',distTol].
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

    def clear_all(self):

        '''
        This method will clear all intersections and tasks of the current manager,
        so that it can be used once again from scratch.
        The geometry objects will remain.
        '''

        self.intCurves = {}
        self.tasks = []

    #=====================================================
    # AD METHODS

    def forwardAD(self):
        '''
        This step will execute forward AD for all steps stored in self.tasks.
        '''

        print ''
        print '================================================='
        print 'Starting forward AD pass'
        print ''

        # Get the number of tasks
        numTasks = int(len(self.tasks))

        # Execute forward AD for every task
        for taskID in range(numTasks):

            # Get the name and the arguments of the task
            task = self.tasks[taskID]
            taskName = task[0]
            taskArg = task[1:]

            print ''
            print 'forwardAD task'
            print taskName
            print taskArg
            print ''

            # Run the corresponding AD code
            if taskName == 'intersect':

                # Get arguments
                distTol = taskArg[0]
                geomList = taskArg[1]

                # Run the AD code
                self.intersect_d(distTol, geomList)

            if taskName == 'remesh_intCurve':

                # Get arguments
                newCurveName = taskArg[0]
                curveName = taskArg[1]
                optionsDict = taskArg[2]

                # Run the AD code
                self.remesh_intCurve_d(newCurveName, curveName, optionsDict)

            if taskName == 'split_intCurve':

                # Get arguments
                curveName = taskArg[0]
                childrenName = taskArg[1]

                # Run the AD code
                self.split_intCurve_d(curveName, childrenName)

            if taskName == 'merge_intCurves':

                # Get arguments
                curveNames = taskArg[0]
                mergedCurveName = taskArg[1]

                # Run the AD code
                self.merge_intCurves_d(curveNames, mergedCurveName)

        print ''
        print 'Finished forward AD pass'
        print '================================================='
        print ''

    def reverseAD(self):
        '''
        This step will execute reverse AD for all steps stored in self.tasks.
        '''

        print ''
        print '================================================='
        print 'Starting reverse AD pass'
        print ''

        # Get the number of tasks
        numTasks = int(len(self.tasks))

        # Execute reverse AD for every task (in reverse order)
        for taskID in reversed(range(numTasks)):

            # Get the name and the arguments of the task
            task = self.tasks[taskID]
            taskName = task[0]
            taskArg = task[1:]

            print ''
            print 'reverseAD task'
            print taskName
            print taskArg
            print ''

            # Run the corresponding AD code
            if taskName == 'intersect':

                # Get arguments
                distTol = taskArg[0]
                geomList = taskArg[1]

                # Run the AD code
                self.intersect_b(distTol, geomList)

            if taskName == 'remesh_intCurve':

                # Get arguments
                newCurveName = taskArg[0]
                curveName = taskArg[1]
                optionsDict = taskArg[2]

                # Run the AD code
                self.remesh_intCurve_b(newCurveName, curveName, optionsDict)

            if taskName == 'split_intCurve':

                # Get arguments
                curveName = taskArg[0]
                childrenName = taskArg[1]

                # Run the AD code
                self.split_intCurve_b(curveName, childrenName)

            if taskName == 'merge_intCurves':

                # Get arguments
                curveNames = taskArg[0]
                mergedCurveName = taskArg[1]

                # Run the AD code
                self.merge_intCurves_b(curveNames, mergedCurveName)

        print ''
        print 'Finished reverse AD pass'
        print '================================================='
        print ''

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

        # Initialize list of intersection curve names
        intCurveNames = []

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

                    # Add curve name to the list
                    intCurveNames.append(curve.name)

                    # Store name of parent components
                    curve.extra_data['parentGeoms'] = [name1, name2]

                    # Add curve to the manager object
                    self.add_curve(curve)

        # Print log
        print 'Computed',numCurves,'intersection curves.'

        # Save the current task and its argument
        self.tasks.append(['intersect', distTol, geomList])

        # Return the names of the intersection curves
        return intCurveNames

    def intersect_d(self, distTol, geomList):

        '''
        This method will execute the forward AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curve in self.intCurves.itervalues():

            # Get pointers to the parent objects
            if curve.extra_data['parentGeoms']:

                geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
                geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

                # Check if the current intersection was generated by the
                # current group of geometries
                if (geom1.name in geomList) and (geom2.name in geomList):
                    
                    # Run the AD intersection code
                    geom1.intersect_d(geom2, curve, distTol)

    def intersect_b(self, distTol, geomList, accumulateSeeds=True):

        '''
        This method will execute the reverse AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curve in self.intCurves.itervalues():

            # Get pointers to the parent objects
            if curve.extra_data['parentGeoms']:

                geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
                geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

                # Check if the current intersection was generated by the
                # current group of geometries
                if (geom1.name in geomList) and (geom2.name in geomList):
                    
                    # Run the AD intersection code
                    geom1.intersect_b(geom2, curve, distTol, accumulateSeeds)

    #=====================================================
    # REMESHING METHODS

    def remesh_intCurve(self, curveName, optionsDict={}):
        '''
        This method will remesh an intersection curve defined under the current
        manager object.

        optionsDict: should be a dictionary whose keys are the arguments of the
        remesh function used by the curve.remesh function. The keys are the values
        of these arguments.
        '''

        if curveName in self.intCurves.keys():

            newCurve = self.intCurves[curveName].remesh(**optionsDict)

            # Store information regarding the parent curve (the one that was remeshed to get the new curve)
            newCurve.extra_data['parentCurve'] = self.intCurves[curveName].name

            # Rename the new curve
            newCurveName = newCurve.name

            # Delete intersection history (since the new curve did not come from an intersection)
            newCurve.extra_data['parentGeoms'] = []

            # Add the new curve to the intersection list
            self.add_curve(newCurve)

            # Save task information
            self.tasks.append(['remesh_intCurve',newCurveName,curveName,optionsDict])

        else:
            raise NameError('Cannot remesh curve '+curveName+'. Curve not defined.')

        # Return the name of the new curve
        return newCurveName

    def remesh_intCurve_d(self, newCurveName, curveName, optionsDict):
        '''
        This method has the forward mode AD for the remesh code.

        newCurveName: string -> Name of the curve generated by the remesh code.

        curveName: string -> Name of the original curve used as input to the remesh code.

        optionsDict: should be a dictionary whose keys are the arguments of the
        remesh function used by the curve.remesh function. The keys are the values
        of these arguments.
        '''

        # Set flag to identify errors
        foundCurves = False

        if curveName in self.intCurves.keys():

            if newCurveName in self.intCurves.keys():

                # Set flag to identify errors
                foundCurves = True

                # Get derivative seeds
                coord = self.intCurves[curveName].get_forwardADSeeds()

                # Call AD code
                newCoord = self.intCurves[curveName].remesh_d(coord, **optionsDict)

                # Store propagated seeds
                self.intCurves[newCurveName].set_forwardADSeeds(newCoord)

        if not foundCurves:
            raise NameError('Cannot use remesh_d with curves '+curveName+' and '+newCurveName+'. Curves not defined.')

    def remesh_intCurve_b(self, newCurveName, curveName, optionsDict, clean=True, accumulateSeeds=True):
        '''
        This method has the reverse mode AD for the remesh code.

        newCurveName: string -> Name of the curve generated by the remesh code.

        curveName: string -> Name of the original curve used as input to the remesh code.

        optionsDict: should be a dictionary whose keys are the arguments of the
        remesh function used by the curve.remesh function. The keys are the values
        of these arguments.
        '''

        # Set flag to identify errors
        foundCurves = False

        if curveName in self.intCurves.keys():

            if newCurveName in self.intCurves.keys():

                # Set flag to identify errors
                foundCurves = True

                # Get derivative seeds
                newCoorb = self.intCurves[newCurveName].get_reverseADSeeds(clean=clean)

                # Call AD code
                coorb = self.intCurves[curveName].remesh_b(newCoorb, **optionsDict)

                # Store propagated seeds
                if accumulateSeeds:
                    self.intCurves[curveName].accumulate_reverseADSeeds(coorb)
                else:
                    self.intCurves[curveName].set_reverseADSeeds(coorb)

        if not foundCurves:
            raise NameError('Cannot use remesh_b with curves '+curveName+' and '+newCurveName+'. Curves not defined.')

    #=====================================================
    # SPLITTING METHODS

    def split_intCurve(self, curveName, optionsDict={}, criteria='sharpness'):
        '''
        This method will split a given intersection curve based on a certain criteria.
        The child curves will be added to the self.intCurves dictionary.

        curveName: string -> Name of the original curve that will be split
        '''
        
        if curveName in self.intCurves.keys():

            # Call split function
            splitCurvesDict = self.intCurves[curveName].split(optionsDict, criteria)

            # Add new curves to the manager's dictionary
            for curve in splitCurvesDict.itervalues():

                self.add_curve(curve)

            # Save this task
            self.tasks.append(['split_intCurve',curveName,splitCurvesDict.keys()])

        else:
            raise NameError('Cannot split curve '+curveName+'. Curve not defined.')

        # Return the names of the new curves
        return splitCurvesDict.keys()

    def split_intCurve_d(self, curveName, childrenNames):
        '''
        This method propagates forward AD seeds from the parent curve to its children curves.
        '''

        # Check if the parent curve is defined
        if curveName not in self.intCurves.keys():
            raise NameError('Cannot use split_intCurve_d with curve '+curveName+'. Curve not defined.')
        else:
            parentCurve = self.intCurves[curveName]

        # Loop over all children to propagate seeds
        for childName in childrenNames:

            # Check if this child actually belongs to this parent
            if childName in parentCurve.extra_data['splitCurves']:

                # Run AD code
                parentCurve.split_d(self.intCurves[childName])

    def split_intCurve_b(self, curveName, childrenNames):
        '''
        This method propagates reverse AD seeds from the children curve to its children curves.
        '''

        # Check if the parent curve is defined
        if curveName not in self.intCurves.keys():
            raise NameError('Cannot use split_intCurve_d with curve '+curveName+'. Curve not defined.')
        else:
            parentCurve = self.intCurves[curveName]

        # Loop over all children to propagate seeds
        for childName in childrenNames:

            # Check if this child actually belongs to this parent
            if childName in parentCurve.extra_data['splitCurves']:

                # Run AD code
                parentCurve.split_b(self.intCurves[childName])

    #=====================================================
    # MERGING METHODS

    def merge_intCurves(self, curveNames, mergedCurveName):

        '''
        This will merge all curves whose names are in curveNames

        curveNames: list of strings -> Name of the curves to be merged
        '''

        # Get the name of the first curve
        mainCurveName = curveNames[0]
        mainCurve = self.intCurves[mainCurveName]

        # Call the mesh function from the main curve
        mergedCurve = mainCurve.merge(self.intCurves, mergedCurveName, curveNames[1:])

        # Add the new curve to the manager's list
        self.add_curve(mergedCurve)

        # Save current task
        self.tasks.append(['merge_intCurves',curveNames,mergedCurveName])

    def merge_intCurves_d(self, curveNames, mergedCurveName):

        '''
        This will run forward mode AD to the merge process
        '''

        # Get pointer to the merged curve
        mergedCurve = self.intCurves[mergedCurveName]

        # Create dictionary with the parent curves
        curveDict = {}
        for curveName in curveNames:
            curveDict[curveName] = self.intCurves[curveName]

        # Call AD code
        mergedCurve.merge_d(curveDict)

    def merge_intCurves_b(self, curveNames, mergedCurveName):

        '''
        This will run reverse mode AD to the merge process
        '''

        # Get pointer to the merged curve
        mergedCurve = self.intCurves[mergedCurveName]

        # Create dictionary with the parent curves
        curveDict = {}
        for curveName in curveNames:
            curveDict[curveName] = self.intCurves[curveName]

        # Call AD code
        mergedCurve.merge_b(curveDict)
