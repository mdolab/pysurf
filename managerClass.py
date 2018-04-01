from __future__ import division
import numpy as np
from mpi4py import MPI
import pysurf
from collections import OrderedDict
import os
from cgnsutilities import cgns_utils as cs

class Manager(object):

    '''
    This is the pySurf manager object. The geometry engines DO NOT NEED a
    derived Manager class.
    This Manager class will manage information passing among the several
    geometries, intersections curves, and meshes defined in the problem.
    The user will use the manager to specify all geometry and meshing
    operations during the forward pass. The manager will store all these
    steps to properly execute the AD modes and compute derivatives.
    '''

    def __init__(self, distTol=1e-7, comm=None):

        # Set up MPI communicator if the user provided None
        if comm is None:
            comm = MPI.COMM_WORLD

        # Assign communicator
        self.comm = comm

        # Save ID of the current proc
        self.myID = self.comm.Get_rank()
        '''
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

        self.commSingle = newComm
        '''
        # Define dictionary that will hold all geometries
        self.geoms = OrderedDict()

        # Define dictionary to hold intersection curves
        self.intCurves = OrderedDict()

        # Define dictionary to hold mesh generators
        self.meshGenerators = {}

        # Define a task list.
        # This list will store all tasks done during the forward pass so
        # that we could repeat the same steps when propagating derivatives.
        # Every element of this list contains a sub-list used to define the task.
        # The first element of the sub-list should be the task type, and the
        # remaining elements could be the arguments used by the task.
        # For instance, one sub-list could be: ['remesh',optionsDict], or
        # ['intersect',distTol].
        self.tasks = []

        # MACH INTERFACE ATTRIBUTES

        # Set dictionary that will contain surface mesh points for different sets.
        self.points = OrderedDict()
        self.updated = {}

        # Store tolerance to merge and identify nodes
        self.distTol = distTol

        # Store solver node indices owned by each primary component and collar mesh
        self.solverPointIDs = {}

        # Store total number of solver points for each point set
        self.numSolverPts = OrderedDict()

        pass

    #---------

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

    #---------

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

    #---------

    def add_meshGenerator(self, meshGen):

        '''
        This method adds a Mesh Generator object to the current Manager's dictionary.
        '''

        self.meshGenerators[meshGen.name] = meshGen

    def remove_meshGenerator(self, meshGenName):

        '''
        This method removes a Mesh Generator object from the current Manager's dictionary.
        '''

        del self.meshGenerators[meshName]

    #---------

    def add_collarMesh(self, mesh):

        '''
        This method adds a Mesh object to the current Manager's dictionary.
        '''

        self.collarMeshes[mesh.name] = mesh

    def remove_collarMesh(self, meshName):

        '''
        This method removes a Mesh object from the current Manager's dictionary.
        '''

        del self.collarMeshes[meshName]

    #---------

    def clean_all(self):

        '''
        This method will clear all intersections, meshes, and tasks of the current manager,
        so that it can be used once again from scratch.
        The geometry objects will remain.
        '''

        self.intCurves = OrderedDict()
        self.meshGenerators = {}
        self.tasks = []

    #=====================================================
    # OPERATION METHODS

    def assign_baseFunction(self, baseFunction):

        '''
        This method assigns an user-defined operation function to the current manager.
        The user should define a function of the form:
        baseFunction(manager)
        This function should receive a manager object and then conduct all necessary
        operations to generate the collar surface meshes using this object. This includes all
        itersection, split, merge, and mesh marching calls.
        ATTENTION: The user should NOT call any mesh extrusion operation (volume mesh generation)
        within baseFunction. We only need the surface nodes of the collar meshes.

        This function will be used throughout the optimization to update the nodal
        coordinates of the collar meshes.

        INPUTS:

        baseFunction : function handle -> Handle to the function that performs the geometry operations.

        Ney Secco 2017-02
        '''

        # Assign function handle to the manager object
        self.baseFunction = baseFunction

    def run_baseFunction(self):

        '''
        This method will execute the base function to update the collar mesh coordiantes.

        ASSUMPTIONS:
        - We assume that the triangulated surfaces are up to date with respect to the design variables.
        This is usually guaranteed since this method is called from self.update.
        '''

        # Only the root proc will work here
        if self.myID == 0:

            # Clean previous data
            self.clean_all()
            self.clean_reverseADSeeds()

            # Call base function to operate on the manager itself
            self.baseFunction(self)

        # Now we need to send mesh information to the other procs.
        # The root proc will prepare a simple dictionary containing the relevant info
        # then the other procs will store this data in their intCurves field.
        # The curve data is used at addPointSet and getSurfacePoints.

        if self.myID == 0:

            # Initialize dictionary
            dummyMeshDict = OrderedDict()

            # Append data to the dictionary
            for meshName in self.meshGenerators:

                # Save mesh name in the dictionary with dummy data
                dummyMeshDict[meshName] = 0.0

        else:

            # Other procs just initialize variable to receive data
            dummyMeshDict = None

        # Broadcast the curve info to other procs
        dummyMeshDict = self.comm.bcast(dummyMeshDict, root=0)

        # Now the other procs store the data
        if self.myID != 0:
            self.meshGenerators = dummyMeshDict

    def initialize(self, directory, backgroundMeshInfo=None, fileNameTag=None, staticMeshFiles=None):

        '''
        This method will do the initialization step. This includes:

        - run base function to generate surface collar meshes
        - extrude all meshes with pyHyp
        - generate CGNS files with isolated meshes for pyWarp inputs
        - generate combined CGNS file with all meshes for ADflow

        Remember to use this outside of optimization scripts. This is useful to
        generate inputs that will be used to initialize other MACH modules during the
        actual optimization.

        INPUTS/OUTPUTS:

        Please refer to self.extrude_meshes to verify the inputs/outputs description.

        Ney Secco 2017-02
        '''

        # Add a slash if the directory does not have it
        if directory[-1] != '/':
            directory = directory + '/'

        # Only the root proc will work here
        if self.myID == 0:

            # Run base function to generate surface meshes
            self.run_baseFunction()

        # Extrude meshes
        combinedFileName = self.extrude_meshes(directory, backgroundMeshInfo, fileNameTag, staticMeshFiles)

        # Set symmetry planes to zero (only the root proc will use cgns_utils
        if self.myID == 0:

            # Import cgns_utils
            from cgnsutilities import cgns_utils as cs

            # Load the CGNS file
            grid = cs.readGrid(combinedFileName)

            # Force symmetry planes to zero
            grid.symmZero(backgroundMeshInfo['sym'])

            # Write updated mesh
            grid.writeToCGNS(combinedFileName)

        # Return combined file name to use in ADflow
        return combinedFileName

    def reinitialize(self, directory, fileNameTag=None):

        '''
        This method will do the pre-optimization step. This includes:

        - run base function to generate surface collar meshes
        - generate CGNS files with isolated meshes for pyWarp inputs
        - generate combined CGNS file with all meshes for ADflow

        Remember to use this at the beginning of an optimization script.

        INPUTS/OUTPUTS:

        Please refer to self.extrude_meshes to verify the inputs/outputs description.

        IMPORTANT: You should use the same fileNameTag used for the initialization!

        Ney Secco 2017-02
        '''

        # Add a slash if the directory does not have it
        if directory[-1] != '/':
            directory = directory + '/'

        # Run base function to generate surface meshes
        self.run_baseFunction()

        # Regenerate filename that has combined meshes
        combinedFileName = generate_combined_filename(directory, fileNameTag)

        # Return combined file name to use in ADflow
        return combinedFileName

    #=====================================================
    # AD METHODS

    def forwardAD(self):
        '''
        This step will execute forward AD for all steps stored in self.tasks.
        '''

        # Only the root proc will work here
        if self.myID == 0:

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
                    intCurveNames = taskArg[1]

                    # Run the AD code
                    self._intersect_d(distTol, intCurveNames)

                if taskName == 'remesh_intCurve':

                    # Get arguments
                    newCurveName = taskArg[0]
                    curveName = taskArg[1]
                    optionsDict = taskArg[2]

                    # Run the AD code
                    self._remesh_intCurve_d(newCurveName, curveName, optionsDict)

                if taskName == 'split_intCurve':

                    # Get arguments
                    curveName = taskArg[0]
                    childrenName = taskArg[1]

                    # Run the AD code
                    self._split_intCurve_d(curveName, childrenName)

                if taskName == 'merge_intCurves':

                    # Get arguments
                    curveNames = taskArg[0]
                    mergedCurveName = taskArg[1]

                    # Run the AD code
                    self._merge_intCurves_d(curveNames, mergedCurveName)

                if taskName == 'march_intCurve_surfaceMesh':

                    # Get arguments
                    curveName = taskArg[0]

                    # Run the AD code
                    self._march_intCurve_surfaceMesh_d(curveName)

            print ''
            print 'Finished forward AD pass'
            print '================================================='
            print ''

    def reverseAD(self):
        '''
        This step will execute reverse AD for all steps stored in self.tasks.
        '''

        # Only the root proc will work here
        if self.myID == 0:

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
                    intCurveNames = taskArg[1]

                    # Run the AD code
                    self._intersect_b(distTol, intCurveNames)

                if taskName == 'remesh_intCurve':

                    # Get arguments
                    newCurveName = taskArg[0]
                    curveName = taskArg[1]
                    optionsDict = taskArg[2]

                    # Run the AD code
                    self._remesh_intCurve_b(newCurveName, curveName, optionsDict)

                if taskName == 'split_intCurve':

                    # Get arguments
                    curveName = taskArg[0]
                    childrenName = taskArg[1]

                    # Run the AD code
                    self._split_intCurve_b(curveName, childrenName)

                if taskName == 'merge_intCurves':

                    # Get arguments
                    curveNames = taskArg[0]
                    mergedCurveName = taskArg[1]

                    # Run the AD code
                    self._merge_intCurves_b(curveNames, mergedCurveName)

                if taskName == 'march_intCurve_surfaceMesh':

                    # Get arguments
                    curveName = taskArg[0]

                    # Run the AD code
                    self._march_intCurve_surfaceMesh_b(curveName)

            print ''
            print 'Finished reverse AD pass'
            print '================================================='
            print ''

    def clean_reverseADSeeds(self):

        '''
        This function will clean reverse AD seeds of all objects associated with this manager.
        '''

        for geom in self.geoms.itervalues():
            geom.clean_reverseADSeeds()

        for intCurve in self.intCurves.itervalues():
            intCurve.clean_reverseADSeeds()

        if self.myID == 0:

            for meshGen in self.meshGenerators.itervalues():
                meshGen.meshObj.clean_reverseADSeeds()

    #=====================================================
    # INTERSECTION METHODS

    def intersect(self, geomList=None, distTol=None):

        '''
        This method intersects all geometries contained in the current Manager,
        provided that their names are in geomList.
        All geometry objects should be of same type.

        if geomList==None, all geometries will be intersected.

        distTol is a distance tolerance used to merge nearby nodes when
        generating the intersection finite element data. If the user provides None,
        then we use the default value of the manager class
        '''

        # Assign tolerance if user provided none
        if distTol is None:
            distTol = self.distTol

        # Only the root proc will work here
        if self.myID == 0:

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
                    newIntersections = geomObjList[ii].intersect(geomObjList[jj],distTol=distTol)

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
            self.tasks.append(['intersect', distTol, intCurveNames])

            # Return the names of the intersection curves
            return intCurveNames

    def _intersect_d(self, distTol, intCurveNames):

        '''
        This method will execute the forward AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curveName in intCurveNames:

            # Get current curve object
            curve = self.intCurves[curveName]

            # Get pointers to the parent objects
            geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
            geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

            # Run the AD intersection code
            geom1.intersect_d(geom2, curve, distTol)

    def _intersect_b(self, distTol, intCurveNames, accumulateSeeds=True):

        '''
        This method will execute the reverse AD code for every intersection
        curve.
        '''

        # Run the derivative code for every curve
        for curveName in intCurveNames:

            # Get current curve object
            curve = self.intCurves[curveName]

            # Get pointers to the parent objects
            geom1 = self.geoms[curve.extra_data['parentGeoms'][0]]
            geom2 = self.geoms[curve.extra_data['parentGeoms'][1]]

            # Run the AD intersection code
            geom1.intersect_b(geom2, curve, distTol, accumulateSeeds)

    #=====================================================
    # REMESHING METHODS

    def remesh_intCurve(self, curveName, optionsDict={}, inheritParentGeoms=True):
        '''
        This method will remesh an intersection curve defined under the current
        manager object.

        optionsDict: should be a dictionary whose keys are the arguments of the
        remesh function used by the curve.remesh function. The keys are the values
        of these arguments.

        inheritParentGeoms: This is a flag to indicate if the new curve should have
        the same parents as the original curve. This can make it easier to generate
        the surface meshes for intersections.
        '''

        # Only the root proc will work here
        if self.myID == 0:

            if curveName in self.intCurves.keys():

                newCurve = self.intCurves[curveName].remesh(**optionsDict)

                # Store information regarding the parent curve (the one that was remeshed to get the new curve)
                newCurve.extra_data['parentCurve'] = self.intCurves[curveName].name

                # Rename the new curve
                newCurveName = newCurve.name

                # Assign intersection history
                if inheritParentGeoms:
                    if self.intCurves[curveName].extra_data['parentGeoms'] is not None:
                        newCurve.extra_data['parentGeoms'] = self.intCurves[curveName].extra_data['parentGeoms'][:]
                    else:
                        newCurve.extra_data['parentGeoms'] = None
                else:
                    newCurve.extra_data['parentGeoms'] = []

                # Add the new curve to the intersection list
                self.add_curve(newCurve)

                # Save task information
                self.tasks.append(['remesh_intCurve',newCurveName,curveName,optionsDict])

            else:
                raise NameError('Cannot remesh curve '+curveName+'. Curve not defined.')

            # Return the name of the new curve
            return newCurveName

    def _remesh_intCurve_d(self, newCurveName, curveName, optionsDict):
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

                # Get pointer to the remeshed curve
                newCurve = self.intCurves[newCurveName]

                # Call AD code
                self.intCurves[curveName].remesh_d(newCurve, **optionsDict)

        if not foundCurves:
            raise NameError('Cannot use remesh_d with curves '+curveName+' and '+newCurveName+'. Curves not defined.')

    def _remesh_intCurve_b(self, newCurveName, curveName, optionsDict, clean=True, accumulateSeeds=True):
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

                # Get pointer to the remeshed curve
                newCurve = self.intCurves[newCurveName]

                # Call AD code
                self.intCurves[curveName].remesh_b(newCurve, clean, accumulateSeeds, **optionsDict)

        if not foundCurves:
            raise NameError('Cannot use remesh_b with curves '+curveName+' and '+newCurveName+'. Curves not defined.')

    #=====================================================
    # SPLITTING METHODS

    def split_intCurve(self, curveName, optionsDict={}, criteria='sharpness', inheritParentGeoms=True):
        '''
        This method will split a given intersection curve based on a certain criteria.
        The child curves will be added to the self.intCurves dictionary.

        curveName: string -> Name of the original curve that will be split

        inheritParentGeoms: boolean -> This is a flag to indicate if the new curve should have
        the same parents as the original curve. This can make it easier to generate
        the surface meshes for intersections.
        '''

        # Only the root proc will work here
        if self.myID == 0:

            if curveName in self.intCurves.keys():

                # Call split function
                splitCurvesDict = self.intCurves[curveName].split(optionsDict, criteria)

                # Add new curves to the manager's dictionary
                for curve in splitCurvesDict.itervalues():

                    # Assign parents if necessary
                    if inheritParentGeoms:
                        if self.intCurves[curveName].extra_data['parentGeoms'] is not None:
                            curve.extra_data['parentGeoms'] = self.intCurves[curveName].extra_data['parentGeoms'][:]
                        else:
                            curve.extra_data['parentGeoms'] = None

                    self.add_curve(curve)

                # Save this task
                self.tasks.append(['split_intCurve',curveName,splitCurvesDict.keys()])

            else:
                raise NameError('Cannot split curve '+curveName+'. Curve not defined.')

            # Return the names of the new curves
            return splitCurvesDict.keys()

    def _split_intCurve_d(self, curveName, childrenNames):
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

    def _split_intCurve_b(self, curveName, childrenNames):
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

    def merge_intCurves(self, curveNames, mergedCurveName, inheritParentGeoms=True):

        '''
        This will merge all curves whose names are in curveNames

        curveNames: list of strings -> Name of the curves to be merged

        inheritParentGeoms: boolean -> This is a flag to indicate if the new curve should have
        the same parents as the original curve. This can make it easier to generate
        the surface meshes for intersections. We will take the parents of the first
        curve that was merged
        '''

        # Only the root proc will work here
        if self.myID == 0:

            # Get the name of the first curve
            mainCurveName = curveNames[0]
            mainCurve = self.intCurves[mainCurveName]

            # Call the mesh function from the main curve
            mergedCurve = mainCurve.merge(self.intCurves, mergedCurveName, curveNames[1:])

            # Check if we need to inherit parent geometry surfaces
            if inheritParentGeoms:
                if mainCurve.extra_data['parentGeoms'] is not None:
                    mergedCurve.extra_data['parentGeoms'] = mainCurve.extra_data['parentGeoms'][:]
                else:
                    mergedCurve.extra_data['parentGeoms'] = None

            # Add the new curve to the manager's list
            self.add_curve(mergedCurve)

            # Save current task
            self.tasks.append(['merge_intCurves',curveNames,mergedCurveName])

    def _merge_intCurves_d(self, curveNames, mergedCurveName):

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

    def _merge_intCurves_b(self, curveNames, mergedCurveName):

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

    #=====================================================
    # SURFACE MESHING METHODS

    def march_intCurve_surfaceMesh(self, curveName, options0={}, options1={}, meshName=None):

        '''
        This method will generate the surface mesh on both sides of the
        intersection curve.
        '''

        # Only the root proc will work here
        if self.myID == 0:

            # Create a mesh name if the user provided none
            if meshName is None:
                meshName = 'mesh_'+curveName

            # Get pointer to the seed curve
            curve = self.intCurves[curveName]

            # Get geometry objects to march the mesh on
            parentGeoms = curve.extra_data['parentGeoms']

            if parentGeoms is None:
                raise NameError('The curve does not have parent geometries. Cannot march meshes. Try using manager.set_intCurve_parentGeoms')

            # Create hypsurf objects for the two sides of the mesh

            # Create first mesh
            meshGen0 = pysurf.hypsurf.HypSurfMesh(curve,
                                                  self.geoms[parentGeoms[0]],
                                                  options0,
                                                  meshName+'_0')

            #meshGen0.test_all()
            #quit()

            meshGen0.createMesh()

            # Flip the curve
            curve.flip()

            # March the second mesh
            meshGen1 = pysurf.hypsurf.HypSurfMesh(curve,
                                                  self.geoms[parentGeoms[1]],
                                                  options1,
                                                  meshName+'_1')

            meshGen1.createMesh()

            # Unflip the curve
            curve.flip()

            # Store meshes under the manager
            self.add_meshGenerator(meshGen0)
            self.add_meshGenerator(meshGen1)

            # Get names of the new meshes
            meshNames = [meshGen0.name, meshGen1.name]

            # Store these names into the curve object
            curve.extra_data['childMeshes'] = meshNames

            # Store this task info
            self.tasks.append(['march_intCurve_surfaceMesh',curveName])

            # Return names of the new meshes
            return meshNames

    def _march_intCurve_surfaceMesh_d(self, curveName):

        # Get pointers to the curve and meshes
        curve = self.intCurves[curveName]
        meshGen0 = self.meshGenerators[curve.extra_data['childMeshes'][0]]
        meshGen1 = self.meshGenerators[curve.extra_data['childMeshes'][1]]

        # Run AD code for the first mesh
        meshGen0.compute_forwardAD()

        # Flip the curve
        curve.flip()

        # Run AD code for the second mesh
        meshGen1.compute_forwardAD()

        # Unflip the curve
        curve.flip()

    def _march_intCurve_surfaceMesh_b(self, curveName):

        # Get pointers to the curve and meshes
        curve = self.intCurves[curveName]
        meshGen0 = self.meshGenerators[curve.extra_data['childMeshes'][0]]
        meshGen1 = self.meshGenerators[curve.extra_data['childMeshes'][1]]

        # Run AD code for the first mesh
        meshGen0.compute_reverseAD()

        # Flip the curve()
        curve.flip()

        # Run AD code for the second mesh
        meshGen1.compute_reverseAD()

        # Unflip the curve
        curve.flip()

    #=====================================================
    # MESH EXPORTATION METHODS

    def export_surface_meshes(self, directory, fileNameTag=None):

        '''
        This function will export all structured surface meshes into
        plot3d files. The files will be separated by primary geometries and
        also by collar meshes.

        INPUTS:

        directory: string -> Directory where will place all mesh files.

        fileNameTag: string -> Optional tag to append to the file names.
        '''

        # Only the root proc will work here
        if self.myID == 0:

            # Print log
            print ''
            print 'Exporting surface meshes'

            # Add a slash if the directory does not have it
            if directory[-1] != '/':
                directory = directory + '/'

            # Initialize counters
            primaryID = 0
            collarID = 0

            '''
            # First we will export the primary geometry meshes
            for geom in self.geoms.itervalues():

                # Check if the component has an associated surface mesh
                if geom.meshObj is not None:

                    # Generate file name
                    fileName = generate_primary_surface_mesh_filename(directory, geom.name, fileNameTag)

                    # Export mesh
                    geom.meshObj.export_plot3d(fileName)

                    # Increment counter
                    primaryID = primaryID + 1

                    # Print log
                    print 'Exported primary mesh for',geom.name
            '''

            # Initialize dictionary to hold family names.
            # pySurf will use family names to associate the meshes given by the solver
            familyList = {}

            # Now we will export the collar meshes
            for curve in self.intCurves.itervalues():

                # Verify if this curve was used to generate collar meshes
                if curve.extra_data['childMeshes'] is not None:

                    # Generate file name
                    fileName = generate_collar_surface_mesh_filename(directory, curve.name, fileNameTag)

                    # Merge different meshes that make up a single collar
                    meshList = []
                    for meshName in curve.extra_data['childMeshes']:

                        # Get reference to a mesh object with part of the collar
                        meshObj = self.meshGenerators[meshName].meshObj
                        meshList.append(meshObj)

                    mergedMesh = pysurf.mesh_tools.merge_meshes('mergedMesh', meshList)

                    # Export mesh
                    mergedMesh.export_plot3d(fileName)

                    # Increment counter
                    collarID = collarID + 1

                    # Save family names
                    familyList[curve.name] = {(ii+1):name for (ii,name) in enumerate(curve.extra_data['childMeshes'])}

                    # Print log
                    print '    Exported collar mesh for',curve.name

            # Print log
            print 'Exported all meshes!'
            print ''

        else:

            familyList = None

        # Send dictionary from root proc to all other procs
        familyList = self.comm.bcast(familyList, root=0)

        # Return family names that should be used for pyHyp extrusion
        return familyList

    def generate_default_pyWarpMulti_options(self, fileNameTag=None, extraZones=None):

        '''
        This method will give a dictionary with default set of options that can be used
        to initialize a pyWarpMulti instance for the current manager.

        directory and fileNameTag should be the same ones used during self.initialize.

        extraZones: list of strings -> Names of cgns zones that should be deformed but have
                    different names from the components and collar meshes. Nevertheless, remember
                    that this zone should have wall boundary conditions with family names that
                    match the component names.
        '''

        # First we gather names of primary components and intersection curves that have meshes

        # Initialize list
        zoneNames = []

        # Loop over all primary meshes to gather their names
        for geom in self.geoms.itervalues():

            # Append name to the list
            zoneNames.append(geom.name)

        # Loop over all collar meshes to gather their coordinates
        for curve in self.intCurves.itervalues():

            # Check if the current object has an associated mesh
            if curve.extra_data['childMeshes'] is not None:

                # Append name to the list
                zoneNames.append(curve.name)

        # Add extra zone names provided by the user
        if extraZones is not None:
            zoneNames.extend(extraZones)

        # Create default options for every name
        optionsDict = {}
        for zoneName in zoneNames:

            # Append nametag if necessary
            if fileNameTag is not None:
                zoneName = zoneName + '_' + fileNameTag

            optionsDict[zoneName] = {
                'warpType':'unstructured',
                'aExp': 3.0,
                'bExp': 5.0,
                'LdefFact':100.0,
                'alpha':0.25,
                'errTol':0.0001,
                'evalMode':'fast',
                'symmTol':1e-6,
                'useRotations':True,
                'bucketSize':8,
            }

        # The root processor will send its dictionary to everybody else
        optionsDict = self.comm.bcast(optionsDict, root=0)

        return optionsDict

    #=====================================================
    # GENERAL INTERFACE METHODS

    def getSurfacePoints(self, ptSetName):

        '''
        This function returns the surface mesh points of all meshes contained
        in the manager object (including both embeded points by the solver and collar meshes).

        ATTENTION: These points are in Solver ordering, so this can only be used after self.addPointSet

        INPUTS:

        ptSetName: string -> Name of the point set we want points from. This is the same used for self.addPointSet.
        '''

        # Loop over all primary meshes to gather their surface mesh coordinates.
        # Those coordinates are embedded in the FFDs and assign their to the
        # full solver points vector.
        for geom in self.geoms.itervalues():

            # Check if we assigned manipulators to this component
            if geom.name in self.solverPointIDs[ptSetName]:

                # Get indices of the full solver vector that correspond to this component
                geomPointIDs = self.solverPointIDs[ptSetName][geom.name]

                # Get points of the current mesh
                self.points[ptSetName][geomPointIDs] = geom.manipulatorPts[ptSetName]

        # Loop over every block of the collar mesh
        for meshName in self.meshGenerators:

            # Get points of the current collar mesh.
            # Just the root has the collar mesh points initially, so we need to
            # broadcast them to all procs.
            if self.myID == 0:
                collarManagerPts = self.meshGenerators[meshName].meshObj.get_points()
            else:
                collarManagerPts = None
            collarManagerPts = self.comm.bcast(collarManagerPts, root=0)

            # Get indices of the collar mesh points in the solver vector
            collarIDs = self.solverPointIDs[ptSetName][meshName]['collarPointIDs']

            # Extract slice of solver vector with the collar points
            # This is just to get the right size
            collarSolverPts = self.points[ptSetName][collarIDs,:]

            # Take the values from the manager collar and replace the
            # appropriate values of the solver collar
            self._convertManagerToSolver(collarSolverPts,
                                         collarManagerPts,
                                         self.solverPointIDs[ptSetName][meshName]['indexSolverPts'],
                                         self.solverPointIDs[ptSetName][meshName]['indexManagerPts'],
                                         self.solverPointIDs[ptSetName][meshName]['numSurfRepetition'])

            # Assign the points back to the full solver vector
            self.points[ptSetName][collarIDs,:] = collarSolverPts

        # Return the set of points
        return self.points[ptSetName]

    def getSurfaceForwardADSeeds(self, ptSetName):

        '''
        This function returns the forward AD seeds of all meshes contained
        in the manager object (including both primary and collar meshes).

        INPUTS:

        ptSetName: string -> Name of the point set we want points from. THis is the same used for self.addPointSet.

        RETURNS:

        ptsd: float[nPts,3] -> Forward derivative seeds of all surface points, following the solver ordering.
        '''

        # Initialize array of derivatives
        ptsd = np.zeros(self.points[ptSetName].shape)

        # Loop over all primary meshes to gather their surface mesh coordinates.
        # Those coordinates are embedded in the FFDs and assign their to the
        # full solver points vector.
        for geom in self.geoms.itervalues():

            # Check if we assigned manipulators to this component
            if geom.name in self.solverPointIDs[ptSetName]:

                # Get indices of the full solver vector that correspond to this component
                geomPointIDs = self.solverPointIDs[ptSetName][geom.name]

                # Get points of the current mesh
                ptsd[geomPointIDs,:] = geom.manipulatorPtsd[ptSetName]

        # Loop over every block of the collar mesh
        for meshName in self.meshGenerators:

            # Get points of the current collar mesh.
            # Just the root has the collar mesh points initially, so we need to
            # broadcast them to all procs.
            if self.myID == 0:
                collarManagerPtsd = self.meshGenerators[meshName].meshObj.get_forwardADSeeds()
            else:
                collarManagerPtsd = None
            collarManagerPtsd = self.comm.bcast(collarManagerPtsd, root=0)

            # Get indices of the collar mesh points in the solver vector
            collarIDs = self.solverPointIDs[ptSetName][meshName]['collarPointIDs']

            # Extract slice of solver vector with the collar points.
            # These values might be replaced by _convertManagerToSolver
            collarSolverPtsd = ptsd[collarIDs,:]

            # Take the values from the manager collar and replace the
            # appropriate values of the solver collar
            self._convertManagerToSolver(collarSolverPtsd,
                                         collarManagerPtsd,
                                         self.solverPointIDs[ptSetName][meshName]['indexSolverPts'],
                                         self.solverPointIDs[ptSetName][meshName]['indexManagerPts'],
                                         self.solverPointIDs[ptSetName][meshName]['numSurfRepetition'])

            # Assign the points back to the full solver vector
            ptsd[collarIDs,:] = collarSolverPtsd

        # Return the set of points
        return ptsd

    def setSurfaceReverseADSeeds(self, ptsb, ptSetName):

        '''
        This function sets the reverse AD seeds of all meshes contained
        in the manager object (including both primary and collar meshes).

        INPUTS:

        ptsb: float[nPts,3] -> Reverse derivative seeds of all surface points, following the solver ordering.

        ptSetName: string -> Name of the point set we want points from. THis is the same used for self.addPointSet.
        '''

        # Loop over all primary meshes to gather their surface mesh coordinates.
        # Those coordinates are embedded in the FFDs and assign their to the
        # full solver points vector.
        for geom in self.geoms.itervalues():

            # Check if we assigned manipulators to this component
            if geom.name in self.solverPointIDs[ptSetName]:

                # Get indices of the full solver vector that correspond to this component
                geomPointIDs = self.solverPointIDs[ptSetName][geom.name]

                # Get points of the current mesh
                geom.manipulatorPtsb[ptSetName][:,:] = ptsb[geomPointIDs,:]

                # Burn derivatives so that they are not used again
                ptsb[geomPointIDs,:] = 0.0

        # Loop over every block of the collar mesh
        for meshName in self.meshGenerators:

            # Get indices of the collar mesh points in the solver vector
            collarIDs = self.solverPointIDs[ptSetName][meshName]['collarPointIDs']

            # Extract slice of solver vector with the collar points.
            collarSolverPtsb = ptsb[collarIDs,:]

            # The root proc should send the collar mesh size to the other procs
            if self.myID == 0:
                collarSize = self.meshGenerators[meshName].meshObj.numPts
            else:
                collarSize = None
            collarSize = self.comm.bcast(collarSize, root=0)

            # Initialize array that will hold the derivative seeds.
            # These values might be replaced by _convertSolverToManagerb.
            collarManagerPtsb = np.zeros((collarSize,3))

            # Take the values from the manager collar and replace the
            # appropriate values of the solver collar.
            self._convertSolverToManagerb(collarSolverPtsb,
                                          collarManagerPtsb,
                                          self.solverPointIDs[ptSetName][meshName]['indexSolverPts'],
                                          self.solverPointIDs[ptSetName][meshName]['indexManagerPts'],
                                          self.solverPointIDs[ptSetName][meshName]['numSurfRepetition'])

            # Accumulate derivatives into the root proc
            collarManagerPtsb = self.comm.reduce(collarManagerPtsb, MPI.SUM, root=0)

            # After the conversion, the derivatives from collarSolverPtsb that were used are
            # replaced by zeros. So now we can eliminate the derivative seeds from the original
            # vector of seeds.
            ptsb[collarIDs,:] = collarSolverPtsb

            # Get points of the current collar mesh.
            # Just the root has the collar mesh points initially, so we need to
            # broadcast them to all procs.
            if self.myID == 0:
                self.meshGenerators[meshName].meshObj.set_reverseADSeeds(collarManagerPtsb)


    #=====================================================
    # MACH INTERFACE METHODS

    def addPointSet(self, coor, ptSetName, conn, faceSizes, famIDs, famName2FamID, origConfig=True, **kwargs):

        '''
        This function receives a set of coordinates from a solver (coor) and then makes the mapping
        between this set of coordinates and the coordinates currently owned by this manager.
        Therefore, for every surface update, we can recreate an array with the order expected by the solver.

        Parameters
        ----------
        coor : array, size (M, 3)
            Solver nodes on this processor
        ptSetName : string
            Name assigned to the set of points
        conn : int array, size  (sum(faceSizes))
            Connectivity of the nodes on this processor.
            All connectivities are assembled in a single flattened array,
            then we will use faceSizes to slice this array for each cell.
        faceSizes : int Array size (N)
            Treat the conn array as a flat list with the faceSizes giving
            connectivity offset for each element.
        famIDs : int Array size (N)
            family ID of each element.
        famName2FamID : Dictionary of strings to ints
            Dictionary to convert the family names to integers of famID.
            Each component names should match to one family name, since we will take
            the points of that family name and assign to the corresponding manipulator.

        Ney Secco 2017-11
        '''

        # Check if we already have this point set
        if ptSetName in self.points.keys():
            raise NameError('The point set',ptSetName,'is already defined under this manager.')

        # Print log
        if self.myID == 0:
            print ''
            print 'Adding point set',ptSetName,'to the manager.'

        # Store the total number of solver nodes
        self.numSolverPts[ptSetName] = coor.shape[0]

        # Store the coordinates given by the Solver in a separate dictionary
        self.points[ptSetName] = coor.copy()

        # First we need to flag the families of each point in coor.
        # We know the family of each cell, so we need to use the connectivities
        # to transfer this information to the points.

        # Create an array to hold the family ID of each point.
        pointFamIDs = np.zeros(self.numSolverPts[ptSetName], dtype='int64')

        # Initialize counter to slice the conn array
        startID = 0

        # Initialize dictionary to hold IDs of nodes that belong to each component
        self.solverPointIDs[ptSetName] = {}

        # Loop over the cells
        for cellID in np.arange(len(faceSizes)):

            # Get family name of the current cell
            cellFamID = famIDs[cellID]

            # Get size of the current cell
            faceSize = faceSizes[cellID]

            # Slice the conn array to get the point IDs of the current cell
            pointIDs = conn[startID:startID+faceSize]

            # Update counter for the next iteration
            startID = startID + faceSize

            # Loop over each point of the cell
            for pointID in pointIDs:

                # Assign same family of the cell to its points
                pointFamIDs[pointID] = cellFamID

        # Now we loop over each main component to find the solver points that correspond to each
        # one of them
        for geom in self.geoms.itervalues():

            # Check if the geometry object has a manipulator
            if geom.manipulator is not None:

                # Find points of the same family as geom
                geomPointIDs = np.where(pointFamIDs == famName2FamID[geom.name])[0]

                # Store IDs for future reference
                self.solverPointIDs[ptSetName][geom.name] = geomPointIDs

                # Assign the selected points to the manipulator object
                geom.manipulator_addPointSet(coor[geomPointIDs,:], ptSetName)

                # Print log
                if self.myID == 0:
                    print ' Assigned points from solver family "'+geom.name+'" to manipulator'

        # Now we need to find the collar mesh nodes
        # The family name used by the wall BC of the collar meshes should be the
        # same as the curve name that generated it.

        # Loop over every block of the collar mesh
        for meshName in self.meshGenerators:

            # Get IDs of the solver points that correspond to the collar mesh
            collarPointIDs = np.where(pointFamIDs == famName2FamID[meshName])[0]

            # Get solver points that have the collar family
            solverCollarPts = coor[collarPointIDs,:]

            # Get points of the current collar mesh.
            # Just the root has the collar mesh points initially, so we need to
            # broadcast them to all procs.
            if self.myID == 0:
                collarPts = self.meshGenerators[meshName].meshObj.get_points()
            else:
                collarPts = None
            collarPts = self.comm.bcast(collarPts, root=0)

            # Now we need to find which points given by the solver in this proc
            # correspond to the pySurf collar mesh nodes.
            indexSolverPts, indexManagerPts, numSurfRepetition = self._setSolverToManagerMapping(collarPts,
                                                                                                 solverCollarPts,
                                                                                                 distTol=None)

            # We also store the collar mesh mapping under self.solverPointIDs.
            # However, we need to store more data for the collar meshes.
            # So we use a dictionary for that.
            # collarPointIDs is used to extract the collar nodes from the full solver nodes array.
            # indexSolverPts has the indices of the extracted collar nodes that will map to each
            # element from the indexManagerPts.
            # numSurfRepetition is the number of repeated nodes in indexSolverPts.
            self.solverPointIDs[ptSetName][meshName] = {'collarPointIDs':collarPointIDs,
                                                        'indexSolverPts':indexSolverPts,
                                                        'indexManagerPts':indexManagerPts,
                                                        'numSurfRepetition':numSurfRepetition}

            # Print log
            if self.myID == 0:
                print ' Assigned points from solver family "'+meshName

        # Flag this point set as up to date
        self.updated[ptSetName] = True

        # Print log
        if self.myID == 0:
            print 'Done'
            print ''

    def getValues(self):

        '''
        This function returns a dictionary with the current design variable (DV) values.
        This can be used to get a baseline dictionary to assign new DV values with
        self.setDesignVars later on.

        ASSUMPTIONS:
        - We assume that different component do NOT have DVs with the same name.

        Ney Secco 2017-02
        '''

        # Initialize dictionary of design variables
        dvDict = {}

        # Loop over the primary geometries to find design variables
        for geom in self.geoms.itervalues():

            # Check if the geometry object has an associated geometry manipulator
            if geom.manipulator is not None:

                # Get design variables from the manipulator
                curr_dvDict = geom.manipulator_getDVs()

                # Loop over every entry of the new DV dictionary
                for key in curr_dvDict:

                    # Add new entries to the dictionary
                    dvDict[key] = curr_dvDict[key]

        # Return DV dictionary
        return dvDict

    def getVarNames(self):

        '''
        Returns a list of design variable names.
        '''

        # Initialize list of design variable names
        names = []

        # Loop over the primary geometries to find design variables
        for geom in self.geoms.itervalues():

            # Check if the geometry object has an associated geometry manipulator
            if geom.manipulator is not None:

                # Get design variable names from the manipulator
                currNames = geom.manipulator.getVarNames()

                # Add names to the list
                names = names + currNames

        # Return list of names
        return names

    def setDesignVars(self, dvDict):

        '''
        This function will set new values to the design variables.
        IT WILL NOT UPDATE THE SURFACE COORDINATES. However, it will flag
        all point set as outdated. The user should call self.update(ptSetName)
        to get the updated set of points. Note that this is also done by the
        CFD solver when you call it.

        Any additional keys in the DV dictionary are simply ignored.

        Note: you can use self.getValues to get a baseline dictionary with the correct
        structure, so you can change the desired DVs.

        Ney Secco 2017-02
        '''

        # Only the root processor will update its design variables
        if self.myID == 0:

            # Print log
            print ''
            print 'Setting new values for design variables to the manager.'

        # Loop over the primary geometries to find design variables
        for geom in self.geoms.itervalues():

            # Check if the geometry object has an associated geometry manipulator
            if geom.manipulator is not None:

                # Loop over every entry of the new DV dictionary
                for key in dvDict:

                    # Assign new DVs to the manipulator. Remember that we assume that
                    # the manipulator will ignore keys in dvDict that are not defined
                    # as design variables.
                    geom.manipulator.setDesignVars(dvDict)

        # All procs should flag all point sets as outdated
        for ptSetName in self.points:
            self.updated[ptSetName] = False

        # Print log
        if self.myID == 0:
            print 'Done'
            print ''

    def getNDV(self):

        '''
        This function given the total number of design variables under the current manager.

        Ney Secco 2017-02
        '''

        # Initialize dictionary of design variables
        dvDict = {}

        # Initialize DV counter
        NDV = 0

        # Loop over the primary geometries to find design variables
        for geom in self.geoms.itervalues():

            # Check if the geometry object has an associated geometry manipulator
            if geom.manipulator is not None:

                # Get design variables from the manipulator
                curr_dvDict = geom.manipulator_getDVs()

                # Increment the DV counter
                NDV = NDV + len(curr_dvDict)

        # Return number of DVs
        return NDV

    def pointSetUpToDate(self, ptSetName):

        '''
        This function just returns the state of the given point set.
        If False, it means that the nodal coordinates were not updated
        after the DV change. Thus the user should call self.update(ptSetName)

        Ney Secco 2017-02
        '''

        return self.updated[ptSetName]

    def update(self, ptSetName=None, childDelta=True, config=None):

        '''
        This function will update all surface coordinates under ptSetName based on
        the current values of design variables. The user should call self.setDesignVars
        before calling this function.

        If ptSetName is None, pySurf will update the triangulated surfaces and collar meshes,
        but it will not return anything useful.

        Ney Secco 2017-02
        '''

        # Only the root proc will run the manager functions
        if self.myID == 0:

            # Print log
            print ''
            print 'Updating manager manipulators for point set:',ptSetName

        # Loop over all primary geometry objects to update their manipulators
        for geom in self.geoms.itervalues():

            # Check if the geometry object has an associated geometry manipulator
            if geom.manipulator is not None:

                # Update the manipulator. Remember that this will update the associated
                # surface mesh and triangulated mesh
                geom.manipulator_update(ptSetName)

        # Now we need to repeat the tasks associated with the collar mesh generation
        self.run_baseFunction()

        # Print log
        if self.myID == 0:
            print 'Done'
            print ''

        # Set default value for the output
        pts = None

        # Update point set if requested by the user.
        if ptSetName is not None:
            pts = self.getSurfacePoints(ptSetName)
            self.updated[ptSetName] = True

        # Return the new set of solver points.
        # If the user did not set any point set, then this will be the default
        # value of self.points
        return pts

    def totalSensitivityProd(self, xDvDot, ptSetName, comm=None, config=None):

        '''
        This method executes the forward AD to compute the derivatives of the
        surface mesh with respect to the design variables.

        INPUTS:

        xDvDot: dictionary -> Dictionary containing derivative seeds of the design variables.
        It should follow the same structure as dvDict. You can use self.getValues to get
        a baseline dictionary and them change seeds.

        ptSetName: string -> Name of the point set that should be used to propagate derivatives

        OUTPUTS:

        xsDot: float[nPts,3] -> Derivative seeds of the surface mesh points. NOTE: This function
        will also update all derivative seeds stored throughout the manager.

        Ney Secco 2017-02
        '''

        # The root processor will do the main job
        if self.myID == 0:

            # Check if the current point set is updated
            if not self.updated[ptSetName]:
                raise NameError('The point set',ptSetName,'is outdated. Cannot propagate derivatives.')

        # First update the geometry manipulators. This will propagate derivative seeds from
        # design variables to all triangulated surface nodes, discrete curves, and primary structured
        # surface meshes associated with the geometry manipulators.
        # Note that all this derivative passing will be done directly to their corresponding objects.
        for geom in self.geoms.itervalues():
            if geom.manipulator is not None:
                geom.manipulator_forwardAD(xDvDot, ptSetName)

        # The root processor will do the main job
        if self.myID == 0:

            # Now we can propagate derivatives throughout the geometry operations done by the manager.
            # This will update derivative seeds of the surface collar meshes.
            self.forwardAD()

        # Now we need to gather the derivative seeds of all surface meshes, in solver ordering
        xsDot = self.getSurfaceForwardADSeeds(ptSetName)

        # Return derivative seeds
        return xsDot

    def totalSensitivity(self, xsBar, ptSetName, comm=None, config=None, clean=True):

        '''
        This method executes the reverse AD to compute the derivatives of the
        surface mesh with respect to the design variables.

        ATTENTION:
        This code will change the values in xsBar. So make sure you make a copy if you
        need the original values later on!

        INPUTS:

        xsBar: float[nPts,3] -> Dictionary containing derivative seeds of the design variables.
        It should follow the same structure as dvDict. You can use self.getValues to get
        a baseline dictionary and them change seeds.

        ptSetName: string -> Name of the point set that should be used to propagate derivatives

        OUTPUTS:

        xDvBar: dictionary -> Dictionary containing derivative seeds of the design variables.
        It will follow the same structure as dvDict.
        ATTENTION: If comm is provided, the derivatives should be reduced at the end.
                   That is, all procs should receive the sum of all seeds.

        Ney Secco 2017-02
        '''

        # Clean previous seeds
        self.clean_reverseADSeeds()

        # Initialize reverse seeds for design variables.
        # We will do this by getting a baseline dictionary with self.getValues and then
        # replacing all values with zeros.
        xDvBar = self.getValues()
        for key in xDvBar:
            xDvBar[key] = xDvBar[key]*0.0 # We do this to keep same array structure

        # Check if the current point set is updated
        if not self.updated[ptSetName]:
            raise NameError('The point set',ptSetName,'is outdated. Cannot propagate derivatives.')

        # Assign the derivative seeds to the surface meshes, assuming they are in solver ordering
        self.setSurfaceReverseADSeeds(xsBar, ptSetName)

        # Only the root proc will do the next seed propagation for its design variables
        if self.myID == 0:

            # Now we can propagate derivatives throughout the geometry operations done by the manager.
            # This will update derivative seeds of the triangulated surface meshes and discrete curves.
            self.reverseAD()

        # Then update the geometry manipulators. This will propagate derivative seeds from
        # all triangulated surface nodes, discrete curves, and primary structured
        # surface meshes to the design varibles associated with the geometry manipulators.
        # Note that all this derivative passing will be done directly to their corresponding objects.
        # Here we give None as communicator otherwise pyGeo would reduce the derivative seeds and then
        # we would do another allreduce here, yielding wrong results.
        for geom in self.geoms.itervalues():
            if geom.manipulator is not None:
                geom.manipulator_reverseAD(xDvBar, ptSetName, comm=None, clean=clean)

        # Do the reduction if we have a communicator object
        if comm is not None:
            for dvkey in xDvBar:
                xDvBar[dvkey] = self.comm.allreduce(xDvBar[dvkey], op=MPI.SUM)

        # Return design variable seeds
        return xDvBar

    def addVariablesPyOpt(self, optProb, globalVars=True, localVars=True,
                          sectionlocalVars=True, ignoreVars=None, freezeVars=None):
        """
        Add the current set of variables to the optProb object.

        Parameters
        ----------
        optProb : pyOpt_optimization class
            Optimization problem definition to which variables are added

        globalVars : bool
            Flag specifying whether global variables are to be added

        localVars : bool
            Flag specifying whether local variables are to be added

        ignoreVars : list of strings
            List of design variables the user DOESN'T want to use
            as optimization variables.

        freezeVars : list of string
            List of design variables the user WANTS to add as optimization
            variables, but to have the lower and upper bounds set at the current
            variable. This effectively eliminates the variable, but it the variable
            is still part of the optimization.
        """

        # The root processor will do the main job
        if self.myID == 0:

            print '#======================================'
            print 'Adding manipulator DVs to pyOpt'

            # Call the manipulators of each geometry object
            for geom in self.geoms.itervalues():
                if geom.manipulator is not None:
                    geom.manipulator.addVariablesPyOpt(optProb, globalVars, localVars,
                                                       sectionlocalVars, ignoreVars, freezeVars)
                    print 'Added variables from',geom.name

            print 'Done'
            print '#======================================'


    #=====================================================
    # DEBUG TOOLS

    def give_randomADSeeds_MACHinterface(self, ptSetName, fixedSeed=True):

        '''
        This method generates a set of random AD seeds to test the MACH interface functions above.
        This basically consists of a normalized set of forward AD seeds for the design variables,
        and a set of reverse AD seeds for the surface mesh points.
        This function only works after the manager object is initialized or reinitialized.
        '''

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        #======================
        # GENERATE FORWARD AD SEEDS FOR DESIGN VARIABLES

        # Copy DV dictionary to use as baseline
        xDvDot = self.getValues()

        # Assign normalized seeds to every entry
        for key in xDvDot:

            if isinstance(xDvDot[key],np.ndarray):
                xDvDot[key] = np.random.random_sample(xDvDot[key].shape)
                xDvDot[key] = xDvDot[key]/np.sqrt(np.sum(xDvDot[key]**2))
            else:
                xDvDot[key] = 1.0

        # The root proc will send its seeds to everyone
        xDvDot = self.comm.bcast(xDvDot, root=0)

        #======================
        # GENERATE REVERSE AD SEEDS FOR SURFACE MESH POINTS

        # Copy coordinate array, in solver ordering, as a baseline
        xsBar = self.getSurfacePoints(ptSetName)

        # Generate random seeds
        xsBar = np.random.random_sample(xsBar.shape)
        xsBar = xsBar/np.sqrt(np.sum(xsBar**2))

        #======================
        # RETURNS
        return xDvDot, xsBar

    #=====================================================
    # INTERNAL METHODS

    def _setSolverToManagerMapping(self, managerPts, solverPts, distTol=None):

        '''
        This method creates a mapping between the solver-provided surface mesh coordinates, and
        the surface mesh coordinates currently stored in the manager object.

        INPUTS:

        solverPts: array[numSolverPtsx3] -> Local array containing the coordinates given by the solver,
        in the current proc.

        managerPts: array[numManagerPtsx3] -> Local array containing the coordinates given by the manager,
        in the current proc. If you using this for collar meshes, remember to broadcast the coordinates
        from the root proc to all other procs.

        distTol: Distance tolerance to flag that a given surface node does not
                 belong to the current Manager surface definition.

        OUTPUTS:

        indexSolverPts: int[nMapping] -> Array containing indices of solverPts that are
                                         linked to managerPts nodes.

        indexManagerPts: int[nMapping] -> Array containing indices of solverPts that are
                                          linked to managerPts nodes.
                                          That is solverPts[indexSolverPts[ii]] = managerPts[indexManagerPts[ii]]

        numSurfRepetitions: int[nMapping] -> Number of times a given solverPts is repeated in managerPts.

        Ney Secco 2018-01
        '''

        # IMPORTS
        from scipy.spatial import cKDTree

        # Assign tolerance if user provided none
        if distTol is None:
            distTol = self.distTol

        # Get number of manager surface nodes
        numManagerPts = managerPts.shape[0]

        # Get number of solver surface nodes
        numSolverPts = solverPts.shape[0]

        # Initialize the tree with the manager nodes.
        # The manager might have repeated nodes (shared block edges for instance)
        tree = cKDTree(managerPts)

        # We need to be careful because some nodes in the manager vector may be repeated (shared by multiple blocks),
        # and we can't leave these repeated nodes out of the mapping. So we allow the KDTree to search for multiple
        # candidates, then we can assign all possible repeated mappings.
        # Here we set the maximum number we expect that a node may be repeated.
        maxRep = 6

        # Now use the KDTree to find which indices from managerPts is correlated to a given
        # node from solverPts.
        # That is solverPts[ii] = managerPts[indexMap[ii]]
        # If a given Solver node does not match with any node in the tree, then its
        # indexMap value will be numManagerPts (which is out of the range of managerPts)
        # We also allow the KDTree to search for the best k candidates, so we can
        # take care of repeated nodes.
        dist, indexMap = tree.query(solverPts, distance_upper_bound=distTol, k=maxRep)

        # Convert indexMap to a numpy array to facilitate next operations
        indexMap = np.array(indexMap)

        # At this point, indexMap is [numSolverPts x maxRep]. Therefore,
        # indexMap[i,j] gives the j-th candidate index in managerPts that corresponds
        # to the i-th node in solverPts.
        # So now we need to analyze column by column to take into account the multiple
        # candidates.

        # First let's initialize 1D arrays that we will increment for every candidate
        # analysis. Please see the comments over indexSolverPtsCurr and indexManagerPtsCurr to
        # understand the role of these arrays
        indexSolverPts = []
        indexManagerPts = []

        for candID in range(maxRep):

            # Now we need to remove entries from the indexMap that correspond to Solver nodes
            # that did not match with any Manager node in this proc. We can find these nodes
            # because indexMap[ii,candID]==numManagerPts.
            # Let's say we have numManagerPts = 5, and our current indexMap is
            # indexMap[:,candID] = [0 1 5 3 2 5 4]
            # This means that the 0th solverPts node is linked to the 0th managerPts node,
            # the 1st solverPts node is linked to the 1st managerPts node, the 2nd solverPts node
            # did not find any match in managerPts, the 3rd solverPts node is linked to the 3rd managerPts
            # node, the 4th solverPts node is linked to the 2nd managerPts node, and so on...
            # In the end, the important index relations are:
            #  solverPts -> managerPts
            #        0th -> 0th
            #        1th -> 1th
            #        3rd -> 3rd
            #        4th -> 2nd
            #        6th -> 4th
            # So we will create two arrays to express these relationships.
            # The first column will be indexSolverPtsCurr, and the second one will be indexManagerPtsCurr
            # These arrays will be concatenated into indexSolverPts and indexManagerPts to gather results
            # for all candidate orders.

            # Find indices of indexMap that are smaller than numManagerPts. These are
            # the indices that managed to find a match in the current proc, for the current candidate level.
            # Note that this also discard the candidates that are beyond the distance tolerance.
            indexSolverPtsCurr = np.where(indexMap[:,candID] < numManagerPts)[0]

            # Now take the corresponding links
            indexManagerPtsCurr = indexMap[indexSolverPtsCurr, candID]

            # Concatenate these arrays into the major arrays that gather results for all
            # candidates
            indexSolverPts = np.hstack([indexSolverPts, indexSolverPtsCurr])
            indexManagerPts = np.hstack([indexManagerPts, indexManagerPtsCurr])

        # Convert arrays to integers
        indexSolverPts = map(int, indexSolverPts)
        indexManagerPts = map(int, indexManagerPts)

        #---------------------------
        # Here is another very important detail.
        # We know that the Manager only repeats a surface node if it is shared by multiple blocks.
        # Therefore, the number of surface nodes in the Manager is exactly the same as the original
        # CGNS file, even if we are working with multiple procs. The Solver, on the other hand, may
        # duplicate surface nodes when working in parallel.
        # Our current mapping will map all repeated Solver nodes to all repeated Manager nodes,
        # regardless if they were generated by partitioning or shared block edges.
        # If we use our mapping to take coordinate values from the Manager vector (managerPts) and insert
        # them directly into the corresponding spot of the Solver vector (solverPts), we will
        # be fine since the repeated nodes will assign the same coordinate values.
        # However, if we are dealing with sensitivities, the same repeated nodes in the Manager vector
        # may have different sensitivities. Therefore, if we just insert values in the scatter operation
        # some sensitivity values may be lost, since the repeated nodes will always overwrite their values.
        # We can solve this by doing and additive scatter, and them taking the average of the added values!
        # If we define the coordinate of a Solver surface node as the average of all Manager nodes
        # that we linked with our mapping, then we will get the correct coordinate value since we will
        # take the average of repeated nodes. The nice thing is that this operation have a well-defined
        # differentiated version: just take the average of the derivative seeds!
        # So here we will count how many Manager nodes are linked to each Solver node, so we can
        # take the average of the additive scatters later on.

        # Initialize counter of repeated Manager surface nodes
        numSurfRepetitions = np.zeros(numSolverPts)

        # Now loop over the indices that Solver will map to in order to count the number of repetitions
        for solverIndex in indexSolverPts:

            # Note that solverIndex is the index that the current Manager node is linked to.
            # So we can increment the repetition counter.
            numSurfRepetitions[solverIndex] = numSurfRepetitions[solverIndex] + 1

        # Save the number of repetitions for future use
        # We need to make it a 2D array for matrix operations
        numSurfRepetitions = np.array([numSurfRepetitions]).T

        '''
        # Check if some nodes were not mapped
        unmatchedNodes = np.where(numSurfRepetitions == 0)[0]
        if len(unmatchedNodes > 0):
            print 'Some solver nodes did not find the corresponding copy in the manager object'
            print 'Here is the list of solver nodes:'
            print solverPts[unmatchedNodes,:]
            raise ValueError('Nodes not found')
        '''

        # Return the mapping info
        return indexSolverPts, indexManagerPts, numSurfRepetitions

    def _convertManagerToSolver(self, solverPts, managerPts, indexSolverPts, indexManagerPts, numSurfRepetitions):

        '''
        This uses the mapping to update surface nodes in the solver vector with the values from
        the manager vector.
        You need to run self._setSolverToManagerMapping() first to establish the mapping.

        INPUTS:

        solverPts: array[numSolverPtsx3] -> Local array containing the coordinates given by the solver,
        in the current proc.

        managerPts: array[numManagerPtsx3] -> Local array containing the coordinates given by the manager,
        in the current proc. If you using this for collar meshes, remember to broadcast the coordinates
        from the root proc to all other procs.

        indexSolverPts: int[nMapping] -> Array containing indices of solverPts that are
                                         linked to managerPts nodes.

        indexManagerPts: int[nMapping] -> Array containing indices of solverPts that are
                                          linked to managerPts nodes.
                                          That is solverPts[indexSolverPts[ii]] = managerPts[indexManagerPts[ii]]

        numSurfRepetitions: int[nMapping] -> Number of times a given solverPts is repeated in managerPts.

        OUTPUTS:

        This function modifies solverPts.

        Ney Secco 2018-01
        '''

        # First we remove previous data from the solver points that will be overwritten
        solverPts[indexSolverPts,:] = 0.0

        # Loop over every solver node to get its contribution from the manager nodes.
        for (solverPtID, managerPtID) in zip(indexSolverPts, indexManagerPts):

            # Assign value
            solverPts[solverPtID,:] = solverPts[solverPtID,:] + managerPts[managerPtID,:]/numSurfRepetitions[solverPtID]


    def _convertSolverToManagerb(self, solverPtsb, managerPtsb, indexSolverPts, indexManagerPts, numSurfRepetitions):

        '''
        This uses the mapping to update surface nodes in the solver vector with the values from
        the manager vector.
        You need to run self._setSolverToManagerMapping() first to establish the mapping.

        INPUTS:

        solverPtsb: array[numSolverPtsx3] -> Local array containing the derivative seeds given by the solver,
        in the current proc.

        managerPtsb: array[numManagerPtsx3] -> Local array containing the derivative seeds given by the manager,
        in the current proc. If you using this for collar meshes, remember to broadcast the coordinates
        from the root proc to all other procs.

        indexSolverPts: int[nMapping] -> Array containing indices of solverPts that are
                                         linked to managerPts nodes.

        indexManagerPts: int[nMapping] -> Array containing indices of solverPts that are
                                          linked to managerPts nodes.
                                          That is solverPts[indexSolverPts[ii]] = managerPts[indexManagerPts[ii]]

        numSurfRepetitions: int[nMapping] -> Number of times a given solverPts is repeated in managerPts.

        OUTPUTS:

        This function modifies managerPtsb and deletes the used seeds from solverPtsb.

        Ney Secco 2018-01
        '''

        # First we remove previous data from the solver points that will be overwritten
        managerPtsb[indexManagerPts,:] = 0.0

        # Loop over every solver node to get its contribution from the manager nodes.
        for (solverPtID, managerPtID) in zip(indexSolverPts, indexManagerPts):

            # Assign value
            managerPtsb[managerPtID,:] = managerPtsb[managerPtID,:] + solverPtsb[solverPtID,:]/numSurfRepetitions[solverPtID]

        # Burn the derivatives so that it won't be used again
        for solverPtID in indexSolverPts:
            solverPtsb[solverPtID,:] = 0.0

#=================================================
# AUXILIARY FUNCTIONS

def generate_primary_surface_mesh_filename(directory, geomName, fileNameTag=None):

    '''
    This just generates a filename for the surface mesh
    '''

    if fileNameTag is None:
        fileName = directory + geomName + '.xyz'
    else:
        fileName = directory + geomName + '_' + fileNameTag + '.xyz'

    return fileName

def generate_primary_volume_mesh_filename(directory, geomName, fileNameTag=None):

    '''
    This just generates a filename for the surface mesh
    '''

    if fileNameTag is None:
        fileName = directory + geomName + '.cgns'
    else:
        fileName = directory + geomName + '_' + fileNameTag + '.cgns'

    return fileName

def generate_collar_surface_mesh_filename(directory, curveName, fileNameTag=None):

    '''
    This just generates a filename for the surface mesh
    '''

    if fileNameTag is None:
        fileName = directory + curveName + '.xyz'
    else:
        fileName = directory + curveName + '_' + fileNameTag + '.xyz'

    return fileName

def generate_collar_volume_mesh_filename(directory, curveName, fileNameTag=None):

    '''
    This just generates a filename for the surface mesh
    '''

    if fileNameTag is None:
        fileName = directory + curveName + '.cgns'
    else:
        fileName = directory + curveName + '_' + fileNameTag + '.cgns'

    return fileName

def generate_nearfield_filename(directory, fileNameTag=None):

    '''
    This just generates a filename for the combined near-field meshes
    '''

    if fileNameTag is None:
        fileName = directory + 'near_field_meshes.cgns'
    else:
        fileName = directory + 'near_field_meshes_' + fileNameTag + '.cgns'

    return fileName

def generate_background_filename(directory, fileNameTag=None):

    '''
    This just generates a filename for the background mesh
    '''

    if fileNameTag is None:
        fileName = directory + 'autobg.cgns'
    else:
        fileName = directory + 'autobg_' + fileNameTag + '.cgns'

    return fileName

def generate_combined_filename(directory,fileNameTag=None):

    # Generate the name of the combined file in all procs
    if fileNameTag is None:
        combinedFileName = directory + 'aeroInput.cgns'
    else:
        combinedFileName = directory + 'aeroInput_' + fileNameTag + '.cgns'

    return combinedFileName
