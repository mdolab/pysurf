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

    def __init__(self, comm=None):

        # Set up MPI communicator if the user provided None
        if comm is None:
            comm = MPI.COMM_WORLD

        # Assign communicator
        self.comm = comm

        # Save ID of the current proc
        self.myID = self.comm.Get_rank()

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

        # Initialize field that will hold mapping between Solver and Manager
        self.indexSolverPts = None

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

    def initialize(self, directory, backgroundMeshInfo=None):

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

        # Only the root proc will work here
        if self.myID == 0:

            # Run base function to generate surface meshes
            self.run_baseFunction()

        # Extrude meshes
        combinedFileName = self.extrude_meshes(directory, backgroundMeshInfo)

        # Return combined file name to use in ADflow
        return combinedFileName

    def reinitialize(self):

        '''
        This method will do the pre-optimization step. This includes:

        - run base function to generate surface collar meshes
        - extrude all meshes with pyHyp
        - generate CGNS files with isolated meshes for pyWarp inputs
        - generate combined CGNS file with all meshes for ADflow

        Remember to use this at the beginning of an optimization script.

        INPUTS/OUTPUTS:
        
        Please refer to self.extrude_meshes to verify the inputs/outputs description.

        Ney Secco 2017-02
        '''

        # Run base function to generate surface meshes
        self.run_baseFunction()

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

        for meshGen in self.meshGenerators.itervalues():
            meshGen.meshObj.clean_reverseADSeeds()

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

    def march_intCurve_surfaceMesh(self, curveName, options0={}, options1={}, meshName=None, extrusionOptions={}):

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

            # Save the extrusion options into the first mesh object
            meshGen0.meshObj.extrusionOptions = extrusionOptions

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

    def export_meshes(self, directory, fileNameTag=''):

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

            # First we will export the primary geometry meshes
            for geom in self.geoms.itervalues():

                # Check if the component has an associated surface mesh
                if geom.meshObj is not None:

                    # Generate file name
                    fileName = generate_primary_surface_mesh_filename(directory, geom.name, primaryID, fileNameTag)

                    # Export mesh
                    geom.meshObj.export_plot3d(fileName)

                    # Increment counter
                    primaryID = primaryID + 1

                    # Print log
                    print 'Exported primary mesh for',geom.name

            # Now we will export the collar meshes
            for curve in self.intCurves.itervalues():

                # Verify if this curve was used to generate collar meshes
                if curve.extra_data['childMeshes'] is not None:

                    # Generate file name
                    fileName = generate_collar_surface_mesh_filename(directory, curve.name, collarID, fileNameTag)

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

                    # Print log
                    print 'Exported collar mesh for',curve.name

            # Print log
            print 'Exported all meshes!'
            print ''

    def extrude_meshes(self, directory, backgroundMeshInfo=None, fileNameTag=''):

        '''
        This function will use pyHyp to extrude all surface meshes into
        volume meshes.
        It will also use cgns_utils to generate the background mesh (if the user
        provided None) and then merge all blocks in a single CGNS file.

        INPUTS:

        directory: string -> Directory where will place all mesh files.

        backgroundMeshInfo: list or dict -> List containing filenames of the background
        meshes that should be appended to the combined CGNS file. The user could also
        provide a dictionary of options to use cgns_utils simpleOCart to generate a new
        background mesh. The dictionary fields should be:
           dh          Uniform cartesian spacing size
           hExtra      Extension in "O" dimension
           nExtra      Number of nodes to use for extension
           sym         Normal for possible sym plane
           mgcycle     Minimum MG cycle to enforce
        If the user provides None, then nothing will be added to the combined CGNS file.

        fileNameTag: string -> Optional tag to append to the file names.

        OUTPUTS:

        combinedFileName: string -> Filename of the CGNS file that contains all blocks
        (including the background mesh). This file should be used as input to ADflow.
        The blocks in the CGNS file always follow the ordering:
        primary meshes, collar meshes, background meshes.
        This order is important to match the coordinates between ADflow and pySurf.

        volFileList: list of strings -> List with names of CGNS files that contains the
        volume mesh of each geometry component and collar mesh. These files should be given
        to pyWarpMulti to initialize multiple instances corresponding to each mesh group.

        Ney Secco 2017-02
        '''

        # Import pyHyp
        from pyhyp import pyHyp

        # Export the surface meshes once again just to make sure we have
        # the correct files available
        self.export_meshes(directory, fileNameTag)

        # Print log
        print ''
        print 'Extruding surface meshes'

        # Add a slash if the directory does not have it
        if directory[-1] != '/':
            directory = directory + '/'

        # Initialize counters
        primaryID = 0
        collarID = 0

        # Initialize list of volume meshes
        volFileList = []

        # First we will export the primary geometry meshes
        for geom in self.geoms.itervalues():
            
            # Check if the component has an associated surface mesh
            if geom.meshObj is not None:

                # Generate file names
                surfFileName = generate_primary_surface_mesh_filename(directory, geom.name, primaryID, fileNameTag)
                volFileName = generate_primary_volume_mesh_filename(directory, geom.name, primaryID, fileNameTag)

                # Get extrusion options
                extrusionOptions = geom.meshObj.extrusionOptions
                
                # Give correct file name
                extrusionOptions['inputFile'] = surfFileName
                extrusionOptions['fileType'] = 'plot3d'

                # Extrude mesh
                hyp = pyHyp(options=extrusionOptions, comm=self.commSingle)
                hyp.run()
                hyp.writeCGNS(volFileName)

                # Increment counter
                primaryID = primaryID + 1

                # Append filename to the list
                volFileList.append(volFileName)

                # Print log
                print 'Extruded primary mesh for',geom.name

        # Now we will export the collar meshes
        for curve in self.intCurves.itervalues():
            
            # Verify if this curve was used to generate collar meshes
            if curve.extra_data['childMeshes'] is not None:

                # Generate file name
                surfFileName = generate_collar_surface_mesh_filename(directory, curve.name, collarID, fileNameTag)
                volFileName = generate_collar_volume_mesh_filename(directory, curve.name, collarID, fileNameTag)

                # Get extrusion options from the first mesh object
                meshName = curve.extra_data['childMeshes'][0]
                meshObj = self.meshGenerators[meshName].meshObj
                extrusionOptions = meshObj.extrusionOptions

                # Give correct file name (this one already contains all blocks of the collar mesh
                # since they were joined in the export_plot3d function).
                extrusionOptions['inputFile'] = surfFileName
                extrusionOptions['fileType'] = 'plot3d'

                # Extrude mesh
                hyp = pyHyp(options=extrusionOptions, comm=self.commSingle)
                hyp.run()
                hyp.writeCGNS(volFileName)

                # Increment counter
                collarID = collarID + 1

                # Append filename to the list
                volFileList.append(volFileName)

                # Print log
                print 'Extruded collar mesh for',curve.name

        # Print log
        print 'Extruded all meshes!'
        print ''

        ### ADDING BACKGROUND MESHES

        # Only the root proc will work here
        if self.myID == 0:

            ### First let's combine the near-field meshes in a single file

            # Initialize list of nearfield meshes
            nearfieldMeshes = []
            
            # Load files with CGNS utils
            for fileName in volFileList:
                nearfieldMeshes.append(cs.readGrid(fileName))

            # Combine grids in a single one
            nearfieldCombo = cs.combineGrids(nearfieldMeshes)

            # Save it in a single file
            nearfieldFilename = directory + 'near_field_meshes.cgns'
            nearfieldCombo.writeToCGNS(nearfieldFilename)

        # Check for background meshes
        if backgroundMeshInfo is None:

            # The user did not provide any background mesh
            bgMeshNames = ''

        elif isinstance(backgroundMeshInfo, dict):

            # Print log
            print ''
            print 'Generating background mesh with cgns_utils'

            # The user provided a set of options to run cgns_utils simpleOCart
            # to generate a new background mesh

            # Assign the background mesh filename for the next operations
            bgMeshNames = directory + 'Auto_background.cgns'

            # Define default set of options
            bg_options = {
                'dh':0.1,
                'hExtra':5.0,
                'nExtra':9,
                'sym':'z',
                'mgcycle':3,
            }

            # Replace default options by the user-provided ones
            for key in backgroundMeshInfo:
                if key not in bg_options.keys():
                    raise NameError('key',key,'not recognized by the background mesh generator.')
                else:
                    bg_options[key] = backgroundMeshInfo[key]

            # Run cartesian mesh generator from cgns_utils
            nearfieldFilename = directory + 'near_field_meshes.cgns'
            nearfieldGrid = cs.readGrid(nearfieldFilename)
            nearfieldGrid.simpleOCart(bg_options['dh'],
                                      bg_options['hExtra'],
                                      bg_options['nExtra'],
                                      bg_options['sym'],
                                      bg_options['mgcycle'],
                                      bgMeshNames,
                                      comm=self.comm)

            # Print log
            print 'Background mesh generated and saved as:'
            print bgMeshNames
            print ''

            # Transform file name to list so we could combine all meshes later on
            bgMeshNames = [bgMeshNames]

        # Now we can run the cgns_utils command to join all blocks in a single file

        # First we generate the name of the combined file in all procs
        combinedFileName = directory + 'aeroInput.cgns'

        # Only the root proc will work here
        if self.myID == 0:

            ### First let's combine the near-field meshes in a single file

            # Initialize list of all meshes
            allMeshes = []
            
            # We can already use all nearfield meshes
            allMeshes = allMeshes + nearfieldMeshes

            # Add background meshes 
            for fileName in bgMeshNames:
                allMeshes.append(cs.readGrid(fileName))

            # Combine grids in a single one
            allMeshCombo = cs.combineGrids(allMeshes)

            # Save it in a single file
            allMeshCombo.writeToCGNS(combinedFileName)

            '''
            combinedFileName = directory + 'aeroInput.cgns'
            os.system('cgns_utils combine '+ ' '.join(volFileList) + ' ' + bgMeshNames + ' ' + combinedFileName)
            '''

            # Check block-to-block connectivities
            os.system('cgns_utils connect ' + combinedFileName)

            # Print log
            print 'Combined all meshes! The combined CGNS file is:'
            print combinedFileName
            print 'This one should be used as input to ADflow.'
            print ''

        # Return the name of the combined file
        return combinedFileName

    #=====================================================
    # GENERAL INTERFACE METHODS
        
    def getSurfacePoints(self, solverOrder=True):

        '''
        This function returns the surface mesh points of all meshes contained
        in the manager object (including both primary and collar meshes).

        ATTENTION: These points are in Manager ordering.
        If you want to convert them to Solver ordering, you have
        to do:
        pts = self._convertManagerToSolver(pts)

        if solverOrder is True, the the points will be returned in solver ordering.
        Otherwise, they will be in manager ordering.
        '''

        # Only the root processor will get the points
        if self.myID == 0:

            # Initialize list of points with a single entry ([0, 0, 0])
            # We will remove this point at the end
            pts = np.zeros((1,3))

            # Loop over all primary meshes to gather their coordinates
            for geom in self.geoms.itervalues():

                # Get points of the current mesh
                currPts = geom.meshObj.get_points()

                # Append them to the list
                pts = np.vstack([pts, currPts])

            # Loop over all collar meshes to gather their coordinates
            for curve in self.intCurves.itervalues():

                # Verify if this curve was used to generate collar meshes
                if curve.extra_data['childMeshes'] is not None:

                    # Loop over every block of the collar mesh
                    for meshName in curve.extra_data['childMeshes']:

                        # Get points of the current mesh
                        currPts = self.meshGenerators[meshName].meshObj.get_points()

                        # Append them to the list
                        pts = np.vstack([pts, currPts])

            # Remove the initial dummy point
            pts = pts[1:,:]

        else:
            
            # Assign dummy variable in other procs
            pts = None

        # Convert the ordering if necessary
        if solverOrder:
            pts = self._convertManagerToSolver(pts)

        # Return the set of points
        return pts

    def getSurfaceForwardADSeeds(self, solverOrder=True):

        '''
        This function returns the forward AD seeds of all meshes contained
        in the manager object (including both primary and collar meshes).

        ATTENTION: These points are in Manager ordering.
        If you want to convert them to Solver ordering, you have
        to do:
        pts = self._convertManagerToSolver(pts)
        '''

        # Only the root processor will get the seeds
        if self.myID == 0:

            # Initialize list of points with a single entry ([0, 0, 0])
            # We will remove this point at the end
            ptsd = np.zeros((1,3))

            # Loop over all primary meshes to gather their coordinates
            for geom in self.geoms.itervalues():

                # Get points of the current mesh
                currPtsd = geom.meshObj.get_forwardADSeeds()

                # Append them to the list
                ptsd = np.vstack([ptsd, currPtsd])

            # Loop over all collar meshes to gather their coordinates
            for curve in self.intCurves.itervalues():

                # Verify if this curve was used to generate collar meshes
                if curve.extra_data['childMeshes'] is not None:

                    # Loop over every block of the collar mesh
                    for meshName in curve.extra_data['childMeshes']:

                        # Get points of the current mesh
                        currPtsd = self.meshGenerators[meshName].meshObj.get_forwardADSeeds()

                        # Append them to the list
                        ptsd = np.vstack([ptsd, currPtsd])

            # Remove the initial dummy point
            ptsd = ptsd[1:,:]

        else:
            
            # Assign dummy variable in other procs
            ptsd = None

        # Convert the ordering if necessary
        if solverOrder:
            ptsd = self._convertManagerToSolver(ptsd)

        # Return the set of points
        return ptsd

    def setSurfaceReverseADSeeds(self, ptsb, solverOrder=True):

        '''
        This function sets the reverse AD seeds of all meshes contained
        in the manager object (including both primary and collar meshes).

        INPUTS:

        ptsb: float[nPts,3] -> Reverse derivative seeds of all surface points,
        following the solver or manager ordering, depending on solverOrder flag.

        solverOrder: logical -> If solverOrder==True, we treat ptsb as if in solver ordering.
        Otherwise, we assume ptsb is in manager ordering.
        If solverOrder==False, only the root proc values will be useful.
        The values on the other procs will be discarded.
        '''

        # Convert ordering if necessary
        if solverOrder:
            ptsb = self._convertSolverToManagerb(ptsb)

        # Only the root processor will assign seeds
        if self.myID == 0:

            # Initialize offset variable to help us slice ptsb
            offset = 0

            # Loop over all primary meshes to set their seeds
            for geom in self.geoms.itervalues():

                # Get number of points in the current mesh
                numPts = geom.meshObj.numPts

                # Slice the global derivative seed vector
                curr_ptsb = ptsb[offset:offset+numPts]

                # Assign the derivative seeds
                geom.meshObj.set_reverseADSeeds(curr_ptsb)

                # Update the offset variable
                offset = offset + numPts

            # Loop over all collar meshes to set their seeds
            for curve in self.intCurves.itervalues():

                # Verify if this curve was used to generate collar meshes
                if curve.extra_data['childMeshes'] is not None:

                    # Loop over every block of the collar mesh
                    for meshName in curve.extra_data['childMeshes']:

                        # Get number of points in the current mesh
                        numPts = self.meshGenerators[meshName].meshObj.numPts

                        # Slice the global derivative seed vector
                        curr_ptsb = ptsb[offset:offset+numPts]

                        # Assign the derivative seeds
                        self.meshGenerators[meshName].meshObj.set_reverseADSeeds(curr_ptsb)

                        # Update the offset variable
                        offset = offset + numPts

    #=====================================================
    # MACH INTERFACE METHODS

    def addPointSet(self, coor, ptSetName, origConfig=True, distTol=1e-6, **kwargs):

        '''
        This function will receive an array of coordinates, and then assign
        these coordinates to the corresponding FFD, under the set ptSetName.

        ADflow will call this function to provide all surface points of the
        structured meshes, so we need to assign the correct points to their
        corresponding FFD. For instance, the surface mesh points of the wing should
        be assigned to the wing object FFD.

        In this function we have the following assumptions:

        - coor is a Nx3 array containing all relevant surface nodes (walls of all
        primary geometries). We assume that the nodes in coor follows the same CGNS ordering
        used to generate the combined file in self.extrude_meshes. That is, we have all primary meshes,
        collar meshes, and then background meshes. If you use a file generated by self.extrude_meshes
        as an input to ADflow, then you should be fine. This function will just verify if the nodes in
        coor have the expected ordering compared to the manager's objects.

        Ney Secco 2017-02
        '''

        # Print log
        print ''
        print 'Adding point set',ptSetName,'to the manager.'

        # Store the coordinates given by ADflow in a separate dictionary
        if ptSetName in self.points.keys():
            raise NameError('The point set',ptSetName,'is already defined under this manager.')
        else:

            # Save the coordinates in solver ordering (they are basically the same ones assigned to the FFDs)
            self.points[ptSetName] = coor.copy()

            # Flag that they are up to date
            self.updated[ptSetName] = True

        # Now we can create the mapping between the solver ordering and the manager ordering
        # Note that this will populate self.managerPoints.
        self._setSolverToManagerMapping(coor, distTol)


        # Print log
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

        # Only the root processor will gather the design variables
        if self.myID == 0:

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

        else:

            # Assign dummy variable for the MPI process
            dvDict = None

        # The root processor will send the dictionary to all others
        dvDict = self.comm.bcast(dvDict, root=0)

        # Return DV dictionary
        return dvDict

    def setDesignVars(self, dvDict):

        '''
        This function will set new values to the design variables.
        IT WILL NOT UPDATE THE SURFACE COORDINATES. However, it will flag
        all point set as outdated. The user should call self.update(ptSetName)
        to get the updated set of points.

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

                    # Get design variables from the manipulator
                    curr_dvDict = geom.manipulator_getDVs()

                    # Loop over every entry of the new DV dictionary
                    for key in dvDict:

                        # Assign new DVs to the manipulator. Remember that we assume that
                        # the manipulator will ignore keys in dvDict that are not defined
                        # as design variables.
                        geom.manipulator.setDesignVars(dvDict)

            # Flag all point sets as outdated
            for ptSetName in self.points:
                self.updated[ptSetName] = False

            # Print log
            print 'Done'
            print ''

    def getNDV(self):

        '''
        This function given the total number of design variables under the current manager.

        Ney Secco 2017-02
        '''

        # Only the root processor will count the design variables
        if self.myID == 0:

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

        else:

            # Assign dummy variable for the MPI process
            NDV = None

        # The root will send the value to all procs
        NDV = self.comm.bcast(NDV, root=0)

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

        Ney Secco 2017-02
        '''

        # Only the root proc will run the manager functions
        if self.myID == 0:

            # Print log
            print ''
            print 'Updating manager surface meshes.'

            # Loop over all primary geometry objects to update their manipulators
            for geom in self.geoms.itervalues():

                # Check if the geometry object has an associated geometry manipulator
                if geom.manipulator is not None:

                    # Update the manipulator. Remember that this will update the associated
                    # surface mesh and triangulated mesh
                    geom.manipulator_update()

            # Now we need to repeat the tasks associated with the collar mesh generation
            self.run_baseFunction()

            # Print log
            print 'Done'
            print ''

        # Gather the updated coordinates, in solver ordering
        pts = self.getSurfacePoints()#.copy()

        # Update corresponding dictionary entry if requested
        if ptSetName is not None:
            self.points[ptSetName] = pts
            self.updated[ptSetName] = True

        # Return the new set of points
        return pts

    def totalSensitivityProd(self, xDvDot, ptSetName):

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
                geom.manipulator_forwardAD(xDvDot)

            # Now we can propagate derivatives throughout the geometry operations done by the manager.
            # This will update derivative seeds of the surface collar meshes.
            self.forwardAD()

        # Now we need to gather the derivative seeds of all surface meshes, in solver ordering
        xsDot = self.getSurfaceForwardADSeeds()

        # Return derivative seeds
        return xsDot

    def totalSensitivity(self,xsBar, ptSetName, comm=None, config=None, clean=True):

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

        Ney Secco 2017-02
        '''

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
        self.setSurfaceReverseADSeeds(xsBar)

        # Only the root proc will do the next seed propagation for its design variables
        if self.myID == 0:

            # Now we can propagate derivatives throughout the geometry operations done by the manager.
            # This will update derivative seeds of the triangulated surface meshes and discrete curves.
            self.reverseAD()

            # Then update the geometry manipulators. This will propagate derivative seeds from
            # all triangulated surface nodes, discrete curves, and primary structured
            # surface meshes to the design varibles associated with the geometry manipulators.
            # Note that all this derivative passing will be done directly to their corresponding objects.
            for geom in self.geoms.itervalues():
                geom.manipulator_reverseAD(xDvBar, clean)

        else:

            # Since nothing happens with the DVs in the other procs, they should remain with zeroes as seeds.
            pass

        # Return design variable seeds
        return xDvBar

    #=====================================================
    # DEBUG TOOLS

    def give_randomADSeeds_MACHinterface(self, fixedSeed=True):

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
        xsBar = self.getSurfacePoints()

        # Generate random seeds
        xsBar = np.random.random_sample(xsBar.shape)
        xsBar = xsBar/np.sqrt(np.sum(xsBar**2))

        #======================
        # RETURNS
        return xDvDot, xsBar

    #=====================================================
    # INTERNAL METHODS

    def _setSolverToManagerMapping(self, solverPts, distTol=1e-6):

        '''
        This method creates a mapping between the solver-provided surface mesh coordinates, and
        the surface mesh coordinates currently stored in the manager object.

        INPUTS:

        solverPts: array[Nx3] -> Local array containing the coordinates given by the solver,
        in the current proc.

        distTol: Distance tolerance to flag that a given surface node does not
                 belong to the current Manager surface definition.

        OUTPUTS:

        This function has no explicit outputs. It sets:
        self.indexSolverPts
        self.indexManagerPts
        self.numSurfRepetitions
        self.startSolverIndex
        self.endSolverIndex
        self.numSolverPts
        self.numManagerPts

        Ney Secco 2017-03
        '''

        # IMPORTS
        from scipy.spatial import cKDTree

        # Each proc should compute how many nodes it has
        numSolverPts = solverPts.shape[0]

        # Now we do an allGather operation so that every proc knows how many nodes each other proc has
        numSolverPtsGather = np.hstack(self.comm.allgather(numSolverPts))

        # Now we use cumulative sums to find out the slice of the global surface node vector that
        # belongs to the current proc
        self.startSolverIndex = 0 + np.sum(numSolverPtsGather[:self.myID])
        self.endSolverIndex = self.startSolverIndex + numSolverPts

        # Send all coordinates to the root processor
        solverPts = self.comm.gather(solverPts.copy(), root=0)

        # Now only the root processor will do the mapping
        if self.myID == 0:

            # Concatenate all coordinates, since the gather operation brings a list of arrays
            solverPts = np.vstack(solverPts)

            # We can skip this if we already have a mapping
            if self.indexSolverPts is None:

                # Get manager coordinates (only the root proc has pySurf coordinates, so this is fine)
                managerPts = self.getSurfacePoints(solverOrder=False)

                # Get number of manager surface nodes
                numManagerPts = managerPts.shape[0]

                # Get number of solver surface nodes
                numSolverPts = solverPts.shape[0]

                # Store the number of points in each ordering
                self.numSolverPts = numSolverPts
                self.numManagerPts = numManagerPts

                # Initialize the tree with the manager nodes.
                # The manager might have repeated nodes (shared block edges for instance)
                tree = cKDTree(managerPts)

                # We need to be careful because some nodes in the manager vector may be repeated (shared by multiple blocks),
                # and we can't leave these repeated nodes out of the mapping. So we allow the KDTree to search for multiple
                # candidates, then we can assign all possible repeated mappings.
                # Here we set the maximum number we expect that a node may be repeated.
                maxRep = 6

                # Now use the KDTree to find which index from managerPts is correlated to a given
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

                    # Now we need to remove entries from the indexMap that correspond to ADflow nodes
                    # that did not match with any pyWarp node in this proc. We can find these nodes
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
                    indexManagerPtsCurr = indexMap[indexSolverPtsCurr,candID]

                    # Concatenate these arrays into the major arrays that gather results for all
                    # candidates
                    indexSolverPts = np.hstack([indexSolverPts,indexSolverPtsCurr])
                    indexManagerPts = np.hstack([indexManagerPts,indexManagerPtsCurr])

                # We don't need the local copy of the global surface vector anymore.
                # We also don't need the unsorted index mapping as well.
                del solverPts
                del indexMap

                # Now we can save the mapping for future use
                self.indexSolverPts = indexSolverPts
                self.indexManagerPts = indexManagerPts

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
                # If we define that the coordinate of a Solver surface node is the average of all Manager nodes
                # that we linked with our mapping, then we will get the correct coordinate value since we will
                # take the average of repeated nodes. The nice thing is that this operation have a well-defined
                # differentiated version: just take the average of the derivative seeds!
                # So here we will count how many Manager nodes are linked to each Solver node, so we can
                # take the average of the additive scatters later on.

                # Initialize counter of repeated Manager surface nodes
                numSurfRepetitions = np.zeros(self.numSolverPts)

                # Now loop over the indices that Solver will map to in order to count the number of repetitions
                for solverIndex in indexSolverPts:

                    # Note that solverIndex is the index that the current Manager node is linked to.
                    # So we can increment the repetition counter.
                    numSurfRepetitions[solverIndex] = numSurfRepetitions[solverIndex] + 1

                # Save the number of repetitions for future use
                # We need to make it a 2D array for matrix operations
                self.numSurfRepetitions = np.array([numSurfRepetitions]).T

    def _convertManagerToSolver(self, managerPts):

        '''
        This uses the mapping to convert surface nodes in Manager ordering to Solver ordering.
        You need to run self._setSolverToManagerMapping() first to establish the mapping.
        
        INPUTS:

        managerPts: array[Nx3] -> Local array containing the coordinates given by the manager,
        in the current proc. Since we are running this in serial, only the managerPts in the
        root proc will be actually used.

        OUTPUTS:



        Ney Secco 2017-03
        '''

        # The root proc should generate all entries of the solver vector.
        # Each proc will slice its data from this general vector later on
        if self.myID == 0:

            # Check if we already generated the mapping
            if self.indexSolverPts is None:
                raise NameError('There is no mapping between pySurf and the Solver. Run manager.addPointsSet first')

            # Define an array to store coordinates in the solver ordering
            solverPts = np.zeros((self.numSolverPts,3))

            # Loop over every solver node to get its contribution from the manager nodes.
            for (solverPtID, managerPtID) in zip(self.indexSolverPts,self.indexManagerPts):

                # Assign value
                solverPts[solverPtID,:] = solverPts[solverPtID,:] + managerPts[managerPtID,:]

            # Now divide the values to compute averaged values
            solverPts = solverPts/self.numSurfRepetitions

        else:

            # Assign a dummy variable to receive data later on
            solverPts = None

        # Now do a broadcast operation to send the full surface vector (in solver ordering) to all procs
        solverPts = self.comm.bcast(solverPts, root=0)

        # Then each proc takes its slice from the full vector
        solverPts = solverPts[self.startSolverIndex:self.endSolverIndex,:]

        return solverPts

    def _convertSolverToManager(self, solverPts):

        '''
        This uses the mapping to convert surface nodes in Solver ordering to Manager ordering.
        You need to run self._setSolverToManagerMapping() first to establish the mapping.
        Remember that we will send the full value back to the manager nodes. This can't be used
        as a reverse AD of the forward scatter. Use _convertSolverToManagerb for this purpose.

        INPUTS:

        solverPts: array[Nx3] -> Local array containing the coordinates given by the solver,
        in the current proc.

        Ney Secco 2017-03
        '''

        # Send all coordinates to the root processor
        solverPts = self.comm.gather(solverPts.copy(), root=0)

        # Now only the root processor will do the conversion
        if self.myID == 0:

            # Concatenate all coordinates, since the gather operation brings a list of arrays
            solverPts = np.vstack(solverPts)

            # Define an array to store to coordinates in the manager ordering
            managerPts = np.zeros((self.numManagerPts,3))

            # Loop over every solver node to send its contribution to the manager nodes
            for (solverPtID, managerPtID) in zip(self.indexSolverPts,self.indexManagerPts):

                # Accumulate contribution
                managerPts[managerPtID,:] = solverPts[solverPtID,:]

            # Return the averaged points
            return managerPts
        
        else:

            return None

    def _convertSolverToManagerb(self, solverPtsb):

        '''
        This uses the mapping to send reverse derivative seeds from the Solver ordering
        to the Manager ordering.
        You need to run self._setSolverToManagerMapping() first to establish the mapping.
        This is the reverse AD mode of self._convertManagerToSolver().

        INPUTS:

        solverPts: array[Nx3] -> Local array containing the coordinate seeds given by the solver,
        in the current proc.

        OUTPUTS:



        Ney Secco 2017-03
        '''

        # Send all coordinate seeds to the root processor
        solverPtsb = self.comm.gather(solverPtsb.copy(), root=0)

        # Now only the root processor will do the mapping
        if self.myID == 0:

            # Concatenate all coordinates, since the gather operation brings a list of arrays
            solverPtsb = np.vstack(solverPtsb)

            # Define an array to store to coordinates in the manager ordering
            managerPtsb = np.zeros((self.numManagerPts,3))

            # First we need to divide the solver seeds to take into account the reverse
            # AD of the division
            solverPtsb = solverPtsb/self.numSurfRepetitions

            # Loop over every solver node to send its contribution to the manager nodes
            for (solverPtID, managerPtID) in zip(self.indexSolverPts,self.indexManagerPts):

                # Accumulate contribution
                managerPtsb[managerPtID,:] = managerPtsb[managerPtID,:] + solverPtsb[solverPtID,:]

            # Return the averaged points
            return managerPtsb
        
        else:

            return None

#=================================================
# AUXILIARY FUNCTIONS

def generate_primary_surface_mesh_filename(directory, geomName, primaryID, fileNameTag=''):

    '''
    This just generates a filename for the surface mesh
    '''

    fileName = directory + fileNameTag + 'Primary_%03d'%primaryID + '.xyz'

    return fileName

def generate_primary_volume_mesh_filename(directory, geomName, primaryID, fileNameTag=''):

    '''
    This just generates a filename for the surface mesh
    '''

    fileName = directory + fileNameTag + 'Primary_vol_%03d'%primaryID + '.cgns'

    return fileName

def generate_collar_surface_mesh_filename(directory, curveName, collarID, fileNameTag=''):

    '''
    This just generates a filename for the surface mesh
    '''
    
    # Generate file name
    fileName = directory + fileNameTag + 'Collar_%03d'%collarID + '.xyz'

    return fileName

def generate_collar_volume_mesh_filename(directory, curveName, collarID, fileNameTag=''):

    '''
    This just generates a filename for the surface mesh
    '''
    
    # Generate file name
    fileName = directory + fileNameTag + 'Collar_vol_%03d'%collarID + '.cgns'

    return fileName
