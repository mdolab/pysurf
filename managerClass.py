from __future__ import division
import numpy as np
from mpi4py import MPI
import pysurf
from collections import OrderedDict

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

    def __init__(self):

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

    def clear_all(self):

        '''
        This method will clear all intersections, meshes, and tasks of the current manager,
        so that it can be used once again from scratch.
        The geometry objects will remain.
        '''

        self.intCurves = {}
        self.meshGenerators = {}
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

        if curveName in self.intCurves.keys():

            newCurve = self.intCurves[curveName].remesh(**optionsDict)

            # Store information regarding the parent curve (the one that was remeshed to get the new curve)
            newCurve.extra_data['parentCurve'] = self.intCurves[curveName].name

            # Rename the new curve
            newCurveName = newCurve.name

            # Assign intersection history
            if inheritParentGeoms:
                newCurve.extra_data['parentGeoms'] = self.intCurves[curveName].extra_data['parentGeoms'][:]
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
        
        if curveName in self.intCurves.keys():

            # Call split function
            splitCurvesDict = self.intCurves[curveName].split(optionsDict, criteria)

            # Add new curves to the manager's dictionary
            for curve in splitCurvesDict.itervalues():

                # Assign parents if necessary
                if inheritParentGeoms:
                    curve.extra_data['parentGeoms'] = self.intCurves[curveName].extra_data['parentGeoms'][:]

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

        # Get the name of the first curve
        mainCurveName = curveNames[0]
        mainCurve = self.intCurves[mainCurveName]

        # Call the mesh function from the main curve
        mergedCurve = mainCurve.merge(self.intCurves, mergedCurveName, curveNames[1:])

        # Check if we need to inherit parent geometry surfaces
        if inheritParentGeoms:
            mergedCurve.extra_data['parentGeoms'] = mainCurve.extra_data['parentGeoms'][:]

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

    def export_meshes(self, directory):

        '''
        This function will export all structured surface meshes into
        plot3d files. The files will be separated by primary geometries and
        also by collar meshes.

        INPUTS:

        directory: string -> Directory where will place all mesh files.
        '''

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
                fileName = generate_primary_surface_mesh_filename(directory, geom.name, primaryID)

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
                fileName = generate_collar_surface_mesh_filename(directory, curve.name, collarID)

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

    def extrude_meshes(self, directory):

        '''
        This function will use pyHyp to extrude all surface meshes into
        volume meshes.

        INPUTS:

        directory: string -> Directory where will place all mesh files.
        '''

        # Import pyHyp
        from pyhyp import pyHyp

        # Export the surface meshes once again just to make sure we have
        # the correct files available
        self.export_meshes(directory)

        # Print log
        print ''
        print 'Extruding surface meshes'

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

                # Generate file names
                surfFileName = generate_primary_surface_mesh_filename(directory, geom.name, primaryID)
                volFileName = generate_primary_volume_mesh_filename(directory, geom.name, primaryID)

                # Get extrusion options
                extrusionOptions = geom.meshObj.extrusionOptions
                
                # Give correct file name
                extrusionOptions['inputFile'] = surfFileName
                extrusionOptions['fileType'] = 'plot3d'

                # Extrude mesh
                hyp = pyHyp(options=extrusionOptions)
                hyp.run()
                hyp.writeCGNS(volFileName)

                # Increment counter
                primaryID = primaryID + 1

                # Print log
                print 'Extruded primary mesh for',geom.name

        # Now we will export the collar meshes
        for curve in self.intCurves.itervalues():
            
            # Verify if this curve was used to generate collar meshes
            if curve.extra_data['childMeshes'] is not None:

                # Generate file name
                surfFileName = generate_collar_surface_mesh_filename(directory, curve.name, collarID)
                volFileName = generate_collar_volume_mesh_filename(directory, curve.name, collarID)

                # Get extrusion options from the first mesh object
                meshName = curve.extra_data['childMeshes'][0]
                meshObj = self.meshGenerators[meshName].meshObj
                extrusionOptions = meshObj.extrusionOptions

                # Give correct file name
                extrusionOptions['inputFile'] = surfFileName
                extrusionOptions['fileType'] = 'plot3d'

                # Extrude mesh
                hyp = pyHyp(options=extrusionOptions)
                hyp.run()
                hyp.writeCGNS(volFileName)

                # Increment counter
                collarID = collarID + 1

                # Print log
                print 'Extruded collar mesh for',curve.name

        # Print log
        print 'Exported all meshes!'
        print ''        

    #=====================================================
    # MACH INTERFACE METHODS

    def addPointsSet(coor, ptSetName, origConfig=True, **kwargs):

        '''
        This function will receive an array of coordinates, and then assign
        these coordinates to the corresponding FFD, under the set ptSetName.

        ADflow will call this function to provide all surface points of the
        structured meshes, so we need to assign the correct points to their
        corresponding FFD. For instance, the surface mesh points of the wing should
        be assigned to the wing object FFD.

        Ney Secco 2017-02
        '''

#=================================================
# AUXILIARY FUNCTIONS

def generate_primary_surface_mesh_filename(directory, geomName, primaryID):

    '''
    This just generates a filename for the surface mesh
    '''

    fileName = directory + 'primary_%03d'%primaryID + '.xyz'

    return fileName

def generate_primary_volume_mesh_filename(directory, geomName, primaryID):

    '''
    This just generates a filename for the surface mesh
    '''

    fileName = directory + 'primary_vol_%03d'%primaryID + '.cgns'

    return fileName

def generate_collar_surface_mesh_filename(directory, curveName, collarID):

    '''
    This just generates a filename for the surface mesh
    '''
    
    # Generate file name
    fileName = directory + 'collar_%03d'%collarID + '.xyz'

    return fileName

def generate_collar_volume_mesh_filename(directory, curveName, collarID):

    '''
    This just generates a filename for the surface mesh
    '''
    
    # Generate file name
    fileName = directory + 'collar_vol_%03d'%collarID + '.cgns'

    return fileName
