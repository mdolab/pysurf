from __future__ import division
import numpy as np
import cgnsAPI, adtAPI, curveSearch
from mpi4py import MPI

#===========================================
# AUXILIARY READING FUNCTIONS
#===========================================

def getCGNSsections(inputFile, comm=MPI.COMM_WORLD):

    '''
    This function opens a CGNS file, reads its sections, and returns a
    dictionary with selected sections. This code also initializes ADT
    for each section.
    '''

    # Read CGNS file
    cgnsAPI.cgnsapi.readcgns(inputFile, comm.py2f())

    # Retrieve data from the CGNS file.
    # We need to do actual copies, otherwise data will be overwritten if we read another
    # CGNS file.
    coor = np.array(cgnsAPI.cgnsapi.coor)
    triaConn = np.array(cgnsAPI.cgnsapi.triaconn)
    quadsConn = np.array(cgnsAPI.cgnsapi.quadsconn)
    barsConn = np.array(cgnsAPI.cgnsapi.barsconn)
    surfTriaPtr = np.array(cgnsAPI.cgnsapi.surftriaptr)
    surfQuadsPtr = np.array(cgnsAPI.cgnsapi.surfquadsptr)
    curveBarsPtr = np.array(cgnsAPI.cgnsapi.curvebarsptr)
    surfNames = cgnsAPI.cgnsapi.surfnames.copy()
    curveNames = cgnsAPI.cgnsapi.curvenames.copy()

    # Now we deallocate variables on the Fortran side
    cgnsAPI.cgnsapi.releasememory()

    # Format strings coming to Fortran into a Python list
    surfNames = formatStringArray(surfNames)
    curveNames = formatStringArray(curveNames)

    # Initialize dictionary
    sectionDict = {}

    # Now we split the connectivity arrays according to the sections
    iSurf = 0
    for surf in surfNames:

        # Get indices to slice connectivity array.
        # Remember to shift indices as Python starts at 0
        iTriaStart = surfTriaPtr[iSurf]-1
        iTriaEnd = surfTriaPtr[iSurf+1]-1
        iQuadsStart = surfQuadsPtr[iSurf]-1
        iQuadsEnd = surfQuadsPtr[iSurf+1]-1

        # Slice connectivities
        currTriaConn = triaConn[:,iTriaStart:iTriaEnd]
        currQuadsConn = quadsConn[:,iQuadsStart:iQuadsEnd]

        # Initialize surface dictionary
        currSurf = {'triaConn':currTriaConn,
                    'quadsConn':currQuadsConn}

        # Add this section to the dictionary
        sectionDict[surf] = currSurf

        # Increment surface counter
        iSurf = iSurf + 1

    iCurve = 0
    for curve in curveNames:

        # Get indices to slice connectivity array.
        # Remember to shift indices as Python starts at 0
        iBarsStart = curveBarsPtr[iCurve]-1
        iBarsEnd = curveBarsPtr[iCurve+1]-1

        # Take the bar connectivity slice and reorder it
        sortedConn = FEsort(barsConn[:,iBarsStart:iBarsEnd].T.tolist())

        # Check for failures
        if (sortedConn is not None) and (len(sortedConn)==1):
            sortedConn = np.array(sortedConn[0])
        else:
            # Just use the original connectivity
            sortedConn = barsConn[:,iBarsStart:iBarsEnd]
            # Print message
            print ''
            print 'Curve','"'+curve+'"','could not be sorted. It might be composed by disconnect curves.'
            print ''

        # Initialize curve dictionary
        currCurve = {'barsConn':sortedConn}

        # Add this entry to the dictionary
        sectionDict[curve] = currCurve

        # Increment curve counter
        iCurve = iCurve + 1

    # Return sections dictionary
    return coor, sectionDict

#========================================


#=================================================================
# AUXILIARY SURFACE FUNCTIONS
#=================================================================

def merge_surface_sections(sectionDict, selectedSections):

    '''
    This function merges the connectivity data of all surface sections
    specified in selectedSections.

    INPUTS:
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                  keys are section names, while the values are also dictionaries with the
                  following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                  for curve sections.
    '''

    # Initialize global connectivities
    triaConn = None
    quadsConn = None

    # Loop over the selected surfaces to gather connectivities
    for sectionName in selectedSections:

        if sectionName not in sectionDict.keys():

            # User wants a section that is not defined. Print a warning message:
            print 'ERROR: Surface',sectionName,'is not defined. Check surface names in your CGNS file.'
            quit()

        elif 'triaConn' in sectionDict[sectionName].keys(): # Then we have a surface section

            # Assign connectivities
            if triaConn is None:
                # Start new connectivities if we have none
                triaConn = sectionDict[sectionName]['triaConn']
                quadsConn = sectionDict[sectionName]['quadsConn']
            else:
                # Append new connectivities
                triaConn = np.hstack([triaConn, sectionDict[sectionName]['triaConn']])
                quadsConn = np.hstack([quadsConn, sectionDict[sectionName]['quadsConn']])

        else:

            # The user provided a name that is not a surface section
            print sectionName,'is not a surface section.'

    return triaConn, quadsConn

#========================================

def initialize_surface(TSurfComponent):

    '''
    This function receives a TSurf component and initializes its surface based
    on the initial connectivity and coordinates
    '''

    # Assign an adtID for the current surface. This ID should be a
    # unique string, so we select a name based on the number of trees
    # defined so far
    TSurfComponent.adtID = 'tree%08d'%adtAPI.adtapi.adtgetnumberoftrees()
    print 'My ID is',TSurfComponent.adtID

    # Now call general function that sets ADT
    update_surface(TSurfComponent)

#=================================================================

def update_surface(TSurfComponent):

    '''
    This function receives a TSurf component and initializes its surface
    ADT based on the current connectivity and coordinates.
    '''

    # Deallocate previous tree
    adtAPI.adtapi.adtdeallocateadts(TSurfComponent.adtID)

    # Set bounding box for new tree
    BBox = np.zeros((3, 2))
    useBBox = False

    # Compute set of nodal normals by taking the average normal of all
    # elements surrounding the node. This allows the meshing algorithms,
    # for instance, to march in an average direction near kinks.
    TSurfComponent.nodal_normals = adtAPI.adtapi.computenodalnormals(TSurfComponent.coor,
                                                                   TSurfComponent.triaConn,
                                                                   TSurfComponent.quadsConn)

    # Create new tree (the tree itself is stored in Fortran level)
    adtAPI.adtapi.adtbuildsurfaceadt(TSurfComponent.coor,
                                     TSurfComponent.triaConn, TSurfComponent.quadsConn,
                                     BBox, useBBox,
                                     TSurfComponent.comm.py2f(),
                                     TSurfComponent.adtID)

#=================================================================

#=================================================================
# AUXILIARY CURVE CLASSES AND FUNCTIONS
#=================================================================

class Curve(object):

    def __init__(self, coor, barsConn):

        self.coor = coor
        self.barsConn = barsConn

    def extract_points(self):

        # Get the number of elements
        numElems = self.barsConn.shape[1]

        # Initialize coordinates matrix
        pts = np.zeros((numElems+1, 3))

        # Get coordinates
        for ii in range(numElems):
            pts[ii,:] = self.coor[:, self.barsConn[0,ii]-1]

        # Get the last point
        pts[-1,:] = self.coor[:, self.barsConn[-1,-1]-1]

        # Return coordiantes
        return pts

    def update_points(self, coor):
        self.coor = coor

    def flip(self):

        # Flip elements
        self.barsConn = self.barsConn[::-1]
        # Flip nodes
        for ii in range(len(self.barsConn)):
            self.barsConn[ii] = self.barsConn[ii][::-1]


    def project(self, xyz, dist2, xyzProj, tangents):

        '''
        This function will take the points given in xyz and project them to the surface.

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected

        dist2 -> float[nPoints] : values of best distance**2 found so far. If the distance of the projected
                                  point is less than dist2, then we will take the projected point. This allows us to find the best
                                  projection even with multiple surfaces.
                                  If you don't have previous values to dist2, just initialize all elements to
                                  a huge number (1e10).

        xyzProj -> float[nPoints,3] : coordinates of projected points found so far. These projections could be on other curves as well. This function will only replace projections whose dist2 are smaller
                                      than previous dist2. This allows us to use the same array while working with multiple curves.
                                      If you don't have previous values, just initialize all elements to zero. Also
                                      remember to set dist2 to a huge number so that all values are replaced.

        tangents -> float[nPoints,3] : tangent directions for the curve at the projected points.

        This function has no explicit outputs. It will just update dist2, xyzProj, and tangents
        '''

        # Call fortran code
        curveSearch.curveproj.mindistancecurve(xyz.T, self.coor, self.barsConn, xyzProj.T, tangents.T, dist2)


    def remesh(self, nNewNodes=None, method='linear', spacing='linear'):

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
                 ['linear']

        OUTPUTS:
        This method has no explicit outputs. It will update self.coor and self.barsConn instead.

        Ney Secco 2016-08
        '''

        # Get connectivities and coordinates of the current Curve object
        coor = self.coor
        barsConn = self.barsConn

        # Get the number of nodes and elements in the curve
        nElem = barsConn.shape[1]
        nNodes = nElem + 1

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

        # Use the original number of nodes if the user did not specify any
        if nNewNodes is None:
            nNewNodes = nNodes

        # The arcLength will be used as our parametric coordinate to sample the new nodes.

        # First we will initialize arrays for the new coordinates and arclengths
        newNodeCoor = np.zeros((3,nNewNodes), order='F')
        newArcLength = np.zeros(nNewNodes, order='F')

        # Now that we know the initial and final arcLength, we can redistribute the
        # parametric coordinates based on the used defined spacing criteria.
        # These statements should initially create parametric coordinates in the interval
        # [0.0, 1.0]. We will rescale it after the if statements.
        if spacing == 'linear':
            newArcLength = np.linspace(0.0, 1.0, nNewNodes)

        elif spacing == 'cosine':
            newArcLength = 0.5*(1.0 - np.cos(np.linspace(0.0, np.pi, nNewNodes)))

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
        self.coor = newNodeCoor

        # Check if the baseline curve is periodic
        if barsConn[0,0] == barsConn[1,-1]:
            periodic = True
        else:
            periodic = False

        # Generate new connectivity (the nodes are already in order so we just
        # need to assign an ordered set to barsConn).
        barsConn = np.zeros((2,nNewNodes-1), order='F')
        barsConn[0,:] = range(1,nNewNodes)
        barsConn[1,:] = range(2,nNewNodes+1)

        # We still need to keep periodicity if the original curve is periodic
        if periodic:
            barsConn[1,-1] = barsConn[0,0]

        self.barsConn = barsConn

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

#=================================================================

def initialize_curves(TSurfComponent, sectionDict, selectedSections):

    '''
    This function initializes all curves given in sectionDict that are
    shown in selectedSections.

    INPUTS:
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                 keys are section names, while the values are also dictionaries with the
                 following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                 for curve sections.

    OUTPUTS:
    This function has no explicit outputs. It assigns TSurfComponent.Curves

    Ney Secco 2016-08
    '''

    # CURVE OBJECTS
    # A DiscreteSurf component may have multiple curve objects

    # Initialize dictionary that will hold all curve objects
    curveObjDict = {}

    # Now we will initialize curve objects
    for sectionName in sectionDict:

        if ('barsConn' in sectionDict[sectionName].keys()) and \
           (sectionName in selectedSections): # Then we have a curve section that should be stored

            # Get data from current section
            barsConn = sectionDict[sectionName]['barsConn']

            # Create Curve object and append entry to the dictionary
            curveObjDict[sectionName] = Curve(TSurfComponent.coor, barsConn)

    # Assign curve objects to ADT component
    TSurfComponent.Curves = curveObjDict

#=================================================================

def extract_curves_from_surface(TSurfComponent, criteria='sharpness'):

    '''
    This function will extract features from the surface of TSurfComponent and
    return Curve objects.

    Ney Secco 2016-08
    '''

    # Get points and connectivities
    coor = TSurfComponent.coor
    triaConn = TSurfComponent.triaConn
    quadsConn = TSurfComponent.quadsConn

    # We need to identify the neighbors of each element, so we can compute edge information
    

#=================================================================

#=================================================================

#=================================================================

#===================================
# GENERAL AUXILIARY FUNCTIONS
#===================================

def formatStringArray(fortranArray):

    '''
    This functions receives a numpy array of characters coming from Fortran
    and returns a list of formatted data.
    '''

    # Get array size
    shape = fortranArray.shape

    # We need to tranpose the array considering the Fortran ordering
    fortranArray = fortranArray.reshape([shape[1], shape[0]], order='F')

    # Initialize new array
    pythonArray = []

    # Now we format and concatenate the elements
    for index in range(shape[0]):

        # format string
        newString = ''.join(fortranArray[:,index]).strip().lower()

        # Add formatted string to the list
        pythonArray.append(newString)

    # Return formatted array
    return pythonArray

#=============================================================

def FEsort(barsConn):

    '''
    This function can be used to sort connectivities coming from the CGNS file
    It assumes that we only have one curve. It will crash if we have
    multiple curves defined in the same FE set.

    barsConnIn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

    newConn will be set to None if the sorting algorithm fails.
    '''

    # We will solve this in an iterative process until the connectivities do not change
    # anymore

    # Initialize the newest connectivities as our current input.
    newConn = barsConn

    # Say that we still need to search
    keep_searching = True

    while keep_searching:

        # If nothing happens, we will get out of the loop
        keep_searching = False

        # Update "old" connectivities
        oldConn = newConn[:]

        # Initialize list of new connectivities with the first element of the old ones.
        # We will also pop this element from oldConn
        newConn = [oldConn.pop(0)]

        # We will keep sorting until all elements of oldConn are assigned and popped out

        while len(oldConn) > 0:

            # Pop another element from oldConn
            oldElement = oldConn.pop(0)

            # We do not know if this element could be linked beforehand
            linked_element = False

            # Now we check if we can fit this element into any new connectivity
            newElemCounter = 0
            for newElement in newConn:

                if oldElement[0] == newElement[0]:

                    # We will flip the old element and place it at the beginning of the
                    # new connectivity
                    newConn[newElemCounter] = oldElement[::-1] + newElement[1:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[0] == newElement[-1]:

                    # Just append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[0]:

                    # Place the old element at the beginning
                    newConn[newElemCounter] = oldElement + newElement[1:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[-1]:

                    # Flip and append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement[::-1]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                # Increment element counter
                newElemCounter = newElemCounter + 1

            if not linked_element: # We have an element that is not connected to anyone

                    # Define a new curve in newConn
                    newConn.append(oldElement)

    # Right now, newConn represent each line by an 1D array of point IDs.
    # e.g. [[2,3,4,5],[1,7,8]]
    # We need to convert this back to FE format, represented by 2D arrays.
    # e.g. [[[2,3,4],[3,4,5]],[[1,7],[7,8]]]

    # Initialize FE array for each curve
    newConnFE = [np.zeros((2,len(curve)-1),dtype=int) for curve in newConn]

    # Fill FE data for each curve
    for curveID in range(len(newConn)):

        # Get current curve
        curve = newConn[curveID]

        # Get FE list of the current curve
        FEcurve = newConnFE[curveID]

        # Assign FE connectivity
        for pointID in range(1,len(curve)):

            # Get point indices
            prevPoint = curve[pointID-1]
            currPoint = curve[pointID]

            # Assign bar FE
            FEcurve[0,pointID-1] = prevPoint
            FEcurve[1,pointID-1] = currPoint


    # Now we do a final check to remove degenerate bars (bars that begin and end at the same point)
    # Take every disconnect curve in newConnFE:
    for curveID in range(len(newConnFE)):

        # We convert the connectivity array to list so we can 'pop' elements
        curveFE = newConnFE[curveID].T.tolist()

        for FE in curveFE:

            # Check if the start and end points are the same
            if FE[0] == FE[1]:

                # Remove FE
                curveFE.remove(FE)

        # Convert connectivity back to numpy array
        newConnFE[curveID] = np.array(curveFE, order='F').T

    # Return the sorted array
    return newConnFE

#=============================================================

def FEsort_dumb(barsConnIn):

    '''
    This function can be used to sort connectivities coming from the CGNS file
    It assumes that we only have one curve. It will crash if we have
    multiple curves defined in the same FE set.

    barsConnIn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

    newConn will be set to None if the sorting algorithm fails.
    '''

    # Find if we have any node that is used only once.
    # This indicates that the curve has free ends.
    # First we need to flatten the array
    flatArray = [inner for outer in barsConnIn for inner in outer]

    # Initialize logical variable to stop iteration
    end_not_found = True

    while end_not_found and len(flatArray) > 0:

        # Remove one node from the array
        node_to_check = flatArray.pop(0)

        # Check for duplicates
        if node_to_check not in flatArray:

            # We found an end!
            end_node = node_to_check
            end_not_found = False

    # If we did not find an end node, the curve is periodic, so we can start in any of them
    if end_not_found:
        end_node = barsConnIn[0,0]

    # Initialize new connectivity list
    newConn = []

    # Create a copy of bars conn as we will pop things out
    barsConn = barsConnIn[:]

    # We should stop the code if it finds nodes connected to three elements or disconnected curves.
    # We will use the noFail flag for this task
    noFail = True

    # Now we will build a new connectivity starting from this end node
    while (len(barsConn) > 0) and noFail:

        # Assume we have a fail
        noFail = False

        # Initialize new element
        newElement = [end_node, 0]

        # Set counter for number of elements checked in barsConn
        elemCounter = 0

        # Loop over the elements in barsConn to find one the has the end node
        for element in barsConn:

            if end_node in element:

                # Assign the second nodes in the new element
                if end_node == element[0]:
                    newElement[1] = element[1]
                elif end_node == element[1]:
                    newElement[1] = element[0]

                # Update end_node
                end_node = newElement[1]

                # Append the new element to the list
                newConn.append(newElement)

                # Pop the old element from the old connectivities
                barsConn.pop(elemCounter)

                # If we assign an element during this loop, we can continue
                noFail = True

                # Jump out of the for loop
                break

            else:

                # Just increment the element counter
                elemCounter = elemCounter + 1

    # Check for fails
    if noFail is False:

        # Set new connectivity to None, indicating failure
        newConn = None

    # Return the sorted array
    return newConn

#=============================================================

def remove_unused_points(TSurfComponent):

    '''
    This function will removed unused points from coor and
    also update all connectivities.
    '''

    # Gather data
    coor = TSurfComponent.coor
    triaConn = TSurfComponent.triaConn
    quadsConn = TSurfComponent.quadsConn

    # Get total number of points and elements
    nPoints = coor.shape[1]
    nTria = triaConn.shape[1]
    nQuads = quadsConn.shape[1]

    # Initialize mask to identify used points
    usedPtsMask = np.zeros(nPoints, dtype=int)

    # First we loop over all elements to create a mask that indicates used points
    for triaID in range(nTria):
        # Flag used points
        usedPtsMask[triaConn[0,triaID]-1] = 1
        usedPtsMask[triaConn[1,triaID]-1] = 1
        usedPtsMask[triaConn[2,triaID]-1] = 1

    for quadID in range(nQuads):
        # Flag used points
        usedPtsMask[quadsConn[0,quadID]-1] = 1
        usedPtsMask[quadsConn[1,quadID]-1] = 1
        usedPtsMask[quadsConn[2,quadID]-1] = 1
        usedPtsMask[quadsConn[3,quadID]-1] = 1

    for curve in TSurfComponent.Curves.itervalues():
        # Get current connectivity
        barsConn = curve.barsConn

        # Get number of elements in this curve
        nBars = barsConn.shape[1]

        for barID in range(nBars):
            # Flag used points
            usedPtsMask[barsConn[0,barID]-1] = 1
            usedPtsMask[barsConn[1,barID]-1] = 1

    # Now we can compute the number of points actually used
    nUsedPts = np.sum(usedPtsMask)

    # Initialize new coordinate array
    cropCoor = np.zeros((3,nUsedPts),order='F')

    # Initialize counter to fill cropCoor
    cropPointID = -1

    # Now we fill the points of the cropped array
    for pointID in range(nPoints):

        # Check if the point is used
        if usedPtsMask[pointID] == 1:

            # Increment counter
            cropPointID = cropPointID + 1

            # Add point to the cropped array
            cropCoor[:,cropPointID] = coor[:,pointID]

            # Now we replace the value in the mask array so that we
            # can use it as a pointer from coor to cropCoor when we
            # update element connectivities.
            # The +1 is necessary because Fortran use 1-based indexing
            usedPtsMask[pointID] = cropPointID + 1

    # Store the new set of points
    TSurfComponent.coor = cropCoor

    # Now we need to update connectivities so that they point to the the correct
    # indices of the cropped array
    for triaID in range(nTria):
        # Use pointer to update connectivity
        triaConn[0,triaID] = usedPtsMask[triaConn[0,triaID]-1]
        triaConn[1,triaID] = usedPtsMask[triaConn[1,triaID]-1]
        triaConn[2,triaID] = usedPtsMask[triaConn[2,triaID]-1]

    for quadID in range(nQuads):
        # Use pointer to update connectivity
        quadsConn[0,quadID] = usedPtsMask[quadsConn[0,quadID]-1]
        quadsConn[1,quadID] = usedPtsMask[quadsConn[1,quadID]-1]
        quadsConn[2,quadID] = usedPtsMask[quadsConn[2,quadID]-1]
        quadsConn[3,quadID] = usedPtsMask[quadsConn[3,quadID]-1]

    for curve in TSurfComponent.Curves.itervalues():
        # Update coordinates
        curve.coor = cropCoor

        # Get current connectivity
        barsConn = curve.barsConn

        # Get number of elements in this curve
        nBars = barsConn.shape[1]

        for barID in range(nBars):
            # Flag used points
            barsConn[0,barID] = usedPtsMask[barsConn[0,barID]-1]
            barsConn[1,barID] = usedPtsMask[barsConn[1,barID]-1]
