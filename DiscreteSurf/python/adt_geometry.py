from __future__ import division
import numpy as np
import cgnsAPI, adtAPI, curveSearch
from mpi4py import MPI

#===========================================
# AUXILIARY FUNCTIONS
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
    coor = cgnsAPI.cgnsapi.coor.copy()
    triaConn = cgnsAPI.cgnsapi.triaconn.copy()
    quadsConn = cgnsAPI.cgnsapi.quadsconn.copy()
    barsConn = cgnsAPI.cgnsapi.barsconn.copy()
    surfTriaPtr = cgnsAPI.cgnsapi.surftriaptr.copy()
    surfQuadsPtr = cgnsAPI.cgnsapi.surfquadsptr.copy()
    curveBarsPtr = cgnsAPI.cgnsapi.curvebarsptr.copy()
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
        if sortedConn is not None:
            sortedConn = np.array(sortedConn).T
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
    This function initializes a single ADT for all triangulated surfaces defined in selectedSections.

    INPUTS:
    section_dict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
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
    
    # Instantiate the Surface class. This step will initialize ADT
    # with the global connectivities
    # AAA surfObj = Surface(coor, quadsConn, triaConn, comm)

    return triaConn, quadsConn

#========================================

def initialize_surface(ADTComponent):

    '''
    This function receives an ADT component and initializes its surface based on the
    initial connectivity and coordinates
    '''

    # Assign an adtID for the current surface. This ID should be a unique string, so
    # we select a name based on the number of trees defined so far
    ADTComponent.adtID = 'tree%08d'%adtAPI.adtapi.adtgetnumberoftrees()
    print 'My ID is',ADTComponent.adtID

    # Now call general function that sets ADT
    update_surface(ADTComponent)

#=================================================================

def update_surface(ADTComponent):

    '''
    This function receives an ADT component and initializes its surface based on the
    current connectivity and coordinates. If deallocate==True, this code will try to
    deallocate the previous ADT. This should not be used when creating the initial ADT,
    as there is nothing to deallocate.
    '''

    # Deallocate previous tree
    adtAPI.adtapi.adtdeallocateadts(ADTComponent.adtID)

    # Set bounding box for new tree
    BBox = np.zeros((3, 2))
    useBBox = False

    # Compute set of nodal normals by taking the average normal of all elements surrounding
    # the node. This allows the meshing algorithms, for instance, to march in an average
    # direction near kinks.
    ADTComponent.nodal_normals = adtAPI.adtapi.computenodalnormals(ADTComponent.coor,
                                                                   ADTComponent.triaConn,
                                                                   ADTComponent.quadsConn)

    # Create new tree (the tree itself is stored in Fortran level)
    adtAPI.adtapi.adtbuildsurfaceadt(ADTComponent.coor, 
                                     ADTComponent.triaConn, ADTComponent.quadsConn,
                                     BBox, useBBox,
                                     ADTComponent.comm.py2f(),
                                     ADTComponent.adtID)

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

        curveSearch.curveproj.mindistancecurve(xyz.T, self.coor, self.barsConn, xyzProj.T, tangents.T, dist2)

#=================================================================

def initialize_curves(ADTComponent, sectionDict):

    '''
    This function initializes a single ADT for all triangulated surfaces defined in selectedSections.

    INPUTS:
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                 keys are section names, while the values are also dictionaries with the
                 following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                 for curve sections.
    
    OUTPUTS:
    This function has no explicit outputs. It assigns ADTComponent.Curves
    '''

    # CURVE OBJECTS
    # A DiscreteSurf component may have multiple curve objects

    # Initialize dictionary that will hold all curve objects
    curveObjDict = {}

    # Now we will initialize curve objects
    for sectionName in sectionDict:

        if 'barsConn' in sectionDict[sectionName].keys(): # Then we have a curve section

            # Get data from current section
            barsConn = sectionDict[sectionName]['barsConn']

            # Create Curve object and append entry to the dictionary
            curveObjDict[sectionName] = Curve(ADTComponent.coor, barsConn)

    # Assign curve objects to ADT component
    ADTComponent.Curves = curveObjDict

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

def FEsort(barsConnIn):

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
