from __future__ import division
import numpy as np
import cgnsAPI, adtAPI, curveSearch
from mpi4py import MPI

def getCGNScomponents(inp, comm=MPI.COMM_WORLD.py2f()):

    # Read CGNS file
    cgnsAPI.cgnsapi.readcgns(inp, comm)

    # Retrieve data from the CGNS file
    coor = cgnsAPI.cgnsapi.coor
    triaConn = cgnsAPI.cgnsapi.triaconn
    quadsConn = cgnsAPI.cgnsapi.quadsconn
    barsConn = cgnsAPI.cgnsapi.barsconn
    surfTriaPtr = cgnsAPI.cgnsapi.surftriaptr
    surfQuadsPtr = cgnsAPI.cgnsapi.surfquadsptr
    curveBarsPtr = cgnsAPI.cgnsapi.curvebarsptr
    surfNames = cgnsAPI.cgnsapi.surfnames
    curveNames = cgnsAPI.cgnsapi.curvenames

    # Format strings coming to Fortran into a Python list
    surfNames = formatStringArray(surfNames)
    curveNames = formatStringArray(curveNames)

    # Initialize dictionary
    componentsDict = {}

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

        # Initialize surface object
        currSurf = Surface(surf, coor,
                           currQuadsConn, currTriaConn,
                           comm)

        # Add this component to the dictionary
        componentsDict[surf] = currSurf

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

        # Initialize surface object
        currCurve = Curve(curve, coor,
                          sortedConn,
                          comm)

        # Add this entry to the geometry dictionary
        componentsDict[curve] = currCurve

        # Increment curve counter
        iCurve = iCurve + 1

    # Return components dictionary
    return componentsDict

#=================================================================
# COMPONENT CLASSES
#=================================================================

class Surface(object):

    def __init__(self, name, coor, quadsConn, triaConn, comm=MPI.COMM_WORLD.py2f()):

        self.comm = comm

        self.type = 'surface'

        self.coor = coor
        self.quadsConn = quadsConn
        self.triaConn = triaConn

        self.name = name

        self.update_points(self.coor)

    def update_points(self, coor):
        self.coor = coor
        self.compute_nodal_normals()

        BBox = np.zeros((3, 2))
        useBBox = False

        adtAPI.adtapi.adtbuildsurfaceadt(self.coor, self.triaConn, self.quadsConn, BBox, useBBox, self.comm, self.name)

    def compute_nodal_normals(self):
        nCoor = self.coor.shape[1]

        nTria = self.triaConn.shape[1]
        nQuads = self.quadsConn.shape[1]
        connect_count = np.zeros((nCoor), dtype='int')

        self.nodal_normals = \
            adtAPI.adtapi.computenodalnormals(self.coor, self.triaConn, self.quadsConn)

    def project(self, xyz, dist2, xyzProj, normProj):

        '''
        This function will take the points given in xyz and project them to the surface.

        INPUTS:
        xyz -> float[nPoints,3] : coordinates of the points that should be projected

        dist2 -> float[nPoints] : values of best distance**2 found so far. If the distance of the projected
                                  point is less than dist2, then we will take the projected point. Otherwise,
                                  the projection will be set to (-99999, -99999, -99999), meaning that we could
                                  not find a better projection at this surface. This allows us to find the best
                                  projection even with multiple surfaces.
                                  If you don't have previous values to dist2, just initialize all elements to
                                  a huge number (1e10).

        xyzProj -> float[nPoints,3] : coordinates of projected points found so far. These projections could be in other
                                      surfaces as well. This function will only replace projections whose dist2 are smaller
                                      than previous dist2. This allows us to use the same array while working with multiple surfaces
                                      If you don't have previous values, just initialize all elements to zero. Also
                                      remember to set dist2 to a huge number so that all values are replaced.

        This function has no explicit outputs. It will just update dist2, xyzProj, and normProj
        '''

        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(xyz.T, self.name,
                                                                                 dist2, xyzProj.T,
                                                                                 self.nodal_normals, normProj.T)

class Curve(object):

    def __init__(self, name, coor, barsConn, comm=MPI.COMM_WORLD.py2f()):

        self.name = name

        self.comm = comm
        self.coor = coor
        self.barsConn = barsConn

        self.type = 'curve'

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

#===================================
# AUXILIARY FUNCTIONS
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

if __name__ == "__main__":

    componentsDict = getCGNScomponents('../examples/cubeAndCylinder/cubeAndCylinder.cgns', MPI.COMM_WORLD.py2f())

    pts = np.array([[.6, .5, 0.1]], order='F')

    for component in componentsDict.itervalues():

        numPts = pts.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        normProj = np.zeros((numPts,3))

        print 'Projecting on:',component.name

        component.project(pts,dist2,xyzProj,normProj)

        print xyzProj
        print normProj
