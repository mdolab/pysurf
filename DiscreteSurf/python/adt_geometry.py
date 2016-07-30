from __future__ import division
import numpy as np
import cgnsAPI, adtAPI
from mpi4py import MPI

class Geometry(object):

    def __init__(self, inp, comm=MPI.COMM_WORLD.py2f()):

        # Assign communicator
        self.comm = comm

        # CHECKING INPUTS

        if isinstance(inp, str):
            # Read CGNS file
            cgnsAPI.cgnsapi.readcgns(inp, self.comm)

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

            # Save current names
            self.surfNames = surfNames
            self.curveNames = curveNames

            # Initialize dictionary
            inputDict = {}

            # Add coordinates to the dictionary
            inputDict['coor'] = coor

            # Now we split the connectivity arrays according to the sections
            iSurf = 0
            for surf in surfNames:

                # Get indices to slice connectivity array.
                # Remember to shift indices as Python starts at 0
                iTriaStart = surfTriaPtr[iSurf]-1
                iTriaEnd = surfTriaPtr[iSurf+1]-1
                iQuadsStart = surfQuadsPtr[iSurf]-1
                iQuadsEnd = surfQuadsPtr[iSurf+1]-1

                # Assemble dictionary entry for this surface
                surfDict = {'triaConn':triaConn[:,iTriaStart:iTriaEnd],
                            'quadsConn':quadsConn[:,iQuadsStart:iQuadsEnd]}

                # Add this entry to the geometry dictionary
                inputDict[surf] = surfDict

                # Increment surface counter
                iSurf = iSurf + 1

            iCurve = 0
            for curve in curveNames:

                # Get indices to slice connectivity array.
                # Remember to shift indices as Python starts at 0
                iBarsStart = curveBarsPtr[iSurf]-1
                iBarsEnd = curveBarsPtr[iSurf+1]-1

                # Take the bar connectivity slice and reorder it
                sortedConn = np.array(FEsort(barsConn[:,iBarsStart:iBarsEnd].T.tolist())).T

                # Assemble dictionary entry for this curve
                curveDict = {'barsConn':sortedConn}

                # Add this entry to the geometry dictionary
                inputDict[curve] = curveDict

                # Increment curve counter
                iCurve = iCurve + 1

        elif isinstance(inp, dict):
            pass
            # TODO: We still need to check the input here

        else:
            raise SyntaxError('Incorrect input to Surface; must be a cgns file string or dictionary containing coordinate and connection info.')

        # INITIALIZING SURFACE AND CURVE OBJECTS
        
        self.components = {}

        for componentName in inputDict.iterkeys():
            
            if componentName == 'coor': # We will not do anything with coor
                pass

            elif 'barsConn' in inputDict[componentName].keys(): # We have a curve component
                pass

            elif 'triaConn' in inputDict[componentName].keys(): # We have a surface component

                # Initialize object
                # TODO: Slice coor to take just the necessary elements
                currSurf = Surface(componentName, inputDict['coor'],
                                   inputDict[componentName]['quadsConn'], inputDict[componentName]['triaConn'],
                                   self.comm)

                # Add surface to the list
                self.components[componentName] = currSurf

    def get_names(self):

        '''
        This function will give the names of all components in the geometry object
        '''

        for componentName in self.components.keys():
            print componentName

    def update_points(self, coor):

        '''
        This function will update nodal coordinates.
        The connectivities will stay the same
        '''
        
        # Call update for each component
        for component in self.components.itervalues():
            component.update_points(coor)

    def project_on_surface(self, xyz, xyzProj, normProj, dist2, componentsList=None):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        componentsList -> list of strings : Names of the surfaces components on which we should
                                            look for projections. If nothing is provided, the
                                            code will use all available surfaces.

        xyzProj -> float[numPts, 3] : If the user already have previous projection candidates, they
                                      should be provided in this array. The code will only replace
                                      values if it finds a better candidade. If the user has no
                                      previous candidate, initialize all elements to zero and also
                                      set large values to dist2.

        normProj -> float[numPts, 3] : Surface normals computed at the user-provided projection candidates.
                                       If the user has no previous candidate, initialize all elements
                                       to zero and also set large values to dist2.

        dist2 -> float[numPts] : distance**2 of the user-provided projection candidates. The code
                                 will use this distance value to check for the best candidate and
                                 then replace values in xyzProj, and normProj acordingly. If no previous
                                 candidates are available, set all elements to a large number (1e10) so that
                                 all information in xyzProj and normProj is replaced.

        OUTPUTS
        This function has no explicit outputs. It will update the variables xyzProj, normProj, and dist2 instead.

        '''

        # Initialize values if the user did not provide any previous projection candidates
        if xyzProj == None or normProj == None or dist2 == None:
            numPts = xyz.shape[0]
            dist2 = np.ones(numPts)*1e10
            xyzProj = np.zeros((numPts,3))
            normProj = np.zeros((numPts,3))

        # If the user do not provide any component, we will project into all of them
        if componentsList == None:
            componentsList = self.surfNames

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for componentName in componentsList:
            self.components[componentName.lower()].project(xyz, dist2, xyzProj, normProj)

        print xyzProj

#=================================================================
# COMPONENT CLASSES
#=================================================================

class Surface(object):

    def __init__(self, name, coor, quadsConn, triaConn, comm=MPI.COMM_WORLD.py2f()):

        self.comm = comm

        self.coor = coor
        self.quadsConn = quadsConn
        self.triaConn = triaConn

        self.BBox = np.zeros((3, 2))
        self.useBBox = False
        self.adtID = name

        self.update_points(self.coor)

    def update_points(self, coor):
        self.coor = coor
        self.normals = self.compute_nodal_normals()
        adtAPI.adtapi.adtbuildsurfaceadt(self.coor, self.triaConn, self.quadsConn, self.BBox, self.useBBox, self.comm, self.adtID)

    def compute_nodal_normals(self):
        nCoor = self.coor.shape[1]
        nTria = self.triaConn.shape[1]
        nQuads = self.quadsConn.shape[1]
        nodal_normals = np.zeros((3, nCoor))
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

        xyz = xyz.T
        n = xyz.shape[1]
        arrDonor = self.nodal_normals

        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(xyz, self.adtID, dist2, xyzProj.T, arrDonor, normProj.T)

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

def FEsort(barsConnIn, flip=False):

    '''
    This function can be used to sort connectivities coming from the CGNS file
    It assumes that we only have one curve. It will crash if we have
    multiple curves defined in the same FE set.

    barsConnIn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]
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

    # Now we will build a new connectivity starting from this end node
    while len(barsConn) > 0:

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

                # Jump out of the for loop
                break

            else:

                # Just increment the element counter
                elemCounter = elemCounter + 1

    # Just change ordering if requested
    if flip:
        # Flip elements
        newConn = newConn[::-1]
        # Flip nodes
        for ii in range(len(newConn)):
            newConn[ii] = newConn[ii][::-1]

    # Return the sorted array
    return newConn

#=============================================================

if __name__ == "__main__":

    geometry = Geometry('../examples/cubeAndCylinder/cubeAndCylinder.cgns', MPI.COMM_WORLD.py2f())

    geometry.get_names()

    pts = np.array([[.5, .5, 0.1]], order='F')

    numPts = pts.shape[0]
    dist2 = np.ones(numPts)*1e10
    xyzProj = np.zeros((numPts,3))
    normProj = np.zeros((numPts,3))

    geometry.project_on_surface(pts,xyzProj,normProj,dist2,['CYLINDER'])

    print xyzProj
    print normProj
