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

                # Assemble dictionary entry for this curve
                curveDict = {'barsConn':barsConn[:,iBarsStart:iBarsEnd]}

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

    def inverse_evaluate(self, xyz, componentsList=None, xyzProj=None, normProj=None, dist2=None):

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
                                      previous candidate, he could use xyzProj = None.

        normProj -> float[numPts, 3] : Surface normals computed at the user-provided projection candidates.

        dist2 -> float[numPts] : distance**2 of the user-provided projection candidates. The code
                                 will use this distance value to check for the best candidate and
                                 then replace values in xyzProj, and normProj acordingly.
        '''

        # Initialize values if the user did not provide any previous projection candidates
        if dist2 == None:
            numPts = xyz.shape[0]
            dist2 = np.ones(numPts)*1e10
            xyzProj = np.zeros((numPts,3))
            vecProj = np.zeros((numPts,3))

        # If the user do not provide any component, we will project into all of them
        if componentsList == None:
            componentsList = self.surfNames

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and vecProj
        for componentName in componentsList:
            self.components[componentName.lower()].inverse_evaluate(xyz, dist2, xyzProj, vecProj)

        # For surfaces, projVectors are surface normals. For curves, projVectors are
        # curve tangents.

        # Return projections
        return xyzProj, vecProj

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

    def inverse_evaluate(self, xyz, dist2, xyzProj, normProj):

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

if __name__ == "__main__":

    geometry = Geometry('../examples/cubeAndCylinder/cubeAndCylinder.cgns', MPI.COMM_WORLD.py2f())

    geometry.get_names()

    pts = np.array([[.5, .5, 0.1]], order='F')
    projPoints, projVectors = geometry.inverse_evaluate(pts,['CYLINDER'])

    print projPoints
    print projVectors
