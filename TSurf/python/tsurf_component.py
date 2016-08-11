from __future__ import division
import os
import numpy as np
from mpi4py import MPI
from ...baseClasses import Component
import adtAPI, cgnsAPI, curveSearch
from tsurf_geometry import getCGNSsections, merge_surface_sections, \
                           initialize_surface, initialize_curves, \
                           update_surface, remove_unused_points, \
                           translate, scale, rotate
import copy

# ALL AUXILIARY FUNCTIONS AND CLASSES ARE DEFINED IN adt_geometry.py

class TSurfComponent(Component):

    def _initialize(self, *arg):

        '''
        This function initializes and TSurfComponent.
        It is called by the __init__ method of the parent Component
        class defined in classes.py.

        The expected arguments for the initialization function are:
        TSurfComponent(fileName, sectionsList, comm)

        REQUIRED INPUTS:
        fileName: string -> Name of the CGNS file that contains the
                  triangulated surface definition.

        OPTIONAL INPUTS:
        sectionsList: list of strings -> List of strings containing
                  the names of the sections in the CGNS file that
                  should be included in the current ADTComponent.
                  If nothing is provided, or if sectionList is None
                  or an empty list, then all sections will be included.

        comm: MPI communicator -> An MPI communicator, such as
              MPI.COMM_WORLD
        '''

        # The first input is the file name
        filename = arg[0]

        # Set default values in case we have no additional arguments
        selectedSections = None # We will update this later if necessary
        self.comm = MPI.COMM_WORLD # Communicate to all processors

        # Check the optional arguments and do the necessary changes
        for optarg in arg[1:]:

            if type(optarg) == MPI.Intracomm: # We have an MPI Communicator
                self.comm = optarg

            elif type(optarg) == list:
                selectedSections = optarg

            elif optarg in (None, []):
                print 'Reading all CGNS sections in ADTComponent assigment.'

        # Read CGNS file
        self.coor, sectionDict = getCGNSsections(filename, self.comm)
        self.name = os.path.splitext(os.path.basename(filename))[0]

        # Select all section names in case the user provided none
        if selectedSections is None:
            selectedSections = sectionDict.keys()
        else:
            self.name = self.name + "__" + "_".join(selectedSections)

        # Now we call an auxiliary function to merge selected surface sections in a single
        # connectivity array
        self.triaConn, self.quadsConn = merge_surface_sections(sectionDict, selectedSections)

        # Initialize curves
        initialize_curves(self, sectionDict, selectedSections)

        # Now we remove unused points
        remove_unused_points(self)

        # Create ADT for the current surface, now that coor, triaConn, and quadsConn are correct
        initialize_surface(self)

    def add_curve(self, name, curve):

        '''
        Adds a given curve instance to the self.Curves dictionary.
        '''

        self.Curves[name] = copy.deepcopy(curve)

    def update(self, coor):

        '''
        This function updates the nodal coordinates used by both surface and curve objects.
        '''

        # Update coordinates
        self.coor = coor

        # Update surface definition
        update_surface(self)


    def translate(self, x, y, z):
        translate(self, x, y, z)
        update_surface(self)

    def scale(self, factor):
        scale(self, factor)
        update_surface(self)

    def rotate(self, angle, axis):
        rotate(self, angle, axis)
        update_surface(self)

    def project_on_surface(self, xyz):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        OUTPUTS:
        xyzProj -> float[numPts,3] : Coordinates of the projected points

        normProj -> float[numPts,3] : Surface normal at projected points
        '''

        '''
        Explanation of reference values:


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
        '''

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        normProj = np.zeros((numPts,3))

        # Call projection function
        procID, elementType, elementID, uvw = adtAPI.adtapi.adtmindistancesearch(xyz.T, self.name,
                                                                                 dist2, xyzProj.T,
                                                                                 self.nodal_normals, normProj.T)

        # Return projections
        return xyzProj, normProj

    def project_on_curve(self, xyz, curveCandidates=None):

        '''
        This function will compute projections and surface Normals

        INPUTS:
        xyz -> float[numPts, 3] : Coordinates of the points that should be projected.

        OUTPUTS:
        xyzProj -> float[numPts,3] : Coordinates of the projected points

        tanProj -> float[numPts,3] : Curve tangent at projected points
        '''

        '''
        Explanation of reference values:

        componentsList -> list of strings : Names of the surfaces components on which we should
                                            look for projections. If nothing is provided, the
                                            code will use all available surfaces.

        xyzProj -> float[numPts, 3] : If the user already have previous projection candidates, they
                                      should be provided in this array. The code will only replace
                                      values if it finds a better candidade. If the user has no
                                      previous candidate, initialize all elements to zero and also
                                      set large values to dist2.

        tanProj -> float[numPts, 3] : Curve tangent computed at the user-provided projection candidates.
                                      If the user has no previous candidate, initialize all elements
                                      to zero and also set large values to dist2.

        dist2 -> float[numPts] : distance**2 of the user-provided projection candidates. The code
                                 will use this distance value to check for the best candidate and
                                 then replace values in xyzProj, and normProj acordingly. If no previous
                                 candidates are available, set all elements to a large number (1e10) so that
                                 all information in xyzProj and normProj is replaced.
        '''

        # Use all curves if None is provided by the user
        if curveCandidates is None:
            curveCandidates = self.Curves.keys()

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        tanProj = np.zeros((numPts,3))

        # Check if the candidates are actually defined
        curveKeys = self.Curves.keys()
        for curve in curveCandidates:
            if curve not in curveKeys:
                print 'ERROR: Curve',curve,'is not defined. Check the curve names in your CGNS file.'
                print '       Also check if you included this curve name when selecting CGNS sections.'
                quit()

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curveName in self.Curves:
            if curveName in curveCandidates:
                self.Curves[curveName].project(xyz, dist2, xyzProj, tanProj)

        # Return projections
        return xyzProj, tanProj

#=============================================================
