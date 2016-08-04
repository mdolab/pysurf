from __future__ import division
import numpy as np
from mpi4py import MPI
from ...classes import Component
import adtAPI, cgnsAPI, curveSearch
from adt_geometry import Surface, Curve, getCGNSsections

class ADTComponent(Component):

    def _initialize(self, comm=MPI.COMM_WORLD.py2f(), *arg):
        filename = arg[0]

        self.comm = comm

        sectionDict = getCGNSsections(filename, self.comm)
        try:
            sections = arg[1]
        except IndexError:
            sections = sectionDict.keys()

        self.Surfaces = {}
        self.Curves = {}

        for section in sectionDict.itervalues():

            if section.type == 'surface' and section.name in sections: # We have a surface section

                self.Surfaces[section.name] = section

            elif section.type == 'curve' and section.name in sections: # We have a curve section

                self.Curves[section.name] = section

    def project_on_surface(self, xyz, surfCandidates=None):

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
        '''

        # Use all surfaces if None is provided by the user
        if surfCandidates is None:
            surfCandidates = self.Surfaces.keys()

        # Initialize reference values (see explanation above)
        numPts = xyz.shape[0]
        dist2 = np.ones(numPts)*1e10
        xyzProj = np.zeros((numPts,3))
        normProj = np.zeros((numPts,3))

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for surface in self.Surfaces.itervalues():
            if surface.name in surfCandidates:
                surface.project(xyz, dist2, xyzProj, normProj)

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

        # Call inverse_evaluate for each component in the list, so that we can update
        # dist2, xyzProj, and normProj
        for curve in self.Curves.itervalues():
            if curve.name in curveCandidates:
                curve.project(xyz, dist2, xyzProj, tanProj)

        # Return projections
        return xyzProj, tanProj
