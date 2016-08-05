from __future__ import division
import numpy as np
from mpi4py import MPI

class Component(object):

    def __init__(self, *arg, **kwargs):
        print arg
        self._initialize(*arg)

    def get_names(self):

        '''
        This function will give the names of all components in the geometry object
        '''

        print 'SURFACES:'
        for component in self.Surfaces:
            print component

        print 'CURVES:'
        for component in self.Curves:
            print component

    def update_points(self, componentName, coor):

        '''
        This function will update nodal coordinates.
        The connectivities will stay the same
        '''

        if componentName in self.Surfaces.keys():
            self.Surfaces[componentName].update_points(coor)

        if componentName in self.Curves.keys():
            self.Curves[componentName].update_points(coor)

    def project_on_surface():
        pass

    def project_on_curve():
        pass
