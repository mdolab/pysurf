from __future__ import division
import numpy as np
from mpi4py import MPI

class Component(object):

    def __init__(self, *arg, **kwargs):

        # Call the initialization method defined in each
        # child class
        self._initialize(*arg)

    def get_names(self):

        '''
        This function will give the names of all components in the geometry object
        '''

        print 'SURFACES:'
        for surface in self.Surfaces:
            print surface

        print 'CURVES:'
        for curve in self.Curves:
            print curve

    def project_on_surface():
        pass

    def project_on_curve():
        pass
