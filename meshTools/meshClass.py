from __future__ import division
import numpy as np
from mpi4py import MPI

class Mesh(object):

    def __init__(self, *arg, **kwargs):

        # Call the initialization method defined in each
        # child class
        self._initialize(*arg)

    def generate_mesh(self):
        '''
        This function will be replaced by a function with the 
        same name defined in a child class
        '''
        pass

    def get_mesh(self):
        '''
        This function returns the surface mesh coodinates

        OUTPUT:
        coor: float[3,nu,nv] -> X, Y, and Z coordinates of the surface mesh
        '''
        return self.coor

    def export(self, *arg, **kwargs):
        '''
        This function will be replaced by a function with the 
        same name defined in a child class
        '''
        pass
        
