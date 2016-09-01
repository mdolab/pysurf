from __future__ import division
import numpy as np
from mpi4py import MPI

class Geometry(object):

    '''
    This is the base Geometry class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.

    The required attributes for a Geometry object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.

    self.curves : dictionary{curveName:curveObject} -> Dictionary with component curves.
    '''


    def __init__(self, *arg, **kwargs):
        '''
        Call the initialization method defined in each
        derived class
        '''

        self.name = ''
        self.comm = None
        self.curves = {}

        if 'comm' in kwargs:
            self.comm = kwargs['comm']
        else:
            self.comm = MPI.COMM_WORLD

        self._initialize(*arg, **kwargs)

    def _initialize(self, *arg):
        '''
        Virtual method
        '''
        pass

    def translate(self, xyz):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def scale(self, factor):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def rotate(self, angle, axis, point=None):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        INPUTS:
        angle: float -> Rotation angle, in degrees.

        axis: integer -> Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3] -> Coordinates of the rotation center. If point=None,
                           the origin should be used.
        '''
        pass

    def project_on_surface(self, xyz):
        '''
        This function projects a point in space to the component surface.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def project_on_curve(self, xyz, curveCandidates=None):
        '''
        This function projects a point in space to the component curves.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

class Curve(object):


    '''
    This is the base Curve class.
    All geometry engines should have their own
    derived Geometry classes that contains the
    same attributes and methods defined here.
    Remember that a Geometry object may hold Curve objects.

    The required attributes for a Curve object are:

    self.name : string -> Geometry name.

    self.comm : MPI communicator.
    '''

    def __init__(self, *arg, **kwargs):
        '''
        Call the initialization method defined in each
        child class
        '''
        self._initialize(*arg)

    def extract_points(self):
        '''
        This method should return coordinates along the curve.
        '''
        pass

    def update_dvs(self, coor):
        self.coor = coor

    def flip(self):
        '''
        This flip should flit the curve direction.
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def translate(self, xyz):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def scale(self, factor):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.
        '''
        pass

    def rotate(self, angle, axis, point=None):
        '''
        This function should be overwritten by the derived class, but please
        keep the same inputs, as other parts of the code may depend on it.

        INPUTS:
        angle: float -> Rotation angle, in degrees.

        axis: integer -> Rotation axis [0 for x, 1 for y, 2 for z]

        point: array[3] -> Coordinates of the rotation center. If point=None,
                           the origin should be used.
        '''
        pass

    def project(self,xyz):
        pass
