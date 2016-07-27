from __future__ import division
import numpy as np
from DiscreteSurf.python import cgnsAPI, adtAPI
from mpi4py import MPI

class Surface(object):

    def __init__(self, filename, comm):
        self.comm = comm

        # Read CGNS file
        cgnsAPI.cgnsapi.readcgns(filename, self.comm)

        quadsConn = cgnsAPI.cgnsapi.quadsconn
        triaConn = cgnsAPI.cgnsapi.triaconn
        self.coor = cgnsAPI.cgnsapi.coor
        BBox = np.zeros((3, 2))
        useBBox = False
        self.adtID = 'test_surf'

        adtAPI.adtapi.adtbuildsurfaceadt(self.coor, triaConn, quadsConn, BBox, useBBox, self.comm, self.adtID)

    def inverse_evaluate(self, xyz):
        n = xyz.shape[1]
        allxfs = np.zeros((3, n), order='F')
        allNorms = np.zeros((3, n), order='F')
        arrDonor = np.zeros((1, self.coor.shape[1]), order='F')
        dist2 = np.ones((n), order='F') * 1e9

        adtAPI.adtapi.adtmindistancesearch(xyz, self.adtID, dist2, allxfs, allNorms, arrDonor)

        self.proj_points = allxfs
        self.normals = allNorms

    def get_points(self):
        return self.proj_points

    def get_normals(self):
        return self.normals

    def get_tangent(self):
        # get tangent not implemented in ADT yet for python
        return None

if __name__ == "__main__":

    surf = Surface('DiscreteSurf/python/cube2.cgns', MPI.COMM_WORLD.py2f())

    pts = np.array([[.5, .5, 1.]], order='F').T
    surf.inverse_evaluate(pts)

    print
    print 'Normals:'
    print surf.get_normals()
