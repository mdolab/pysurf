from __future__ import division
import numpy as np
from DiscreteSurf.python import cgnsAPI, adtAPI
from mpi4py import MPI

class Surface(object):

    def __init__(self, inp, comm=MPI.COMM_WORLD.py2f()):

        self.comm = comm

        if isinstance(inp, str):
            # Read CGNS file
            cgnsAPI.cgnsapi.readcgns(inp, self.comm)

            self.coor = cgnsAPI.cgnsapi.coor
            self.quadsConn = cgnsAPI.cgnsapi.quadsconn
            self.triaConn = cgnsAPI.cgnsapi.triaconn

        elif isinstance(inp, dict):
            self.coor = inp['coor']
            self.quadsConn = inp['quadsConn']
            self.triaConn = inp['triaConn']

        else:
            raise SyntaxError('Incorrect input to Surface; must be a cgns file string or dictionary containing coordinate and connection info.')

        self.BBox = np.zeros((3, 2))
        self.useBBox = False
        self.adtID = 'test_surf'

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

    def inverse_evaluate(self, xyz):
        xyz = xyz.T
        n = xyz.shape[1]
        allxfs = np.zeros((3, n), order='F')
        arrDonor = self.nodal_normals
        dist2 = np.ones((n), order='F') * 1e9

        procID, elementType, elementID, uvw, arrInterpol = adtAPI.adtapi.adtmindistancesearch(xyz, self.adtID, dist2, allxfs, arrDonor)

        self.proj_points = allxfs
        self.normals = arrInterpol

    def get_points(self):
        return self.proj_points.T

    def get_normals(self):
        return self.normals.T

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
