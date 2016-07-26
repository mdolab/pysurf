import CGNSf2py
import adtAPI
from mpi4py import MPI
import numpy as np

# Set CGNS file
fileName = 'cube2.cgns'

# Set communicator
comm = MPI.COMM_WORLD.py2f()

# Read CGNS file
CGNSf2py.cgnsf2py.cgnsf2py_routine(fileName,comm)

quadsConn = CGNSf2py.cgnsf2py.quadsconn
triaConn = CGNSf2py.cgnsf2py.triaconn
coor = CGNSf2py.cgnsf2py.coor
BBox = np.zeros((3, 2))
useBBox = False
adtID = 'test_surf'

adtAPI.adtapi.adtbuildsurfaceadt(coor, triaConn, quadsConn, BBox, useBBox, comm, adtID)

coorPts = np.array([[.5, .5, 1.]]).T

n = coorPts.shape[1]
allxfs = np.zeros((3, n), order='F')
allNorms = np.zeros((3, n), order='F')
arrDonor = np.zeros((1, coor.shape[1]), order='F')
dist2 = np.ones((n), order='F') * 1e9

adtAPI.adtapi.adtmindistancesearch(coorPts, adtID, dist2, allxfs, allNorms, arrDonor)

print
print 'Global coords:'
print allxfs
print
print 'Normals:'
print allNorms
