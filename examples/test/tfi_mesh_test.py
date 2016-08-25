# IMPORTS
from __future__ import division
import pysurf
from mpi4py import MPI
import numpy as np

# TESTING FUNCTION

nu = 7
nv = 7

X = np.zeros((nu,nv))
Y = np.zeros((nu,nv))
Z = np.zeros((nu,nv))

X[:,0] = np.arange(nu)
X[:,-1] = np.arange(nu)
X[-1,:] = nu-1

Y[0,:] = np.arange(nv)
Y[-1,:] = np.arange(nv)
Y[:,-1] = nv-1

curve1 = np.vstack([X[:,0], Y[:,0], Z[:,0]]).T
curve2 = np.vstack([X[:,-1], Y[:,-1], Z[:,-1]]).T
curve3 = np.vstack([X[0,:], Y[0,:], Z[0,:]]).T
curve4 = np.vstack([X[-1,::-1], Y[-1,::-1], Z[-1,::-1]]).T

leftCurve, bottomCurve, rightCurve, topCurve = pysurf.tfi_mesh.link_curves(curve1,
                                                                           curve2,
                                                                           curve3,
                                                                           curve4)

refComp = pysurf.TSurfComponent('../inputs/plate.cgns')

mesh = pysurf.tfi_mesh.TFIMesh(leftCurve, bottomCurve, rightCurve, topCurve, refComp)

mesh.generate_mesh()

print mesh.get_mesh()
