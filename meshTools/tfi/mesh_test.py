# IMPORTS
import numpy as np
from tfiMeshTools import computeDerivatives, computeMesh, exportPlot3d
import os

nu = 10
nv = 41

X = np.zeros((nu, nv))
Y = np.zeros((nu, nv))
Z = np.zeros((nu, nv))

X[:, 0] = np.arange(nu)
X[:, -1] = np.arange(nu)
X[-1, :] = nu - 1

Y[0, :] = np.arange(nv)
Y[-1, :] = np.arange(nv)
Y[:, -1] = nv - 1

s0u = 0.3 * np.ones((nu, 2))
s0v = 0.3 * np.ones((2, nv))

N0u = np.zeros((X.shape[0], 6))
N0u[:, 2] = 1.0
N0u[:, 5] = 1.0

N0v = np.zeros((6, X.shape[1]))
N0v[2, :] = 1.0
N0v[5, :] = 1.0

derivDict = computeDerivatives(X, Y, Z, s0u, s0v, N0u, N0v)

print(derivDict)

computeMesh(X, Y, Z, derivDict)

exportPlot3d(X, Y, Z, "output.xyz")

# print(X)
# print(Y)
# print(Z)

os.system("tec360 layout_mesh.lay")
