from __future__ import division
import adtAPI
import numpy as np

# Define connectivities
coor = np.array([[0.0, 0.0, 0.0],
                 [3.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [3.0, 1.0, 0.0]],order='F').T
triaConn = np.array([[1,2,3],
                    [2,4,3]],order='F').T
quadsConn = np.zeros((0,4),order='F').T
#quadsConn = np.array([[1,2,4,3]],order='F').T

# Define perturbation
stepSize = 1e-6

coord = np.array(np.random.rand(coor.shape[0],coor.shape[1]),order='F')
coord = coord/np.array([np.linalg.norm(coord,axis=1)]).T
nodal_normalsb = np.array(np.random.rand(coor.shape[0],coor.shape[1]),order='F')
nodal_normalsb = nodal_normalsb/np.array([np.linalg.norm(nodal_normalsb,axis=1)]).T

# Use original function
nodal_normals0 = adtAPI.adtapi.adtcomputenodalnormals(coor,
                                                      triaConn,
                                                      quadsConn)

# Compute derivatives of the normal vectors with forward mode
nodal_normals_AD, nodal_normalsd_AD = adtAPI.adtapi.adtcomputenodalnormals_d(coor,
                                                                             coord,
                                                                             triaConn,
                                                                             quadsConn)
# Save original values
nodal_normals0temp = np.array(nodal_normals0,order='F')
nodal_normalsbtemp = np.array(nodal_normalsb,order='F')

# Compute derivatives of the coordinates with reverse mode
coorb_AD = adtAPI.adtapi.adtcomputenodalnormals_b(coor,
                                                  triaConn,
                                                  quadsConn,
                                                  nodal_normals0temp,
                                                  nodal_normalsbtemp)

# Apply perturbation
coor = coor + stepSize*coord

# Use original function at the perturbed point
nodal_normals = adtAPI.adtapi.adtcomputenodalnormals(coor,
                                                     triaConn,
                                                     quadsConn)

# Compute directional derivative
nodal_normalsd_FD = (nodal_normals - nodal_normals0)/stepSize

# Compute dot product test
dotProd = 0.0
print dotProd
dotProd = dotProd + np.sum(coord*coorb_AD)
print dotProd
dotProd = dotProd - np.sum(nodal_normalsd_AD*nodal_normalsb)
print dotProd

# Print results
print 'nodal_normals0'
print nodal_normals0
print 'nodal_normals_AD'
print nodal_normals_AD
print 'nodal_normalsd_AD'
print nodal_normalsd_AD
print 'nodal_normalsd_FD'
print nodal_normalsd_FD
print 'dot product test (should be 0.0)'
print dotProd
