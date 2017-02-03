"""
This script tests the projection algorithm for a TSurf geometry object.

First it checks the projection on the original cube, then it translate the cube
3 units in the z-direction and again checks the projection.

"""

# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np

# Load geometry
'''
curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('../manager_tests/crm/curve_te_low.plt_')
#curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('test_curve.plt')
curveName = curve_te_low.keys()[0]
curve = curve_te_low[curveName]
curve.flip()

print curve.barsConn
barsConn = curve.barsConn[:,1:8]
coor = []
coor.append(curve.coor[:,barsConn[0,0]])
for ii in range(barsConn.shape[1]):
    coor.append(curve.coor[:,barsConn[-1,ii]])
coor = np.array(coor).T

curve = pysurf.TSurfCurve(curveName, coor, barsConn)
curve.export_tecplot('test_curve')
print coor
'''

initCurveName = 'test_curve'
'''
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.0],
                 [0.7, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.5, 0.5, 0.0]],order='F').T
'''
'''
coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.3],
                 [0.7, 0.0, 0.6],
                 [1.0, 0.0, 0.3],
                 [1.0, 0.3, 0.0],
                 [1.0, 0.7, -0.3],
                 [1.0, 1.0, -0.3],
                 [0.5, 0.5, 0.0]],order='F').T
'''

coor = np.array([[0.0, 0.0, 0.0],
                 [0.3, 0.0, 0.3]],order='F').T

'''
coor = np.array([[ 36.55844879,  36.57482147,  36.59119415,  36.60756302,  36.6239357,
   36.64030838,  36.65667725,  36.67304993],
 [  0.19090894,   0.28636342,   0.38181788,   0.47727233,   0.57272685,
    0.6681813,    0.76363575,   0.85909021],
 [  2.89357114,  2.9156549,    2.93773842,   2.95982218,   2.9819057,
    3.00398946,   3.02607298,   3.04815674]],order='F')

coor = np.array(np.round(coor,1),order='F')
'''
'''
barsConn = np.array([[4,5],
                     [6,7],
                     [2,3],
                     [3,4],
                     [1,2],
                     [5,6],
                     [7,8],
                     [8,1]],dtype='int32',order='F').T
'''
'''
barsConn = np.array([[1,2],
                     [2,3],
                     [3,4],
                     [4,5],
                     [5,6],
                     [6,7],
                     [7,8],
                     [8,1]],dtype='int32',order='F').T
'''
barsConn = np.array([[1,2]],dtype='int32',order='F').T

curve = pysurf.TSurfCurve(initCurveName, coor, barsConn)
curve.export_tecplot('test_curve')

# Define generic function to compute projections

def computeCurveProjections(xyz, xyzd, coord, xyzProjb, tanProjb, coor=None):

    # Replace baseline coordinates if the user provided it
    if coor is not None:
        curve.set_points(coor)

    # Get number of points                                                                        
    nPoints = xyz.shape[1]

    # Initialize references if user provided none                                                 
    dist2 = np.ones((nPoints),order='F')*1e10
    xyzProj = np.zeros((3,nPoints),order='F')
    tanProj = np.zeros((3,nPoints),order='F')
    elemIDs = np.zeros((nPoints),dtype='int32',order='F')
    '''
    numPts = xyz.shape[0]
    dist2 = np.ones(numPts)*1e10
    xyzProj = np.zeros((numPts,3))
    tanProj = np.zeros((numPts,3))
    elemIDs = np.zeros((numPts),dtype='int32')
    '''
    # Call projection algorithm
    curveMask = curve.project(xyz, dist2, xyzProj, tanProj, elemIDs)

    # Set derivative seeds
    curve.set_forwardADSeeds(coord)
    
    # Call derivatives code in forward mode
    xyzProjd = np.zeros(xyzProj.shape,order='F')
    tanProjd = np.zeros(tanProj.shape,order='F')
    curve.project_d(xyz, xyzd, xyzProj, xyzProjd, tanProj, tanProjd, elemIDs, curveMask)

    # Call derivatives code in backward mode
    xyzb = np.zeros(xyz.shape,order='F')
    curve.project_b(xyz, xyzb, xyzProj, xyzProjb, tanProj, tanProjb, elemIDs, curveMask)
    coorb = curve.get_reverseADSeeds()

    # Print results
    print
    print 'Curve projection results:'
    print xyzProj
    print tanProj
    print

    # Return results
    return xyzProj, xyzProjd, tanProj, tanProjd, xyzb, coorb

# BACK TO MAIN PROGRAM

# Define points
#xyz = np.array([[32.0, 1.5, 2.1]],order='F').T
#xyz = np.array([[0.8, 0.0, 0.1]],order='F').T
xyz = np.array([[0.3, 0.0, 0.1]],order='F').T

# Define derivatives and normalize them to a given step size
stepSize = 1e-7

xyzd = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
xyzd = xyzd/np.sqrt(np.sum(xyz**2))
curve.set_randomADSeeds(mode='forward')
coord = curve.get_forwardADSeeds()

xyzProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')
tanProjb = np.array(np.random.rand(xyz.shape[0],xyz.shape[1]),order='F')

xyzd[:,:] = 0.0
xyzd[0,0] = 1.0
#coord[:,:] = 1.0
#xyzProjb[:,:] = 0.0
#xyzProjb[2,0] = 1.0
#tanProjb[:,:] = 0.0
#tanProjb[0,0] = 1.0

xyzProjb_copy = xyzProjb.copy()
tanProjb_copy = tanProjb.copy()

#############################
# CURVE PROJECTION
#############################

# Call projection function at the original point
xyzProj0, xyzProjd_AD, tanProj0, tanProjd_AD, xyzb_AD, coorb_AD = computeCurveProjections(xyz,
                                                                                          xyzd, coord,
                                                                                          xyzProjb, tanProjb)

print 'derivatives'
print xyzProjd_AD
print tanProjd_AD
print xyzb_AD
print coorb_AD
print coor
# Create coordinates of perturbed points
coor = curve.get_points() + stepSize*coord

# Call projection function at the perturbed point
xyzf = xyz+stepSize*xyzd
xyzProj, dummy0, tanProj, dummy1, dummy2, dummy3 = computeCurveProjections(xyzf,
                                                                           xyzd, coord,
                                                                           xyzProjb, tanProjb,
                                                                           coor)

# Compute directional derivatives with finite differencing
xyzProjd_FD = (xyzProj - xyzProj0)/stepSize
tanProjd_FD = (tanProj - tanProj0)/stepSize

# Compute dot product test
dotProd = 0.0
dotProd = dotProd + np.sum(xyzd*xyzb_AD)
print dotProd
dotProd = dotProd + np.sum(coord*coorb_AD)
print dotProd
dotProd = dotProd - np.sum(xyzProjd_AD*xyzProjb_copy)
print dotProd
dotProd = dotProd - np.sum(tanProjd_AD*tanProjb_copy)

# Dot product test
print 'dot product (should be 0.0)'
print dotProd

# Compare directional derivatives
print 'CURVE PROJECTION'
print ''
print 'xyzProjd_AD'
print xyzProjd_AD
print 'xyzProjd_FD'
print xyzProjd_FD
print 'tanProjd_AD'
print tanProjd_AD
print 'tanProjd_FD'
print tanProjd_FD
