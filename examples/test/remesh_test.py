# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import pysurf
import unittest


# Create simple curve and store within the class
coor = np.array([[0.0, 0.0, 0.0],
                 [0.5, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 0.5, 0.0],
                 [1.0, 1.0, 0.0]],order='F').T

# Create curve object
curve = pysurf.tsurf_tools.create_curve_from_points(coor, 'test', periodic=False)
curve.export_tecplot('original')

# Remesh the curve
optionsDict = {
    'nNewNodes':13,
    'spacing':'linear',
    'initialSpacing':0.005,
    'finalSpacing':0.005,
}
remeshedCurve = curve.remesh(**optionsDict)
remeshedCurve.export_tecplot('original_remesh')

# Assign derivative seeds
coord = curve.set_randomADSeeds(mode='forward')
remeshedCoorb = remeshedCurve.set_randomADSeeds(mode='reverse')

# Adjust seeds
coord[:,:] = 0.0
coord[1,2] = -1.0
curve.set_forwardADSeeds(coord)

# Forward AD
curve.remesh_d(remeshedCurve, **optionsDict)

# Store relevant seeds
remeshedCoord = remeshedCurve.get_forwardADSeeds()

# Reverse AD
curve.remesh_b(remeshedCurve, **optionsDict)

# Store relevant seeds
coorb = curve.get_reverseADSeeds()

# DOT PRODUCT TEST
dotProd = 0.0
dotProd = dotProd + np.sum(coord*coorb)
dotProd = dotProd - np.sum(remeshedCoord*remeshedCoorb)

print('dot product test')
print(dotProd)

# Define step size for finite difference test
stepSize = 1e-7

# Apply perturbations for finite difference tests
curve.set_points(curve.get_points() + stepSize*coord)
curve.export_tecplot('perturbed')

# Remesh the curve
remeshedCurve2 = curve.remesh(**optionsDict)
remeshedCurve2.export_tecplot('perturbed_remesh')

# Compare coordinates to compute derivatives
remeshedCoor0 = remeshedCurve.get_points()
remeshedCoorf = remeshedCurve2.get_points()
remeshedCoord_FD = (remeshedCoorf - remeshedCoor0)/stepSize

# Compare derivatives
max_FD_error = np.max(np.abs(remeshedCoord - remeshedCoord_FD))

print('finite difference test')
print(max_FD_error)

# Export the curve predicted by AD
remeshedCurve.set_points(remeshedCurve.get_points() + stepSize*remeshedCoord)
remeshedCurve.export_tecplot('predicted_remesh')

# Plot errors

def view_mat(mat):
    """ Helper function used to visually examine matrices. """
    import matplotlib.pyplot as plt
    if len(mat.shape) > 2:
        mat = np.sum(mat, axis=2)
    # print "Cond #:", np.linalg.cond(mat)
    im = plt.imshow(mat, interpolation='none')
    plt.colorbar(im, orientation='horizontal')
    plt.show()

view_mat(np.abs(remeshedCoord - remeshedCoord_FD))

print(remeshedCoord)

print(remeshedCoord_FD)
