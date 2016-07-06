# IMPORTS

from __future__ import division
import hypsurf
from numpy import array, cos, sin, linspace, pi, zeros, vstack, ones
from numpy.random import rand

# EXECUTION

# Set problem
n = 41
nums = linspace(2*pi, 0, n)
x = cos(nums)
y = sin(nums)
z = ones(n)
rBaseline = vstack([x, y, z]).T.ravel()
NBaseline = array([array([0, 0, 1])] * n).T
bc1 = 'continuous'
bc2 = 'continuous'
sBaseline = 0.02
numLayers = 20
extension = 10

# Generate surface mesh
X,Y,Z = hypsurf.create_mesh(rBaseline, NBaseline, bc1, bc2, sBaseline, numLayers, extension, \
                            epsE0 = 1.0, theta = 0.0, \
                            alphaP0 = 0.25, num_smoothing_passes = 0, \
                            nuArea = 0.16, num_area_passes = 0, \
                            sigmaSplay = 0.3, \
                            metricCorrection = False, \
                            ratioGuess = 20)

# Plot results
fig = hypsurf.plot_grid(X,Y,Z,show=True)
