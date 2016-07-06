# IMPORTS

from __future__ import division
import hypsurf
from numpy import array

# EXECUTION

# Set problem
rBaseline = array([0,0,0, 1,1,0, 2,2,0, 3,2,0, 4,1,0, 5,0,0])
NBaseline = array([[0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1], [0,0,1]]).T
bc1 = 'splay'
bc2 = 'splay'
sBaseline = 0.2
numLayers = 10
extension = 2

# Generate surface mesh
X,Y,Z = hypsurf.create_mesh(rBaseline, NBaseline, bc1, bc2, sBaseline, numLayers, extension, \
                            epsE0 = 1.0, theta = 0.0, \
                            alphaP0 = 0.25, num_smoothing_passes = 0, \
                            nuArea = 0.16, num_area_passes = 0, \
                            sigmaSplay = 0.3, \
                            metricCorrection = False, \
                            ratioGuess = 20)

print 'X'
print X
print 'Y'
print Y

# Plot results
fig = hypsurf.plot_grid(X,Y,Z,show=True)
