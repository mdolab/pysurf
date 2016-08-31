# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np

# TESTING FUNCTION

cube = pysurf.TSurfComponent('../inputs/cube.cgns', MPI.COMM_WORLD)

pts = np.array([[.6, .5, 1.0],
                [.6, .5, 0.1]], order='F')

'''
pts = np.array([[ 0.75      ,  0.75      ,  0.        ],
                [ 0.73717949,  0.73717949,  0.        ],
                [ 0.72435897,  0.72435897,  0.        ],
                [ 0.71153846,  0.71153846,  0.        ],
                [ 0.69871795,  0.69871795,  0.        ],
                [ 0.68589744,  0.68589744,  0.        ],
                [ 0.67307692,  0.67307692,  0.        ],
                [ 0.66025641,  0.66025641,  0.        ],
                [ 0.6474359 ,  0.6474359 ,  0.        ],
                [ 0.63461538,  0.63461538,  0.        ],
                [ 0.62179487,  0.62179487,  0.        ],
                [ 0.60897436,  0.60897436,  0.        ],
                [ 0.59615385,  0.59615385,  0.        ],
                [ 0.58333333,  0.58333333,  0.        ],
                [ 0.57051282,  0.57051282,  0.        ],
                [ 0.55769231,  0.55769231,  0.        ],
                [ 0.54487179,  0.54487179,  0.        ],
                [ 0.53205128,  0.53205128,  0.        ],
                [ 0.51923077,  0.51923077,  0.        ],
                [ 0.50641026,  0.50641026,  0.        ],
                [ 0.49358974,  0.49358974,  0.        ],
                [ 0.48076923,  0.48076923,  0.        ],
                [ 0.46794872,  0.46794872,  0.        ],
                [ 0.45512821,  0.45512821,  0.        ],
                [ 0.44230769,  0.44230769,  0.        ],
                [ 0.42948718,  0.42948718,  0.        ],
                [ 0.41666667,  0.41666667,  0.        ],
                [ 0.40384615,  0.40384615,  0.        ],
                [ 0.39102564,  0.39102564,  0.        ],
                [ 0.37820513,  0.37820513,  0.        ],
                [ 0.36538462,  0.36538462,  0.        ],
                [ 0.3525641 ,  0.3525641 ,  0.        ],
                [ 0.33974359,  0.33974359,  0.        ],
                [ 0.32692308,  0.32692308,  0.        ],
                [ 0.31410256,  0.31410256,  0.        ],
                [ 0.30128205,  0.30128205,  0.        ],
                [ 0.28846154,  0.28846154,  0.        ],
                [ 0.27564103,  0.27564103,  0.        ],
                [ 0.26282051,  0.26282051,  0.        ],
                [ 0.25      ,  0.25      ,  0.        ]])
'''
xyzProj, normProj = cube.project_on_surface(pts)

print 'Original cube projection:'
print xyzProj
print normProj

xyzProj, normProj = cube.project_on_curve(pts)

print 'Original cube edge projection:'
print xyzProj
print normProj

# Change cube and update points
cube.translate(0,0,3)

# Repeat projections

xyzProj, normProj = cube.project_on_surface(pts)

print 'Modified cube projection:'
print xyzProj
print normProj

xyzProj, normProj = cube.project_on_curve(pts)

print 'Modified cube edge projection:'
print xyzProj
print normProj
