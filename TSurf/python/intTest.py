import numpy as np
import matplotlib.pyplot as plt
from intersectionAPI import intersectionapi
from mpl_toolkits.mplot3d import Axes3D

tt = intersectionapi.testtri

def find_int_and_plot(tri1, tri2):

    u0 = tri1[0, :]
    u1 = tri1[1, :]
    u2 = tri1[2, :]

    v0 = tri2[0, :]
    v1 = tri2[1, :]
    v2 = tri2[2, :]

    b, s, e = tt(u0, u1, u2, v0, v1, v2)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    ax.plot([u0[0], u1[0]], [u0[1], u1[1]], [u0[2], u1[2]], color='b')
    ax.plot([u1[0], u2[0]], [u1[1], u2[1]], [u1[2], u2[2]], color='b')
    ax.plot([u2[0], u0[0]], [u2[1], u0[1]], [u2[2], u0[2]], color='b')

    ax.plot([v0[0], v1[0]], [v0[1], v1[1]], [v0[2], v1[2]], color='b')
    ax.plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]], color='b')
    ax.plot([v2[0], v0[0]], [v2[1], v0[1]], [v2[2], v0[2]], color='b')

    ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], color='r')

    plt.show()

###################################
# Orthogonal and share edge
tri1 = np.array([[-1., -1., -1.],
                 [1., 1., 1.],
                 [-.5, -.5, 1.]])

tri2 = np.array([[-.5, -.5, -.5],
                 [.5, .5, .5],
                 [1., 0., 0.]])

find_int_and_plot(tri1, tri2)

###################################
# Orthogonal and intersect
tri1 = np.array([[-2., 0., 0.],
                 [2., 0., 0.],
                 [0., 2., 0.]])

tri2 = np.array([[-2., 1., -1.],
                 [2., 1., -1.],
                 [0., 1., 2.]])

find_int_and_plot(tri1, tri2)

###################################
# Askew and intersect
tri1 = np.array([[-2., 0., 1.],
                 [2., 0., 0.],
                 [0., 2., 2.]])

tri2 = np.array([[-2., 1., -1.],
                 [2., 1., -1.],
                 [0., 1., 2.]])

find_int_and_plot(tri1, tri2)


###################################
# Coplanar and fully enclosed
tri1 = np.array([[-2., 2., 0.],
                 [2., 2., 0.],
                 [0., -2., 0.]])

tri2 = np.array([[-1., 1., 0.],
                 [1., 1., 0.],
                 [0., -1., 0.]])

find_int_and_plot(tri1, tri2)

###################################
# Coplanar and share a side
tri1 = np.array([[-1., 1., 0.],
                 [1., 1., 0.],
                 [0., -1., 0.]])

tri2 = np.array([[-1., 0., 0.],
                 [1., 0., 0.],
                 [0., 1., 0.]])

find_int_and_plot(tri1, tri2)
