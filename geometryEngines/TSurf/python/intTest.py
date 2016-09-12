import numpy as np
import matplotlib.pyplot as plt
from intersectionAPI import intersectionapi
from mpl_toolkits.mplot3d import Axes3D

tt = intersectionapi.testtri

def find_int_and_plot(tri1, tri2, ax):

    u0 = tri1[0, :]
    u1 = tri1[1, :]
    u2 = tri1[2, :]

    v0 = tri2[0, :]
    v1 = tri2[1, :]
    v2 = tri2[2, :]

    b, s, e = tt(u0, u1, u2, v0, v1, v2)

    return s, e


###################################

for j in range(1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    tri_list = []
    int_list = []

    # for i in range(10):
    #     tri_list.append(np.random.rand(3, 3))

    #tri_list.append(np.array([[1.0166, 0.3295, .6827],
    #                          [1.006, .3446, .6956],
    #                          [.9943, .3299, .6831]]))
    #tri_list.append(np.array([[1.000, 1.000, 1.0000],
    #                          [.3322, .3588, .3094],
    #                          [.6923, .6492, .6467]]).T)

    V0 = [-2.0, 0.0, 0.0] 
    V1 = [2.0, 0.0, 0.0]
    V2 = [0.0, 2.0, 0.0]
    
    U0 = [-2.0, 1.0, -1.0] 
    U1 = [2.0, 1.0, -1.0]
    U2 = [0.0, 1.0, 2.0]
    
    tri_list.append(np.array([V0,V1,V2]))
    tri_list.append(np.array([U0,U1,U2]))


    for tri1 in tri_list:
        for tri2 in tri_list:
            s, e = find_int_and_plot(tri1, tri2, ax)
            if np.any(s != 0.) and np.any(e != 0.):
                int_list.append([s, e])
            break

    for tri in tri_list:
        u0 = tri[0, :]
        u1 = tri[1, :]
        u2 = tri[2, :]

        ax.plot([u0[0], u1[0]], [u0[1], u1[1]], [u0[2], u1[2]], color='b')
        ax.plot([u1[0], u2[0]], [u1[1], u2[1]], [u1[2], u2[2]], color='b')
        ax.plot([u2[0], u0[0]], [u2[1], u0[1]], [u2[2], u0[2]], color='b')

    for row in int_list:
        s = row[0]
        e = row[1]
        ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], color='r')

    print s
    print e

    plt.show()
