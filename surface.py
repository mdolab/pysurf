from __future__ import division
from GeoMACH.BSE import BSEmodel
from scipy.sparse.linalg import bicgstab
import numpy



class Surface(object):

    def __init__(self, surf):
        self.bse = BSEmodel([surf])

        nu, nv = surf.shape[:2]
        mu = max(4, int(nu/4))
        mv = max(4, int(nv/4))
        self.bse.set_bspline_option('num_cp', 0, 'u', mu)
        self.bse.set_bspline_option('num_cp', 0, 'v', mv)
        self.bse.set_bspline_option('num_pt', 0, 'u', nu)
        self.bse.set_bspline_option('num_pt', 0, 'v', nv)

        self.bse.assemble()

        cp = self.bse.vec['cp_str'].array
        pt = self.bse.vec['pt_str'].array
        jac = self.bse.jac['d(pt_str)/d(cp_str)']

        self.bse.vec['pt_str'](0)[:, :, :] = surf
        for ind in xrange(3):
            cp[:, ind] = bicgstab(jac.T.dot(jac), jac.T.dot(pt[:, ind]), tol=1.0e-8)[0]

    def inverse_evaluate(self, xyz):
        self.bse.compute_projection('proj', xyz, ndim=3)
        return None

    def get_points(self):
        self.bse.apply_jacobian('proj', 'd(proj)/d(cp_str)', 'cp_str')
        return self.bse.vec['proj'].array

    def get_normals(self):
        self.bse.apply_jacobian('proj', 'd(proj_du)/d(cp_str)', 'cp_str')
        du = numpy.array(self.bse.vec['proj'].array)
        self.bse.apply_jacobian('proj', 'd(proj_dv)/d(cp_str)', 'cp_str')
        dv = numpy.array(self.bse.vec['proj'].array)
        cross = numpy.cross(du, dv)
        norms = numpy.sum(cross**2, axis=1)**0.5
        for ind in xrange(3):
            cross[:, ind] = cross[:, ind] / norms
        return cross

    def get_tangent(self):
        self.bse.apply_jacobian('proj', 'd(proj_du)/d(cp_str)', 'cp_str')
        du = numpy.array(self.bse.vec['proj'].array)
        du = du / numpy.sum(du**2, axis=1)**0.5
        return du


if __name__ == '__main__':
    import numpy

    # nu, nv = 11, 11
    # array = numpy.zeros((nu, nv, 3))
    # for i in xrange(nu):
    #     for j in xrange(nv):
    #         array[i, j, 0] = i / (nu-1)
    #         array[i, j, 1] = j / (nv-1)
    #         array[i, j, 2] = j / (nv-1)
    # s = Surface(array)
    #
    # xyz = numpy.zeros((1, 3))
    # xyz[0, :] = [0.5, 0.5, 1.0]
    # s.inverse_evaluate(xyz)
    # print s.get_points(None)
    # print s.get_normals(None)

    # Define surface points
    s0 = 100
    s1 = 100
    nu, nv = 10, 4
    curve1_pts = numpy.zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            curve1_pts[i, j, 0] = -5. + i
            curve1_pts[i, j, 1] = 2.5
    for i in range(4):
        curve1_pts[:, i, :] = numpy.array([[ -3.71661445e+00,   1.85840743e+00,   1.72689747e+00],
                                           [ -3.13941694e+00,   1.96221597e+00,   1.37089540e+00],
                                           [ -2.48185790e+00,   2.06828033e+00,   1.04407268e+00],
                                           [ -1.73311129e+00,   2.16645113e+00,   7.70094923e-01],
                                           [ -8.95715540e-01,   2.23938614e+00,   5.82100990e-01],
                                           [  3.42931180e-05,   2.26693406e+00,   5.14241593e-01],
                                           [  8.95798123e-01,   2.23942268e+00,   5.82005704e-01],
                                           [  1.73321982e+00,   2.16650657e+00,   7.69937250e-01],
                                           [  2.48198542e+00,   2.06833871e+00,   1.04388412e+00],
                                           [  3.13954714e+00,   1.96226976e+00,   1.37069305e+00]])


    s = Surface(curve1_pts)
    xyz = numpy.zeros((1, 3))
    xyz[0, :] = [ -3.71661445e+00,   1.85840743e+00,   1.72689747e+00]

    s.inverse_evaluate(xyz)
    print s.get_points()
