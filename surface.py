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

    def get_points(self, dmy):
        self.bse.apply_jacobian('proj', 'd(proj)/d(cp_str)', 'cp_str')
        return self.bse.vec['proj'].array

    def get_normals(self, dmy):
        self.bse.apply_jacobian('proj', 'd(proj_du)/d(cp_str)', 'cp_str')
        du = numpy.array(self.bse.vec['proj'].array)
        self.bse.apply_jacobian('proj', 'd(proj_dv)/d(cp_str)', 'cp_str')
        dv = numpy.array(self.bse.vec['proj'].array)
        cross = numpy.cross(du, dv)
        norms = numpy.sum(cross**2, axis=1)**0.5
        for ind in xrange(3):
            cross[:, ind] = cross[:, ind] / norms
        return cross


if __name__ == '__main__':
    import numpy

    nu, nv = 10, 10
    array = numpy.zeros((nu, nv, 3))
    for i in xrange(nu):
        for j in xrange(nv):
            array[i, j, 0] = i / nu
            array[i, j, 1] = j / nv
    s = Surface(array)

    xyz = numpy.zeros((1, 3))
    xyz[0, :] = [0.5, 0.5, 10]
    s.inverse_evaluate(xyz)
    print s.get_points(None)
    print s.get_normals(None)
