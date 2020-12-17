from GeoMACH.BSE import BSEmodel
from scipy.sparse.linalg import bicgstab
import numpy


class Surface(object):
    def __init__(self, name, surf):

        self.name = name
        self.type = "surface"

        self.bse = BSEmodel([surf])

        nu, nv = surf.shape[:2]
        mu = max(4, int(nu / 4))
        mv = max(4, int(nv / 4))
        self.bse.set_bspline_option("num_cp", 0, "u", mu)
        self.bse.set_bspline_option("num_cp", 0, "v", mv)
        self.bse.set_bspline_option("num_pt", 0, "u", nu)
        self.bse.set_bspline_option("num_pt", 0, "v", nv)

        self.bse.assemble()

        cp = self.bse.vec["cp_str"].array
        pt = self.bse.vec["pt_str"].array
        jac = self.bse.jac["d(pt_str)/d(cp_str)"]

        self.bse.vec["pt_str"](0)[:, :, :] = surf
        for ind in range(3):
            cp[:, ind] = bicgstab(jac.T.dot(jac), jac.T.dot(pt[:, ind]), tol=1.0e-8)[0]

    def project(self, xyz, dist2, xyzProj, normProj):
        self.bse.compute_projection("proj", xyz, ndim=3)

        # Get points
        self.bse.apply_jacobian("proj", "d(proj)/d(cp_str)", "cp_str")
        xyzProjNew = numpy.array(self.bse.vec["proj"].array)

        # Get normals
        self.bse.apply_jacobian("proj", "d(proj_du)/d(cp_str)", "cp_str")
        du = numpy.array(self.bse.vec["proj"].array)
        self.bse.apply_jacobian("proj", "d(proj_dv)/d(cp_str)", "cp_str")
        dv = numpy.array(self.bse.vec["proj"].array)
        cross = numpy.cross(du, dv)
        norms = numpy.sum(cross ** 2, axis=1) ** 0.5
        for ind in range(3):
            cross[:, ind] = cross[:, ind] / norms

        # Get new distances
        dist2New = numpy.linalg.norm(xyz - xyzProjNew, axis=1) ** 2

        # Update only if we find a better candidate
        for ii in range(dist2New.shape[0]):
            if dist2New[ii] < dist2[ii]:
                xyzProj[ii, :] = xyzProjNew[ii, :]
                normProj[ii, :] = cross[ii, :]
                dist2[ii] = dist2New[ii]


class Curve(object):
    def __init__(self, name, curve_2D):

        self.name = name
        self.type = "curve"

        nu = curve_2D.shape[0]
        curve = numpy.zeros((nu, 4, 3))
        for i in range(4):
            curve[:, i, :] = curve_2D

        self.bse = BSEmodel([curve])

        nu, nv = curve.shape[:2]
        mu = max(4, int(nu / 4))
        mv = max(4, int(nv / 4))
        self.bse.set_bspline_option("num_cp", 0, "u", mu)
        self.bse.set_bspline_option("num_cp", 0, "v", mv)
        self.bse.set_bspline_option("num_pt", 0, "u", nu)
        self.bse.set_bspline_option("num_pt", 0, "v", nv)

        self.bse.assemble()

        cp = self.bse.vec["cp_str"].array
        pt = self.bse.vec["pt_str"].array
        jac = self.bse.jac["d(pt_str)/d(cp_str)"]

        self.bse.vec["pt_str"](0)[:, :, :] = curve
        for ind in range(3):
            cp[:, ind] = bicgstab(jac.T.dot(jac), jac.T.dot(pt[:, ind]), tol=1.0e-8)[0]

    def project(self, xyz, dist2, xyzProj, tanProj):
        self.bse.compute_projection("proj", xyz, ndim=3)

        # Get points
        self.bse.apply_jacobian("proj", "d(proj)/d(cp_str)", "cp_str")
        xyzProjNew = self.bse.vec["proj"].array.copy()

        # Get tangent
        self.bse.apply_jacobian("proj", "d(proj_du)/d(cp_str)", "cp_str")
        du = self.bse.vec["proj"].array.copy()
        du = du / numpy.sum(du ** 2, axis=1) ** 0.5

        # Get new distances
        dist2New = numpy.linalg.norm(xyz - xyzProjNew, axis=1) ** 2

        for ii in range(dist2New.shape[0]):
            if dist2New[ii] < dist2[ii]:
                xyzProj[ii, :] = xyzProjNew[ii, :]
                tanProj[ii, :] = du[ii, :]
                dist2[ii] = dist2New[ii]


if __name__ == "__main__":
    import numpy

    # Define surface points
    s0 = 0.0
    s1 = 3.0
    radius = 0.5
    nu, nv = 100, 100
    pts = numpy.zeros((nu, nv, 3))
    for i in range(nu):
        for j in range(nv):
            theta = numpy.pi * j / (nv - 1)
            pts[i, j, 0] = s0 + (s1 - s0) * (i / (nu - 1))
            pts[i, j, 1] = -radius * numpy.cos(theta)
            pts[i, j, 2] = radius * numpy.sin(theta)

    # Create geometry object
    s = Surface("cylinder", pts)

    xyz = numpy.zeros((1, 3))
    xyz[0, :] = [0.5, 0.5, 1.0]

    dist2 = 1e10 * numpy.ones(xyz.shape[0])
    xyzProj = numpy.zeros((xyz.shape[0], 3))
    normProj = numpy.zeros((xyz.shape[0], 3))

    s.project(xyz, dist2, xyzProj, normProj)
    print(xyzProj)
    print(normProj)

    # Define surface points
    s0 = 100
    s1 = 100
    nu = 10
    curve1_pts = numpy.zeros((nu, 3))
    for i in range(nu):
        curve1_pts[i, 0] = -5.0 + i
        curve1_pts[i, 1] = 2.5

    s = Curve("curve", curve1_pts)
    xyz = numpy.zeros((1, 3))
    xyz[0, :] = [-3.71661445e00, 1.85840743e00, 1.72689747e00]

    dist2 = 1e10 * numpy.ones(xyz.shape[0])
    xyzProj = numpy.zeros((xyz.shape[0], 3))
    normProj = numpy.zeros((xyz.shape[0], 3))

    s.project(xyz, dist2, xyzProj, normProj)
    print(xyzProj)
    print(normProj)
