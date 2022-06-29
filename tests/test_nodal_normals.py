from pysurf import adtAPI, adtAPI_cs
import numpy as np
import unittest

# Set random seed for derivative tests
np.random.seed(123)


class TestNodalNormals(unittest.TestCase):

    N_PROCS = 2

    def test_nodal_normals(self):
        # Define connectivities
        coor = np.array(
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [3.0, 1.0, 0.0], [4.0, 0.0, 0.0], [4.0, 1.0, 0.0]],
            order="F",
        ).T
        triaConn = np.array([[1, 2, 3], [2, 4, 3]], order="F").T
        quadsConn = np.array([[2, 5, 6, 4]], order="F").T

        # Define perturbations
        stepSize_FD = 1e-7
        stepSize_CS = 1e-200

        # Set random AD seeds
        coord = np.array(np.random.rand(coor.shape[0], coor.shape[1]), order="F")
        coord = coord / np.array([np.linalg.norm(coord, axis=1)]).T
        nodal_normalsb = np.array(np.random.rand(coor.shape[0], coor.shape[1]), order="F")
        nodal_normalsb = nodal_normalsb / np.array([np.linalg.norm(nodal_normalsb, axis=1)]).T

        # Compute unperturbed normals
        nodal_normals0 = adtAPI.adtapi.adtcomputenodalnormals(coor, triaConn, quadsConn)

        # Compute derivatives of the normal vectors with forward mode
        _, nodal_normalsd = adtAPI.adtapi.adtcomputenodalnormals_d(coor, coord, triaConn, quadsConn)

        # Compute derivatives of the coordinates with reverse mode
        coorb_AD = adtAPI.adtapi.adtcomputenodalnormals_b(coor, triaConn, quadsConn, nodal_normals0, nodal_normalsb)

        # Apply FD perturbation
        coor_FD = coor + stepSize_FD * coord

        # Compute FD derivative
        nodal_normals_pert = adtAPI.adtapi.adtcomputenodalnormals(coor_FD, triaConn, quadsConn)
        nodal_normals_FD = (nodal_normals_pert - nodal_normals0) / stepSize_FD

        # Apply CS perturbation
        coor_CS = coor + stepSize_CS * coord * 1j

        # Compute CS derivative
        nodal_normals_pert = adtAPI_cs.adtapi.adtcomputenodalnormals(coor_CS, triaConn, quadsConn)
        nodal_normals_CS = np.imag(nodal_normals_pert) / stepSize_CS

        # Dot product test
        np.testing.assert_allclose(np.sum(coord * coorb_AD), np.sum(nodal_normalsd * nodal_normalsb), rtol=1e-15)

        # Compare AD to FD
        np.testing.assert_allclose(nodal_normalsd, nodal_normals_FD, atol=3e-8, rtol=1e-14)

        # Compare AD to CS
        np.testing.assert_allclose(nodal_normalsd, nodal_normals_CS, atol=1e-14, rtol=1e-14)


if __name__ == "__main__":
    unittest.main()
