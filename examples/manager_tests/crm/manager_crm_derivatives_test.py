# IMPORTS
import pysurf
from mpi4py import MPI
import numpy as np
import unittest
import os
import pickle
import unittest

class CRMDerivTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(CRMDerivTest, self).__init__(*args, **kwargs)

    def test_CRM_intersection(self):

        # MESH PARAMETERS
        numSkinNodes = 129
        LE_spacing = 0.01
        TE_spacing = 0.001

        # TESTING FUNCTION

        os.system('rm *.plt')

        # Load components
        comp1 = pysurf.TSurfGeometry('../../inputs/initial_full_wing_crm4.cgns',['wing','curve_le'])
        comp2 = pysurf.TSurfGeometry('../../inputs/fuselage_crm4.cgns',['fuse'])

        '''
        # I do not know why the rename is not working. It is flipping the elements of both components.
        # I.e.: elements of component 1 go to component 2 and vice-versa.
        name1 = 'wing'
        name2 = 'body'

        comp1.rename(name1)
        comp2.rename(name2)
        '''

        name1 = comp1.name
        name2 = comp2.name

        # Load TE curves and append them to the wing component
        curve_te_upp = pysurf.tsurf_tools.read_tecplot_curves('curve_te_upp.plt_')
        curve_te_low = pysurf.tsurf_tools.read_tecplot_curves('curve_te_low.plt_')
        curveName = curve_te_upp.keys()[0]
        comp1.add_curve(curve_te_upp[curveName])
        curveName = curve_te_low.keys()[0]
        comp1.add_curve(curve_te_low[curveName])

        # Translate the wing
        comp1.translate(0.0, 0.0, 1.5)

        # Create manager object and add the geometry objects to it
        manager0 = pysurf.Manager()
        manager0.add_geometry(comp1)
        manager0.add_geometry(comp2)

        distTol = 1e-7

        #======================================================
        # FORWARD PASS

        def forward_pass(manager):

            '''
            This function will apply all geometry operations to the given manager.
            '''

            # INTERSECT

            # Call intersection function
            intCurveNames = manager.intersect(distTol=distTol)
            intCurveName = intCurveNames[0]

            # SPLIT

            # Split curves based on TE and LE curves
            optionsDict = {'splittingCurves' : [comp1.curves['curve_le'],
                                                comp1.curves['curve_te_upp'],
                                                comp1.curves['curve_te_low']]}

            splitCurveNames = manager.split_intCurve(intCurveName,
                                                     optionsDict,
                                                     criteria='curve')

            # REMESH

            # Find the highest z-coordinate of the entire intersection (vertical position)
            maxZ = -99999
            for curve in splitCurveNames:
                curr_maxZ = np.max(manager.intCurves[curve].coor[:,2])
                maxZ = max(maxZ, curr_maxZ)

            # Now we can identify and remesh each curve properly
            for curveName in splitCurveNames:

                # Get pointer to the curve object
                curve = manager.intCurves[curveName]

                # The trailing edge curve will have less nodes than the other ones
                if curve.numNodes < 20:

                    # This is the trailing edge curve.
                    # Just apply an uniform spacing
                    optionsDict = {'nNewNodes':9}
                    TE_curveName = manager.remesh_intCurve(curveName,optionsDict)

                else:

                    # We have an upper or lower skin curve.
                    # First let's identify if the curve is defined from
                    # LE to TE or vice-versa

                    deltaX = curve.coor[curve.barsConn[0,0],0] - curve.coor[curve.barsConn[-1,-1],0]

                    if deltaX > 0:
                        LE_to_TE = True
                    else:
                        LE_to_TE = False

                    # Compute the highest vertical coordinate of the curve
                    curr_maxZ = np.max(curve.coor[:,2])

                    # Now we can determine if we have upper or lower skin
                    if curr_maxZ < maxZ:

                        # We are at the lower skin

                        if LE_to_TE:

                            optionsDict = {'nNewNodes':numSkinNodes,
                                           'spacing':'hypTan',
                                           'initialSpacing':LE_spacing,
                                           'finalSpacing':TE_spacing}

                            LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

                        else:

                            optionsDict = {'nNewNodes':numSkinNodes,
                                           'spacing':'hypTan',
                                           'initialSpacing':TE_spacing,
                                           'finalSpacing':LE_spacing}

                            LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

                    else:

                        # We are at the upper skin

                        if LE_to_TE:

                            optionsDict = {'nNewNodes':numSkinNodes,
                                           'spacing':'hypTan',
                                           'initialSpacing':LE_spacing,
                                           'finalSpacing':TE_spacing}

                            US_curveName = manager.remesh_intCurve(curveName, optionsDict)

                        else:

                            optionsDict = {'nNewNodes':numSkinNodes,
                                           'spacing':'hypTan',
                                           'initialSpacing':TE_spacing,
                                           'finalSpacing':LE_spacing}

                            US_curveName = manager.remesh_intCurve(curveName, optionsDict)

            # Now we can merge the new curves
            curveNames = [TE_curveName, LS_curveName, US_curveName]
            mergedCurveName = 'intersection'
            manager.merge_intCurves(curveNames, mergedCurveName)

            manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName+'_%03d'%current_pass)

            return mergedCurveName

        # END OF forward_pass
        #======================================================

        current_pass = 0

        # Call the forward pass function to the original manager
        mergedCurveName = forward_pass(manager0)
        current_pass = current_pass + 1

        # DERIVATIVE SEEDS

        # Generate random seeds
        coor1d, curveCoor1d = manager0.geoms[name1].set_randomADSeeds(mode='forward')
        coor2d, curveCoor2d = manager0.geoms[name2].set_randomADSeeds(mode='forward')
        intCoorb = manager0.intCurves[mergedCurveName].set_randomADSeeds(mode='reverse')


        # FORWARD AD

        # Call AD code
        manager0.forwardAD()

        # Get relevant seeds
        intCoord = manager0.intCurves[mergedCurveName].get_forwardADSeeds()

        # REVERSE AD

        # Call AD code
        manager0.reverseAD()

        # Get relevant seeds
        coor1b, curveCoor1b = manager0.geoms[name1].get_reverseADSeeds()
        coor2b, curveCoor2b = manager0.geoms[name2].get_reverseADSeeds()

        # Dot product test
        dotProd = 0.0
        dotProd = dotProd + np.sum(intCoorb*intCoord)
        dotProd = dotProd - np.sum(coor1b*coor1d)
        dotProd = dotProd - np.sum(coor2b*coor2d)
        for curveName in curveCoor1d:
            dotProd = dotProd - np.sum(curveCoor1b[curveName]*curveCoor1d[curveName])
        for curveName in curveCoor2d:
            dotProd = dotProd - np.sum(curveCoor2b[curveName]*curveCoor2d[curveName])

        print('dotProd test (this will be repeated at the end as well)')
        print(dotProd)
        np.testing.assert_almost_equal(dotProd, 0., decimal=13)

        # FINITE DIFFERENCE
        stepSize = 1e-8

        # Apply perturbations to the geometries
        comp1.update(comp1.coor + stepSize*coor1d)
        for curveName in curveCoor1d:
            comp1.curves[curveName].set_points(comp1.curves[curveName].get_points() + stepSize*curveCoor1d[curveName])
        comp2.update(comp2.coor + stepSize*coor2d)
        for curveName in curveCoor2d:
            comp2.curves[curveName].set_points(comp2.curves[curveName].get_points() + stepSize*curveCoor2d[curveName])

        # Create new manager with perturbed components
        manager1 = pysurf.Manager()
        manager1.add_geometry(comp1)
        manager1.add_geometry(comp2)

        # Do the forward pass to the perturbed case
        mergedCurveName = forward_pass(manager1)

        # Get coordinates of the intersection in both cases
        coor0 = manager0.intCurves[mergedCurveName].get_points()
        coor1 = manager1.intCurves[mergedCurveName].get_points()

        # Compute derivatives with FD
        intCoord_FD = (coor1 - coor0)/stepSize

        # Compute difference in derivatives
        FD_error = np.abs(intCoord-intCoord_FD)
        max_FD_error = np.max(FD_error)

        # Export predicted position of the new points
        intCurve = manager0.intCurves[mergedCurveName]
        intCoor = intCurve.get_points() + stepSize*intCoord
        print(np.max(intCoor - coor0))
        np.testing.assert_almost_equal(np.max(intCoor - coor0), 0., decimal=10)
        print(np.max(coor1 - coor0))
        np.testing.assert_almost_equal(np.max(coor1 - coor0), 0., decimal=10)
        intCurve.set_points(intCoor)
        intCurve.export_tecplot(mergedCurveName+'_predicted')

        # Print results
        print('dotProd test')
        print(dotProd)
        np.testing.assert_almost_equal(dotProd, 0., decimal=13)
        print('FD test')
        print(max_FD_error)
        # Pretty loose tolerance here
        np.testing.assert_almost_equal(max_FD_error, 0., decimal=5)

        # import matplotlib.pyplot as plt

        # fig = plt.figure()
        # plt.plot(FD_error[0,:],label='X')
        # plt.plot(FD_error[1,:],label='Y')
        # plt.plot(FD_error[2,:],label='Z')
        # plt.semilogy()
        # plt.legend()
        # plt.show()



if __name__ == "__main__":
    unittest.main()
