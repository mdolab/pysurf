'''
This file contains functions to intersect ADT objects
'''

# IMPORTS
import intersectionAPI as iapi
import tsurf_geometry
import numpy as np

# FUNCTIONS
def compute_intersections(TSurfComponentList,distTol=1e-7):

    '''
    This function will just compute pair-wise intersections
    for all components given in TSurfComponentList.

    distTol is a distance tolerance used to merge nearby nodes when
    generating the intersection finite element data.
    '''

    # Get number of components
    numComponents = len(TSurfComponentList)

    # Stop if user gives only one component
    if numComponents < 2:
        print 'ERROR: Cannot compute intersections with just one component'
        quit()

    # Initialize list of intersections
    Intersections = []

    # Call intersection function for each pair
    for ii in range(numComponents):
        for jj in range(ii+1,numComponents):


            coor, barsConn = _compute_pair_intersection(TSurfComponentList[ii],
                                                       TSurfComponentList[jj],
                                                       distTol)

            # barsConn is a list of curves. Each element of the list brings
            # FE data for a continuous curve. If the intersection generates
            # multiple disconnect cruves, barsConn will be a list with
            # multiple elements.
            # For each curve (and therefore for each barsConn element), we
            # will initialize a Curve object to hold intersection FE data.

            for currConn in barsConn:

                # Create new curve object
                newCurve = tsurf_geometry.Curve(coor, currConn)

                # Initialize curve object and append it to the list
                Intersections.append(newCurve)

    # Return all intersections
    return Intersections

#=================================================================

def _compute_pair_intersection(TSurfComponentA, TSurfComponentB, distTol):

    '''
    This function finds intersection curves between components A and B
    '''

    # Call Fortran code to find intersections
    iapi.intersectionapi.computeintersection(TSurfComponentA.coor,
                                             TSurfComponentA.triaConn,
                                             TSurfComponentA.quadsConn,
                                             TSurfComponentB.coor,
                                             TSurfComponentB.triaConn,
                                             TSurfComponentB.quadsConn,
                                             distTol)

    # Retrieve results from Fortran
    coor = np.array(iapi.intersectionapi.coor)
    barsConn = np.array(iapi.intersectionapi.barsconn)

    barsConn = np.sort(barsConn, axis=0)
    barsConn = np.vstack({tuple(col) for col in barsConn.T}).T

    # Sort FE data. After this step, barsConn may become a list if we
    # have two disconnect intersection curves.
    barsConn = tsurf_geometry.FEsort(barsConn.T.tolist())

    # Return intersection FE data
    return coor, barsConn
