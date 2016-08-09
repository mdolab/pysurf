'''
This file contains functions to intersect ADT objects
'''

# IMPORTS
import intersectionAPI as iapi
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

    # Call intersection function for each pair
    for ii in range(numComponents):
        for jj in range(ii+1,numComponents):
            compute_pair_intersection(TSurfComponentList[ii],
                                      TSurfComponentList[jj],
                                      distTol)

def compute_pair_intersection(TSurfComponentA, TSurfComponentB, distTol):

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

    print coor
    print barsConn

    # Return intersection FE data
    # return coor, barsConn
