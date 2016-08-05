'''
This file contains functions to intersect ADT objects
'''

# IMPORTS
import intersectionAPI as iapi

# FUNCTIONS
def compute_intersections(ADTComponentList):

    '''
    This function will just compute pair-wise intersections
    for all components given in ADTComponentList
    '''

    # Get number of components
    numComponents = len(ADTComponentList)

    # Stop if user gives only one component
    if numComponents < 2:
        print 'ERROR: Cannot compute intersections with just one component'
        quit()

    # Call intersection function for each pair
    for ii in range(numComponents):
        for jj in range(ii+1,numComponents):
            compute_pair_intersection(ADTComponentList[ii],
                                      ADTComponentList[jj])

def compute_pair_intersection(ADTComponentA, ADTComponentB):

    '''
    This function finds intersection curves between components A and B
    '''

    # Call Fortran code to find intersections
    iapi.intersectionapi.computeintersection(ADTComponentA.coor,
                                             ADTComponentA.triaConn,
                                             ADTComponentA.quadsConn,
                                             ADTComponentB.coor,
                                             ADTComponentB.triaConn,
                                             ADTComponentB.quadsConn)
