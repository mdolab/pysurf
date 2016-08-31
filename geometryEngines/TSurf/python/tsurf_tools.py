from __future__ import division
import numpy as np
import cgnsAPI, adtAPI, utilitiesAPI, intersectionAPI, curveSearch
from mpi4py import MPI
import tsurf_component

'''
TODO:

- Move Curve class to tsurf_component.pt
- Add variable angle threshold to curve feature extraction and curve split
- Add code to remove unused points from surfaces and curves

'''

#===========================================
# AUXILIARY READING FUNCTIONS
#===========================================

def getCGNSsections(inputFile, comm=MPI.COMM_WORLD):

    '''
    This function opens a CGNS file, reads its sections, and returns a
    dictionary with selected sections. This code also initializes ADT
    for each section.
    '''

    # Read CGNS file
    cgnsAPI.cgnsapi.readcgns(inputFile, comm.py2f())

    # Retrieve data from the CGNS file.
    # We need to do actual copies, otherwise data will be overwritten if we read another
    # CGNS file.
    coor = np.array(cgnsAPI.cgnsapi.coor)
    triaConn = np.array(cgnsAPI.cgnsapi.triaconn)
    quadsConn = np.array(cgnsAPI.cgnsapi.quadsconn)
    barsConn = np.array(cgnsAPI.cgnsapi.barsconn)
    surfTriaPtr = np.array(cgnsAPI.cgnsapi.surftriaptr)
    surfQuadsPtr = np.array(cgnsAPI.cgnsapi.surfquadsptr)
    curveBarsPtr = np.array(cgnsAPI.cgnsapi.curvebarsptr)
    surfNames = cgnsAPI.cgnsapi.surfnames.copy()
    curveNames = cgnsAPI.cgnsapi.curvenames.copy()

    # Now we deallocate variables on the Fortran side
    cgnsAPI.cgnsapi.releasememory()

    # Format strings coming to Fortran into a Python list
    surfNames = formatStringArray(surfNames)
    curveNames = formatStringArray(curveNames)

    # Initialize dictionary
    sectionDict = {}

    # Now we split the connectivity arrays according to the sections
    iSurf = 0
    for surf in surfNames:

        # Get indices to slice connectivity array.
        # Remember to shift indices as Python starts at 0
        iTriaStart = surfTriaPtr[iSurf]-1
        iTriaEnd = surfTriaPtr[iSurf+1]-1
        iQuadsStart = surfQuadsPtr[iSurf]-1
        iQuadsEnd = surfQuadsPtr[iSurf+1]-1

        # Slice connectivities
        currTriaConn = triaConn[:,iTriaStart:iTriaEnd]
        currQuadsConn = quadsConn[:,iQuadsStart:iQuadsEnd]

        # Initialize surface dictionary
        currSurf = {'triaConn':currTriaConn,
                    'quadsConn':currQuadsConn}

        # Add this section to the dictionary
        sectionDict[surf] = currSurf

        # Increment surface counter
        iSurf = iSurf + 1

    iCurve = 0
    for curve in curveNames:

        # Get indices to slice connectivity array.
        # Remember to shift indices as Python starts at 0
        iBarsStart = curveBarsPtr[iCurve]-1
        iBarsEnd = curveBarsPtr[iCurve+1]-1

        # Slice barsConn to get current curve elements
        slicedConn = barsConn[:,iBarsStart:iBarsEnd]

        # Initialize curve dictionary
        currCurve = {'barsConn':slicedConn}

        # Add this entry to the dictionary
        sectionDict[curve] = currCurve

        # Increment curve counter
        iCurve = iCurve + 1

    # Return sections dictionary
    return coor, sectionDict

#========================================


#=================================================================
# AUXILIARY COMPONENT FUNCTIONS
#=================================================================

def initialize_surface(TSurfComponent):

    '''
    This function receives a TSurf component and initializes its surface based
    on the initial connectivity and coordinates
    '''

    print 'My ID is',TSurfComponent.name

    # Now call general function that sets ADT
    update_surface(TSurfComponent)

#=================================================================

def update_surface(TSurfComponent):

    '''
    This function receives a TSurf component and initializes its surface
    ADT based on the current connectivity and coordinates.
    '''

    # Deallocate previous tree
    adtAPI.adtapi.adtdeallocateadts(TSurfComponent.name)

    # Set bounding box for new tree
    BBox = np.zeros((3, 2))
    useBBox = False

    # Compute set of nodal normals by taking the average normal of all
    # elements surrounding the node. This allows the meshing algorithms,
    # for instance, to march in an average direction near kinks.
    TSurfComponent.nodal_normals = adtAPI.adtapi.computenodalnormals(TSurfComponent.coor,
                                                                   TSurfComponent.triaConn,
                                                                   TSurfComponent.quadsConn)

    # Create new tree (the tree itself is stored in Fortran level)
    adtAPI.adtapi.adtbuildsurfaceadt(TSurfComponent.coor,
                                     TSurfComponent.triaConn, TSurfComponent.quadsConn,
                                     BBox, useBBox,
                                     TSurfComponent.comm.py2f(),
                                     TSurfComponent.name)

#=================================================================

def initialize_curves(TSurfComponent, sectionDict, selectedSections):

    '''
    This function initializes all curves given in sectionDict that are
    shown in selectedSections.

    INPUTS:
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                 keys are section names, while the values are also dictionaries with the
                 following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                 for curve sections.

    OUTPUTS:
    This function has no explicit outputs. It assigns TSurfComponent.Curves

    Ney Secco 2016-08
    '''

    # CURVE OBJECTS
    # A DiscreteSurf component may have multiple curve objects

    # Initialize dictionary that will hold all curve objects
    curveObjDict = {}

    # Now we will initialize curve objects
    for sectionName in sectionDict:

        if ('barsConn' in sectionDict[sectionName].keys()) and \
           (sectionName in selectedSections): # Then we have a curve section that should be stored

            # Get data from current section
            barsConn = sectionDict[sectionName]['barsConn']

            # Create Curve object and append entry to the dictionary
            curveObjDict[sectionName] = tsurf_component.TSurfCurve(TSurfComponent.coor, barsConn, sectionName)

    # Assign curve objects to ADT component
    TSurfComponent.Curves = curveObjDict

#=================================================================

def extract_curves_from_surface(TSurfComponent, feature='sharpness'):

    '''
    This function will extract features from the surface of TSurfComponent and
    return Curve objects.

    This function is intended for pre-processing only, and should not be
    called during an optimization as it might be slow.

    Ney Secco 2016-08
    '''

    # Get points and connectivities
    coor = TSurfComponent.coor
    triaConn = TSurfComponent.triaConn
    quadsConn = TSurfComponent.quadsConn

    # Get number of elements
    nTria = triaConn.shape[1]
    nQuads = quadsConn.shape[1]

    # BUILDING EDGE-TO-ELEMENT CONNECTIVITY

    # We need to identify the neighbors of each element, so we can compute edge information

    # Initialize array that will hold elements that share an edge.
    # This array will have an initial size, and we will resize it
    # if needed. Each column of the array will hold 4 elements, which will
    # identify the edge and the two elements that share it:
    #
    # For instance, assume that we have a bar that connects nodes 14 and 52,
    # And the elements that share this bar are quad 190 and tria 36. Then the entry
    # associated with this bar will be:
    #
    # sharedBarInfo[:,ii] = [14 (start point ID),
    #                        52 (end point ID),
    #                        -190 (first element ID, negative because it is a quad),
    #                        36 (second element ID)]

    initSize = 10000
    sharedBarInfo = np.zeros((4,initSize),order='F')

    # Create function to reallocate sharedBarInfo if necessary
    def reallocate_sharedBarInfo(sharedBarInfo):

        # Get old size
        oldSize = sharedBarInfo.shape[1]

        # Allocate new array
        newSharedBarInfo = np.zeros((4,oldSize+initSize),order='F',dtype=int)

        # Copy values from the old array
        newSharedBarInfo[:,:oldSize] = sharedBarInfo

        # Overwrite the old array
        sharedBarInfo = newSharedBarInfo

        # return new array
        return sharedBarInfo

    # Initially, we will create a dictionary whose keys represent edges and the values will
    # be the two elements that share this edge. The key will always be sorted, so that
    # The edges from the two elements will point to the same place in the dictionary.
    # For instance, if we have triangle 43 with node connectivity [3,56,7], and quad
    # 146 with connectivity [7,56,101,9], then the dictionary will get an entry
    # that looks: {[7,65]:[43,-146]} (Edge 7,56 is shared between tria 43 and quad 146).
    # Quads will be flagged by negative values.
    # An edge that is connected to a single element will have 0 in its values
    #
    # We will pop an entry from the dictionary as soon as it is complete, to avoid getting
    # a big dictionary. In other words, as soon as we detect both elements that share an
    # edge, we transfer that info from edge2Elem to sharedBarInfo (which is preallocated),
    # and then pop the entry from edge2Elem (which can't be preallocated).

    # Initialize dictionary
    edge2Elem = {}

    # Initialize edge counter
    nEdge = 0

    # Loop over all trias (I start at 1 so I could use 0 if no element is detected)
    for triaID in range(1,nTria+1):

        if np.mod(triaID,200) == 0:
            print 'Checking tria',triaID,'of',nTria

        # Get current element connectivity
        currConn = triaConn[:,triaID-1]

        # Determine the bars that define the element.
        # We sort the node order so that they will point to the same
        # place in the dictionary.
        bars = [tuple(np.sort([currConn[0], currConn[1]])),
                tuple(np.sort([currConn[1], currConn[2]])),
                tuple(np.sort([currConn[2], currConn[0]]))]

        # Check each bar
        for bar in bars:

            # Check if the bar is not defined in the dictionary yet
            if bar not in edge2Elem.keys():

                # Then add a new entry for this bar, containing
                # the first element we found
                edge2Elem[bar] = [triaID, 0]

            else:

                # The entry already exists. This means we already found
                # the first element that shares this edge and we just need
                # to add the second one.
                edge2Elem[bar][1] = triaID

                # Now that the information is complete, we can transfer it to the pre-allocated array.
                # But we need to check if we should increase the array size.

                if nEdge == sharedBarInfo.shape[1]:
                    sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

                # Assign new data
                sharedBarInfo[:,nEdge] = [bar[0],bar[1],edge2Elem[bar][0],edge2Elem[bar][1]]

                # Increment number of edges detected
                nEdge = nEdge + 1

                # Pop entry from dictionary
                edge2Elem.pop(bar)

    # Loop over all quads (I start at 1 so I could use 0 if no element is detected)
    for quadID in range(1,nQuads+1):

        if np.mod(quadID,200) == 0:
            print 'Checking quad',quadID,'of',nQuads

        # Get current element connectivity
        currConn = quadsConn[:,quadID-1]

        # Determine the bars that define the element.
        # We sort the node order so that they will point to the same
        # place in the dictionary.
        # We need to use tuples so that the connectivity can be used as
        # a dictionary key. Lists are not accepted.
        bars = [tuple(np.sort([currConn[0], currConn[1]])),
                tuple(np.sort([currConn[1], currConn[2]])),
                tuple(np.sort([currConn[2], currConn[3]])),
                tuple(np.sort([currConn[3], currConn[0]]))]

        # Check each bar
        for bar in bars:

            # Check if the bar is not defined in the dictionary yet
            if bar not in edge2Elem.keys():

                # Then add a new entry for this bar, containing
                # the first element we found.
                # Remember that quads are flagged with negative IDS.
                edge2Elem[bar] = [-quadID, 0]

            else:

                # The entry already exists. This means we already found
                # the first element that shares this edge and we just need
                # to add the second one.
                # Remember that quads are flagged with negative IDS.
                edge2Elem[bar][1] = -quadID

                # Now that the information is complete, we can transfer it to the pre-allocated array.
                # But we need to check if we should increase the array size.
                if nEdge == sharedBarInfo.shape[1]:
                    sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

                # Assign new data
                sharedBarInfo[:,nEdge] = [bar[0],bar[1],edge2Elem[bar][0],edge2Elem[bar][1]]

                # Increment number of edges detected
                nEdge = nEdge + 1

                # Pop entry from dictionary
                edge2Elem.pop(bar)

    # Now we transfer the remaining edged in edge2Elem to sharedBarInfo. These remaining bars
    # are at domain boundaries, so they are linked to a single element only
    for bar in edge2Elem.keys():

        # Check if we should increase the array size.
        if nEdge == sharedBarInfo.shape[1]:
            sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

        # Copy bar to sharedBarInfo
        sharedBarInfo[:, nEdge] = [bar[0],bar[1],edge2Elem[bar][0],edge2Elem[bar][1]]

        # Increment number of edges detected
        nEdge = nEdge + 1

        # Pop entry from dictionary
        edge2Elem.pop(bar)

    # Now we crop sharedBarInfo to have the exact number of edges
    sharedBarInfo = sharedBarInfo[:,:nEdge]


    # CURVE DETECTION

    # Phew! Now sharedBarInfo has all the information we need! It teels which elements
    # share a given edge. Just to recap, here's the sharedBarInfo structure:
    #
    # For instance, assume that we have a bar that connects nodes 14 and 52,
    # And the elements that share this bar are quad 190 and tria 36. Then the entry
    # associated with this bar will be:
    #
    # sharedBarInfo[:,ii] = [14 (start point ID),
    #                        52 (end point ID),
    #                        -190 (first element ID, negative because it is a quad),
    #                        36 (second element ID)]

    # Initialize list of bars that satisfy the selection criteria
    selectedBarsConn = []

    # Loop over the bars to extract features

    for barID in range(sharedBarInfo.shape[1]):

        # Gather data from sharedBarInfo
        node1 = sharedBarInfo[0,barID]
        node2 = sharedBarInfo[1,barID]
        element1 = sharedBarInfo[2,barID]
        element2 = sharedBarInfo[3,barID]

        # Run function to detect a given feature
        featureIsPresent = detect_feature(node1, node2, element1, element2,
                                          coor, triaConn, quadsConn,
                                          feature)

        # Store bar if feature is present
        if featureIsPresent:
            selectedBarsConn.append([node1,node2])

    # Run sorting algorithm on the detected bars
    # Remember that selectedBarsConn now will be a list containing
    # connectivities of multiple disconnect curves
    selectedBarsConn = FEsort(selectedBarsConn)

    # GENERATE CURVE OBJECTS

    # Initialize list of curve objects
    featureCurves = []

    # Initialize curve counter
    curveID = -1

    # Now we initiliaze curve objects for each disconnect curve we detected
    for currCurveConn in selectedBarsConn:

        # Increment curve counter
        curveID = curveID + 1

        # Create name for this curve
        curveName = 'extracted_' + feature + '_%03d'%curveID

        # Create new curve object
        newCurve = tsurf_component.TSurfCurve(coor, currCurveConn, curveName)

        # Initialize curve object and append it to the list
        TSurfComponent.add_curve(curveName, newCurve)

    # Print log
    print 'Number of extracted curves: ',len(selectedBarsConn)

#=================================================================

#=================================================================
# AUXILIARY SURFACE FUNCTIONS
#=================================================================

def merge_surface_sections(sectionDict, selectedSections):

    '''
    This function merges the connectivity data of all surface sections
    specified in selectedSections.

    INPUTS:
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                  keys are section names, while the values are also dictionaries with the
                  following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                  for curve sections.
    '''

    # Initialize global connectivities
    triaConn = None
    quadsConn = None

    # Loop over the selected surfaces to gather connectivities
    for sectionName in selectedSections:

        if sectionName not in sectionDict.keys():

            # User wants a section that is not defined. Print a warning message:
            print 'ERROR: Surface',sectionName,'is not defined. Check surface names in your CGNS file.'
            quit()

        elif 'triaConn' in sectionDict[sectionName].keys(): # Then we have a surface section

            # Assign connectivities
            if triaConn is None:
                # Start new connectivities if we have none
                triaConn = sectionDict[sectionName]['triaConn']
                quadsConn = sectionDict[sectionName]['quadsConn']
            else:
                # Append new connectivities
                triaConn = np.hstack([triaConn, sectionDict[sectionName]['triaConn']])
                quadsConn = np.hstack([quadsConn, sectionDict[sectionName]['quadsConn']])

        else:

            # The user provided a name that is not a surface section
            print sectionName,'is not a surface section.'

    return triaConn, quadsConn

#=================================================================

#=================================================================
# AUXILIARY CURVE FUNCTIONS
#=================================================================

#=================================================================

def split_curves(curveDict, criteria='sharpness'):

    '''
    This function will loop over all curves in curveDict and split
    them based on a given criteria. The available criteria are:
    criteria=['sharpness']

    Ney Secco 2016-08
    '''

    # Loop over every curve to check for splits
    for curveName in curveDict.keys():

        # First we remove the curve component from the dictionary
        curve = curveDict.pop(curveName)

        # Now we run the split function for this single curve
        splitCurves = split_curve_single(curve, curveName, criteria)

        # Now we add the split curves to the original curves dictionary
        curveDict.update(splitCurves)

#=============================================================

def split_curve_single(curve, curveName, criteria):

    '''
    This function receives a single curve object, splits it according to a criteria,
    and then return a new dictionary with split curves.

    ATTENTION: This function assumes that the FE data is sorted.

    INPUTS:
    curve: Curve object -> Curve that will be split

    curveName: string -> Name of the original curve. This name will be used
               to name the split curves.

    criteria: string -> Criteria that will be used to split curves. The options
              available for now are: ['sharpness']

    OUTPUTS:
    splitCurveDict: dictionary[curve objects] -> Dictionary containing split curves.

    Ney Secco 2016-08
    '''

    # READING INPUTS

    # Get coordinates and connectivities
    coor = curve.coor
    barsConn = curve.barsConn

    # Get the number of elements
    nElem = barsConn.shape[1]

    # DETECT KINKS

    # Initialize list of break elements. The break is considered to be
    # at the first point of the element
    breakList = []

    # The next step will depend on the criteria

    if criteria == 'sharpness':

        # This criteria will split the curve if its sharpness (defined here
        # as the change in tangential direction) is beyond a given threshold.

        # Get the tangent direction of the first bar element
        prevTan = coor[:,barsConn[1,0]-1] - coor[:,barsConn[0,0]-1]
        prevTan = prevTan/np.linalg.norm(prevTan)

        # Loop over the remaining bars to find sharp kinks
        for elemID in range(1,nElem):

            # Compute tangent direction of the current element
            currTan = coor[:,barsConn[1,elemID]-1] - coor[:,barsConn[0,elemID]-1]
            currTan = currTan/np.linalg.norm(currTan)

            # Compute change in direction between consecutive tangents
            angle = np.arccos(prevTan.dot(currTan))

            # Check if the angle is beyond a certain threshold
            if angle > 60*np.pi/180.0:

                # Store the current element as a break position
                breakList.append(elemID)

            # Now the current tangent will become the previous tangent of
            # the next iteration
            prevTan = currTan.copy()

    # CHECK IF BREAKS WERE DETECTED
    # We can stop this function earlier if we detect no break points
    if len(breakList) == 0:
        
        # In this case, we just create a dictionary with the original curve
        
        # Initialize dictionary that will contain the split curves
        splitCurvesDict = {}

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # Slice the original connectivity matrix
        splitBarsConn = barsConn[:,:]

        # Generate a name for this new curve
        splitCurveName = curveName

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitCurvesDict[splitCurveName] = splitCurve

        # Return dictionary with the single curve
        return splitCurvesDict


    # CREATE SPLIT CURVE OBJECTS

    # Now that we have the list of split points (given in splitList), we can create multiple
    # curve objects.

    # Count the number of curves that should be between the first and last break point
    nInnerCurves = len(breakList)-1

    # Initialize dictionary that will contain the split curves
    splitCurvesDict = {}

    # First we will define all curves that are between the first and the last
    # break point. The remaining part of the curve will be determined depending
    # if the original curve is periodic or not.
    for splitID in range(nInnerCurves):

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # Slice the original connectivity matrix
        splitBarsConn = barsConn[:,breakList[splitID]:breakList[splitID+1]]

        # Generate a name for this new curve
        splitCurveName = curveName + '_' + '%02d'%(splitID+1)

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # TODO: Add code to remove unused points from coor

        # Append new curve object to the dictionary
        splitCurvesDict[splitCurveName] = splitCurve

    # Now we need to create curves with elements that come before the first break point and after
    # the last break point.
    # We need to treat periodic curves differently. We use the same check to determine
    # the number of split curves. The periodic case will have minus one curve compared
    # to a non-periodic one.
    if barsConn[0,0] == barsConn[1,-1]: # Curve is periodic, so we need to define one curve

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = np.hstack([barsConn[:,breakList[-1]:],
                                   barsConn[:,:breakList[0]]])

        # Generate a name for this new curve
        splitCurveName = curveName + '_' + '%02d'%0

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # TODO: Add code to remove unused points from coor

        # Append new curve object to the dictionary
        splitCurvesDict[splitCurveName] = splitCurve

    else: # Curve is not periodic, so we need to define two curves

        # CURVE 0 : before the first break point

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = barsConn[:,:breakList[0]]

        # Generate a name for this new curve
        splitCurveName = curveName + '_' + '%02d'%0

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # TODO: Add code to remove unused points from coor

        # Append new curve object to the dictionary
        splitCurvesDict[splitCurveName] = splitCurve

        # CURVE 1 : after the first break point

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = barsConn[:,breakList[-1]:]

        # Generate a name for this new curve
        splitCurveName = curveName + '_' + '%02d'%(nInnerCurves+1)

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # TODO: Add code to remove unused points from coor

        # Append new curve object to the dictionary
        splitCurvesDict[splitCurveName] = splitCurve

    # Return the dictionary of new curves
    return splitCurvesDict

#=================================================================

#===================================
# INTERSECTION FUNCTIONS
#===================================

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

            # Compute new intersections for the current pair
            newIntersections = _compute_pair_intersection(TSurfComponentList[ii],
                                                         TSurfComponentList[jj],
                                                         distTol)

            # Append new curve objects to the list
            Intersections = Intersections + newIntersections

    # Print log
    print 'Computed',len(Intersections),'intersection curves.'

    # Return all intersections
    return Intersections

#=================================================================

def _compute_pair_intersection(TSurfComponentA, TSurfComponentB, distTol):

    '''
    This function finds intersection curves between components A and B
    and returns a list of curve objects corresponding to these intersections.
    '''

    # Call Fortran code to find intersections
    intersectionAPI.intersectionapi.computeintersection(TSurfComponentA.coor,
                                                        TSurfComponentA.triaConn,
                                                        TSurfComponentA.quadsConn,
                                                        TSurfComponentB.coor,
                                                        TSurfComponentB.triaConn,
                                                        TSurfComponentB.quadsConn,
                                                        distTol)

    # Retrieve results from Fortran
    coor = np.array(intersectionAPI.intersectionapi.coor)
    barsConn = np.array(intersectionAPI.intersectionapi.barsconn)

    # Release memory used by Fortran so we run another intersection in the future
    intersectionAPI.intersectionapi.releasememory()

    # Initialize list to hold all intersection curves
    Intersections = []

    if len(coor) > 0:

        barsConn = np.sort(barsConn, axis=0)
        barsConn = np.vstack({tuple(col) for col in barsConn.T}).T

        # Sort FE data. After this step, barsConn may become a list if we
        # have two disconnect intersection curves.
        barsConn = FEsort(barsConn.T.tolist())

        # barsConn is a list of curves. Each element of the list brings
        # FE data for a continuous curve. If the intersection generates
        # multiple disconnect curves, barsConn will be a list with
        # multiple elements.
        # For each curve (and therefore for each barsConn element), we
        # will initialize a Curve object to hold intersection FE data.

        # Set intersection counter
        intCounter = -1

        for currConn in barsConn:

            intCounter = intCounter + 1
            
            # Gather name of parent components
            name1 = TSurfComponentA.name
            name2 = TSurfComponentB.name
            
            # Create a name for the curve
            curveName = 'int_'+name1+'_'+name2+'_%02d'%intCounter
            
            # Create new curve object
            newCurve = tsurf_component.TSurfCurve(coor, currConn, curveName)

            # Initialize curve object and append it to the list
            Intersections.append(newCurve)

    # Return intersection FE data
    return Intersections


#=================================================================

#=================================================================

#=================================================================

#===================================
# GENERAL AUXILIARY FUNCTIONS
#===================================

def formatStringArray(fortranArray):

    '''
    This functions receives a numpy array of characters coming from Fortran
    and returns a list of formatted data.
    '''

    # Get array size
    shape = fortranArray.shape

    # We need to tranpose the array considering the Fortran ordering
    fortranArray = fortranArray.reshape([shape[1], shape[0]], order='F')

    # Initialize new array
    pythonArray = []

    # Now we format and concatenate the elements
    for index in range(shape[0]):

        # format string
        newString = ''.join(fortranArray[:,index]).strip().lower()

        # Add formatted string to the list
        pythonArray.append(newString)

    # Return formatted array
    return pythonArray

#=============================================================

def FEsort(barsConn):

    '''
    This function can be used to sort connectivities coming from the CGNS file
    It assumes that we only have one curve. It will crash if we have
    multiple curves defined in the same FE set.

    barsConnIn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

    newConn will be set to None if the sorting algorithm fails.
    '''

    # We will solve this in an iterative process until the connectivities do not change
    # anymore

    # Initialize the newest connectivities as our current input.
    newConn = barsConn

    # Say that we still need to search
    keep_searching = True

    while keep_searching:

        # If nothing happens, we will get out of the loop
        keep_searching = False

        # Update "old" connectivities
        oldConn = newConn[:]

        # Initialize list of new connectivities with the first element of the old ones.
        # We will also pop this element from oldConn
        newConn = [oldConn.pop(0)]

        # We will keep sorting until all elements of oldConn are assigned and popped out

        while len(oldConn) > 0:

            # Pop another element from oldConn
            oldElement = oldConn.pop(0)

            # We do not know if this element could be linked beforehand
            linked_element = False

            # Now we check if we can fit this element into any new connectivity
            newElemCounter = 0
            for newElement in newConn:

                if oldElement[0] == newElement[0]:

                    # We will flip the old element and place it at the beginning of the
                    # new connectivity
                    newConn[newElemCounter] = oldElement[::-1] + newElement[1:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[0] == newElement[-1]:

                    # Just append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[0]:

                    # Place the old element at the beginning
                    newConn[newElemCounter] = oldElement + newElement[1:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[-1]:

                    # Flip and append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement[::-1]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                # Increment element counter
                newElemCounter = newElemCounter + 1

            if not linked_element: # We have an element that is not connected to anyone

                    # Define a new curve in newConn
                    newConn.append(oldElement)

    # Right now, newConn represent each line by an 1D array of point IDs.
    # e.g. [[2,3,4,5],[1,7,8]]
    # We need to convert this back to FE format, represented by 2D arrays.
    # e.g. [[[2,3,4],[3,4,5]],[[1,7],[7,8]]]

    # Initialize FE array for each curve
    newConnFE = [np.zeros((2,len(curve)-1),dtype=int) for curve in newConn]

    # Fill FE data for each curve
    for curveID in range(len(newConn)):

        # Get current curve
        curve = newConn[curveID]

        # Get FE list of the current curve
        FEcurve = newConnFE[curveID]

        # Assign FE connectivity
        for pointID in range(1,len(curve)):

            # Get point indices
            prevPoint = curve[pointID-1]
            currPoint = curve[pointID]

            # Assign bar FE
            FEcurve[0,pointID-1] = prevPoint
            FEcurve[1,pointID-1] = currPoint


    # Now we do a final check to remove degenerate bars (bars that begin and end at the same point)
    # Take every disconnect curve in newConnFE:
    for curveID in range(len(newConnFE)):

        # We convert the connectivity array to list so we can 'pop' elements
        curveFE = newConnFE[curveID].T.tolist()

        for FE in curveFE:

            # Check if the start and end points are the same
            if FE[0] == FE[1]:

                # Remove FE
                curveFE.remove(FE)

        # Convert connectivity back to numpy array
        newConnFE[curveID] = np.array(curveFE, order='F').T

    # Return the sorted array
    return newConnFE

#=============================================================

def remove_unused_points(coor,
                         triaConn=np.zeros((0,0)),
                         quadsConn=np.zeros((0,0)),
                         barsConn=np.zeros((0,0))):

    '''
    This function will removed unused points from coor and
    also update all connectivities.

    ATTENTION:
    All inputs must be arrays.

    OUTPUTS:
    This function returns cropCoor, which is the cropped set of nodes.
    This function also updates triaConn, quadsConn and barsConn.
    '''

    # Get total number of points and elements
    nPoints = coor.shape[1]
    nTria = triaConn.shape[1]
    nQuads = quadsConn.shape[1]
    nBars = barsConn.shape[1]

    # Initialize mask to identify used points
    usedPtsMask = np.zeros(nPoints, dtype=int)

    # First we loop over all elements to create a mask that indicates used points
    for triaID in range(nTria):
        # Flag used points
        usedPtsMask[triaConn[0,triaID]-1] = 1
        usedPtsMask[triaConn[1,triaID]-1] = 1
        usedPtsMask[triaConn[2,triaID]-1] = 1

    for quadID in range(nQuads):
        # Flag used points
        usedPtsMask[quadsConn[0,quadID]-1] = 1
        usedPtsMask[quadsConn[1,quadID]-1] = 1
        usedPtsMask[quadsConn[2,quadID]-1] = 1
        usedPtsMask[quadsConn[3,quadID]-1] = 1

    for barID in range(nBars):
        # Flag used points
        usedPtsMask[barsConn[0,barID]-1] = 1
        usedPtsMask[barsConn[1,barID]-1] = 1

    # Now we can compute the number of points actually used
    nUsedPts = np.sum(usedPtsMask)

    # Initialize new coordinate array
    cropCoor = np.zeros((3,nUsedPts),order='F')

    # Initialize counter to fill cropCoor
    cropPointID = -1

    # Now we fill the points of the cropped array
    for pointID in range(nPoints):

        # Check if the point is used
        if usedPtsMask[pointID] == 1:

            # Increment counter
            cropPointID = cropPointID + 1

            # Add point to the cropped array
            cropCoor[:,cropPointID] = coor[:,pointID]

            # Now we replace the value in the mask array so that we
            # can use it as a pointer from coor to cropCoor when we
            # update element connectivities.
            # The +1 is necessary because Fortran use 1-based indexing
            usedPtsMask[pointID] = cropPointID + 1

    # Now we need to update connectivities so that they point to the the correct
    # indices of the cropped array
    for triaID in range(nTria):
        # Use pointer to update connectivity
        triaConn[0,triaID] = usedPtsMask[triaConn[0,triaID]-1]
        triaConn[1,triaID] = usedPtsMask[triaConn[1,triaID]-1]
        triaConn[2,triaID] = usedPtsMask[triaConn[2,triaID]-1]

    for quadID in range(nQuads):
        # Use pointer to update connectivity
        quadsConn[0,quadID] = usedPtsMask[quadsConn[0,quadID]-1]
        quadsConn[1,quadID] = usedPtsMask[quadsConn[1,quadID]-1]
        quadsConn[2,quadID] = usedPtsMask[quadsConn[2,quadID]-1]
        quadsConn[3,quadID] = usedPtsMask[quadsConn[3,quadID]-1]

    for barID in range(nBars):
        # Use pointer to update connectivity
        barsConn[0,barID] = usedPtsMask[barsConn[0,barID]-1]
        barsConn[1,barID] = usedPtsMask[barsConn[1,barID]-1]

    # Return the new set of nodes
    return cropCoor

#=============================================================

def detect_feature(node1, node2, element1, element2,
                   coor, triaConn, quadsConn,
                   feature):

    '''
    This function checks if the bar that connects node1 and node2, and
    is shared between element1 and element2 has a desired feature. This is
    used by extract_curves_from_surface function.

    INPUTS:

    node1: integer -> node1 index in coor (1-based since this comes from Fortran)

    node2: integer -> node2 index in coor (1-based since this comes from Fortran)

    element1: integer -> element index in triaConn (if positive) or quadsConn
              (if negative). (1-based since this comes from Fortran)

    element2: integer -> element index in triaConn (if positive) or quadsConn
              (if negative). (1-based since this comes from Fortran)

    coor: float[3,nNodes] -> nodal coordinates

    triaConn: integer[3,nTria] -> triangle connectivities. (1-based since this comes from Fortran)

    quadsConn: integer[4,nQuads] -> quad connectivities. (1-based since this comes from Fortran)

    feature: string -> feature that should be detected. Available options are:
             ['sharpness']

    OUTPUTS:

    featureIsDetected: logical -> True if feature is detected in this edge. Otherwise it is False.

    Ney Secco 2016-08
    '''

    # Check the criteria

    if feature == 'sharpness':

        # We need to compute the angle between the normals that share the bar

        # Let's get two inplane vectors for each element first:
        # v11 - first vector of the first element
        # v21 - second vector of the first element
        # v12 - first vector of the second element
        # v22 - second vector of the second element

        # We start with element 2 because if it is zero, then the bar is at the
        # border of the domain and do not has a second element. So we cannot use
        # it for sharpness detection.

        # ELEMENT 2
        if element2 == 0:
            featureIsDetected = False
            return featureIsDetected

        elif element2 > 0: # We have a tria

            # Get element nodes (The -1 is due 1-based indexing in Fortran)
            node1Coor = coor[:,triaConn[0,element2-1]-1]
            node2Coor = coor[:,triaConn[1,element2-1]-1]
            node3Coor = coor[:,triaConn[2,element2-1]-1]

            # Get inplane vectors so that the normal points outside
            v12 = node2Coor - node1Coor
            v22 = node3Coor - node1Coor

        elif element2 < 0: # We have a quad

            # Get element nodes (The -1 is due 1-based indexing in Fortran)
            node1Coor = coor[:,quadsConn[0,-element2-1]-1]
            node2Coor = coor[:,quadsConn[1,-element2-1]-1]
            node3Coor = coor[:,quadsConn[2,-element2-1]-1]
            node4Coor = coor[:,quadsConn[3,-element2-1]-1]

            # Get inplane vectors so that the normal points outside
            v12 = node3Coor - node1Coor
            v22 = node4Coor - node2Coor

        # ELEMENT 1
        if element1 > 0: # We have a tria

            # Get element nodes (The -1 is due 1-based indexing in Fortran)
            node1Coor = coor[:,triaConn[0,element1-1]-1]
            node2Coor = coor[:,triaConn[1,element1-1]-1]
            node3Coor = coor[:,triaConn[2,element1-1]-1]

            # Get inplane vectors so that the normal points outside
            v11 = node2Coor - node1Coor
            v21 = node3Coor - node1Coor

        elif element1 < 0: # We have a quad

            # Get element nodes (The -1 is due 1-based indexing in Fortran)
            node1Coor = coor[:,quadsConn[0,-element1-1]-1]
            node2Coor = coor[:,quadsConn[1,-element1-1]-1]
            node3Coor = coor[:,quadsConn[2,-element1-1]-1]
            node4Coor = coor[:,quadsConn[3,-element1-1]-1]

            # Get inplane vectors so that the normal points outside
            v11 = node3Coor - node1Coor
            v21 = node4Coor - node2Coor

        # Take the cross products to get normals and normalize the normals =P
        # Element 1
        n1 = np.cross(v11,v21)
        n1 = n1/np.linalg.norm(n1)
        # Element 2
        n2 = np.cross(v12,v22)
        n2 = n2/np.linalg.norm(n2)

        # Use dot product to find the angle between the normals
        angle = np.arccos(n1.dot(n2))

        # We have a "sharp" edge if this angle is beyond a certaing threshold
        if angle > 60*np.pi/180:
            featureIsDetected = True
            return featureIsDetected
        else:
            featureIsDetected = False
            return featureIsDetected

    elif feature == 'open_ends':

        # In this case we check if an edge is shared by a single element only,
        # as this characterizes an open end.

        # We only need to check element 2 because if it is zero, then the bar is at the
        # border of the domain.

        # ELEMENT 2
        if element2 == 0:
            featureIsDetected = True
            return featureIsDetected
        else:
            featureIsDetected = False
            return featureIsDetected            

    else:

        print 'ERROR: Feature',feature,'cannot be detected as it is not an option.'
        quit()



#============================================================

def print_arrows(filename,coor,barsConn):

    variable_names = ['CoordinateX','CoordinateY','CoordinateZ','dX','dY','dZ']

    title = 'FE_orientation'

    data_points = []

    for bar in barsConn.T.tolist():
        X0 = coor[0,bar[0]-1]
        Y0 = coor[1,bar[0]-1]
        Z0 = coor[2,bar[0]-1]

        X1 = coor[0,bar[1]-1]
        Y1 = coor[1,bar[1]-1]
        Z1 = coor[2,bar[1]-1]

        dX = X1-X0
        dY = Y1-Y0
        dZ = Z1-Z0

        data_points.append([X0,Y0,Z0,dX,dY,dZ])

    write_tecplot_scatter(filename,title,variable_names,data_points)

#============================================================

def write_tecplot_scatter(filename,title,variable_names,data_points):

    # Open the data file
    fid = open(filename,'w')
    
    # Write the title
    fid.write('title = '+title+'\n')

    # Write tha variable names
    varnames_commas = ','.join(variable_names) # Merge names in a single string separated by commas
    fid.write('variables = '+varnames_commas+',\n') # Write to the file

    # Write data points
    if type(data_points) is list: # Check if user provided a list
        for point in data_points:
            str_points = [str(x) for x in point] # Covert each entry to string
            str_points = ' '.join(str_points) # Merge all entries in a single string separated by whitespace
            fid.write(str_points+'\n') # Write to file
    else: # The user probably provided a numpy array
        for index in range(data_points.shape[0]):
            str_points = [str(x) for x in data_points[index,:]]
            str_points = ' '.join(str_points) # Merge all entries in a single string separated by whitespace
            fid.write(str_points+'\n') # Write to file
    
    # Close file
    fid.close()

#============================================================
#============================================================
#============================================================


def translate(Component, x, y, z):
    '''
    Translate coor based on given x, y, and z values.
    John Jasa 2016-08
    '''

    Component.coor = Component.coor[:, :] + np.atleast_2d(np.array([x, y, z])).T

def scale(Component, factor, point=None):
    '''
    Scale coor about a defined point or the center.
    John Jasa 2016-08
    '''

    coor = Component.coor
    if not point:
        point = np.mean(coor)
    relCoor = coor - point
    relCoor *= factor
    Component.coor = relCoor + point

def rotate(Component, angle, axis, point=None):
    '''
    Rotate coor about an axis at a defined angle (in degrees).
    John Jasa 2016-08
    '''


    coor = Component.coor
    if not point:
        point = np.mean(coor)
    angle = angle * np.pi / 180.
    rotationMat = np.zeros((3, 3))

    if axis==0:
        ind1 = 1
        ind2 = 2
    elif axis==1:
        ind1 = 2
        ind2 = 0
    else:
        ind1 = 0
        ind2 = 1

    rotationMat[ind1, ind1] = np.cos(angle)
    rotationMat[ind1, ind2] = -np.sin(angle)
    rotationMat[axis, axis] = 1
    rotationMat[ind2, ind1] = np.sin(angle)
    rotationMat[ind2, ind2] = np.cos(angle)

    Component.coor = np.einsum('ij, jk -> ik', rotationMat, coor-point) + point
