import numpy as np
from . import cgnsAPI, adtAPI, utilitiesAPI, intersectionAPI
from mpi4py import MPI
from . import tsurf_component

"""
TODO:
- Add variable angle threshold to curve feature extraction and curve split
"""

# ===========================================
# AUXILIARY READING FUNCTIONS
# ===========================================


def getCGNSsections(inputFile, comm=MPI.COMM_WORLD):

    """
    This function opens a CGNS file, reads its sections, and returns a
    dictionary with selected sections. This code also initializes ADT
    for each section.
    """

    ### Read CGNS file in Fortran

    # This will de done in two steps:
    # The first Fortran call will give us a tuple containing the number of nodes, trias, quads, etc...
    # These are the sizes of the allocatable arrays that were defined in the Fortran call.
    # The second Fortran call will use this information to actually retrieve the surface data.
    # We have to use this two steps approach so that Python does not need direct access to the
    # allocatable arrays, since this does not work when we use Intel compilers.

    # First Fortran call to read CGNS file and get array sizes
    arraySizes = cgnsAPI.cgnsapi.readcgns(inputFile, comm.py2f())

    # If an arraySizes value is zero, which mean we have empty allocated arrays, we need to
    # convert them to at least 1 so that f2py does not complain when we retrieve values
    # (It can't understand that we want to retrieve a 0-sized array). So we will do the following
    # trick:
    # If an arraySizes element is 0, we will set it to -1. We will call the Fortran code with the
    # absolute values of arraySizes (so -1 will be converted to 1 in this case) and we will get
    # a corresponding 1-sized array from Fortran, with meaningless values. Then, back in Python,
    # we will erase the data from the arrays that have corresponding -1 in arraySizes.

    # Conver arraySizes from tuple to list so we can modify it
    arraySizes = list(arraySizes)

    # Loop over arraySizes to replace 0 with -1
    for ii in range(len(arraySizes)):
        if arraySizes[ii] == 0:
            arraySizes[ii] = -1

    # Second Fortran call to retrieve data from the CGNS file.
    # We need to do actual copies, otherwise data will be overwritten if we read another CGNS file.
    # We subtract one to make indices consistent with the Python 0-based indices.
    # We also need to transpose it since Python and Fortran use different orderings to store matrices
    # in memory.
    # We cannot subtract 1 from triaConn and quadsConn because the ADT will have pointers
    # to these values in the Fortran level, so we have to keep these variables fixed, otherwise
    # we might lose the memory location. So remember that triaConnF and quadsConnF use Fortran
    # indices (starting at 1).
    cgnsArrays = cgnsAPI.cgnsapi.retrievedata(*np.abs(arraySizes))

    # Split results
    coor = np.array(cgnsArrays[0].T)
    triaConnF = np.array(cgnsArrays[1].T)
    quadsConnF = np.array(cgnsArrays[2].T)
    barsConn = np.array(cgnsArrays[3].T - 1)
    surfTriaPtr = np.array(cgnsArrays[4] - 1)
    surfQuadsPtr = np.array(cgnsArrays[5] - 1)
    curveBarsPtr = np.array(cgnsArrays[6] - 1)
    surfNames = cgnsArrays[7].copy()
    curveNames = cgnsArrays[8].copy()

    # Transform flagged arrays back to empty ones
    if arraySizes[1] == -1:
        triaConnF = np.zeros((0, 3))
    if arraySizes[2] == -1:
        quadsConnF = np.zeros((0, 4))
    if arraySizes[3] == -1:
        barsConnF = np.zeros((0, 2))
    if arraySizes[4] == -1:
        surfTriaPtr = np.zeros(0)
    if arraySizes[5] == -1:
        surfQuadsPtr = np.zeros(0)
    if arraySizes[6] == -1:
        curveBarsPtr = np.zeros(0)
    if arraySizes[7] == -1:
        surfNames = np.zeros((0, 32))
    if arraySizes[8] == -1:
        curveNames = np.zeros((0, 32))

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
        iTriaStart = surfTriaPtr[iSurf]
        iTriaEnd = surfTriaPtr[iSurf + 1]
        iQuadsStart = surfQuadsPtr[iSurf]
        iQuadsEnd = surfQuadsPtr[iSurf + 1]

        # Slice connectivities
        currTriaConnF = triaConnF[iTriaStart:iTriaEnd, :]
        currQuadsConnF = quadsConnF[iQuadsStart:iQuadsEnd, :]

        # Initialize surface dictionary
        currSurf = {"triaConnF": currTriaConnF, "quadsConnF": currQuadsConnF}

        # Add this section to the dictionary
        sectionDict[surf] = currSurf

        # Increment surface counter
        iSurf = iSurf + 1

    iCurve = 0
    for curve in curveNames:

        # Get indices to slice connectivity array.
        iBarsStart = curveBarsPtr[iCurve]
        iBarsEnd = curveBarsPtr[iCurve + 1]

        # Slice barsConn to get current curve elements
        slicedConn = barsConn[iBarsStart:iBarsEnd, :]

        # Initialize curve dictionary
        currCurve = {"barsConn": slicedConn}

        # Add this entry to the dictionary
        sectionDict[curve] = currCurve

        # Increment curve counter
        iCurve = iCurve + 1

    # Return sections dictionary
    return coor, sectionDict


# ========================================


def read_tecplot_curves(fileName, mergeTol=1e-7):

    """
    This function reads curve described in tecplot format.

    mergeTol: float -> Tolerance used to merge nearby nodes after reading tecplot data.

    Ney Secco 2016-11
    """

    # IMPORTS
    from ....utilities import plot3d_interface as p3d

    # Get curve data
    tecCurves = p3d.read_tecplot_curves(fileName)

    # Create curve objects and append them to a dictionary
    curves = {}

    # Now we convert every Tecplot data curve to a TSurf curve
    for curveID in range(len(tecCurves)):

        # Get current tecplot curve
        currTecCurve = tecCurves[curveID]

        # Create new curve object
        currCurve = tsurf_component.TSurfCurve(currTecCurve.coor, currTecCurve.barsConn, currTecCurve.name, mergeTol)

        # Add it to the TSurf curve dictionary
        curves[currTecCurve.name] = currCurve

    # Return dictionary with TSurf curves
    return curves


# =================================================================
# AUXILIARY COMPONENT FUNCTIONS
# =================================================================


def initialize_surface(TSurfGeometry):

    """
    This function receives a TSurf component and initializes its surface based
    on the initial connectivity and coordinates
    """

    print("My ID is", TSurfGeometry.name)

    # Now call general function that sets ADT
    update_surface(TSurfGeometry)


# =================================================================


def update_surface(TSurfGeometry):

    """
    This function receives a TSurf component and initializes its surface
    ADT based on the current connectivity and coordinates.
    """

    # Deallocate previous tree
    adtAPI.adtapi.adtdeallocateadts(TSurfGeometry.name)

    # Set bounding box for new tree
    BBox = np.zeros((2, 3))
    useBBox = False

    # Compute set of nodal normals by taking the average normal of all
    # elements surrounding the node. This allows the meshing algorithms,
    # for instance, to march in an average direction near kinks.
    nodal_normals = adtAPI.adtapi.adtcomputenodalnormals(
        TSurfGeometry.coor.T, TSurfGeometry.triaConnF.T, TSurfGeometry.quadsConnF.T
    )
    TSurfGeometry.nodal_normals = nodal_normals.T

    # Create new tree (the tree itself is stored in Fortran level)
    adtAPI.adtapi.adtbuildsurfaceadt(
        TSurfGeometry.coor.T,
        TSurfGeometry.triaConnF.T,
        TSurfGeometry.quadsConnF.T,
        BBox.T,
        useBBox,
        TSurfGeometry.comm.py2f(),
        TSurfGeometry.name,
    )


# =================================================================


def initialize_curves(TSurfGeometry, sectionDict, selectedSections):

    """
    This function initializes all curves given in sectionDict that are
    shown in selectedSections.

    Parameters
    ----------
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                 keys are section names, while the values are also dictionaries with the
                 following fields:'quadsConn','triaConn' for surface sections and 'barsConn'
                 for curve sections.

    Returns
    -------
    This function has no explicit outputs. It assigns TSurfGeometry.curves

    Ney Secco 2016-08
    """

    # CURVE OBJECTS
    # A DiscreteSurf component may have multiple curve objects

    # Initialize dictionary that will hold all curve objects
    curveObjDict = {}

    # Now we will initialize curve objects
    for sectionName in sectionDict:

        if ("barsConn" in list(sectionDict[sectionName].keys())) and (
            sectionName in selectedSections
        ):  # Then we have a curve section that should be stored

            # Get data from current section
            barsConn = sectionDict[sectionName]["barsConn"]

            # Create Curve object and append entry to the dictionary
            curveObjDict[sectionName] = tsurf_component.TSurfCurve(TSurfGeometry.coor, barsConn, sectionName)

    # Assign curve objects to ADT component
    TSurfGeometry.curves = curveObjDict


# =================================================================


def extract_curves_from_surface(TSurfGeometry, feature="sharpness"):

    """
    This function will extract features from the surface of TSurfGeometry and
    return Curve objects.

    This function is intended for pre-processing only, and should not be
    called during an optimization as it might be slow.

    Ney Secco 2016-08
    """

    # Get points and connectivities
    coor = TSurfGeometry.coor
    triaConnF = TSurfGeometry.triaConnF
    quadsConnF = TSurfGeometry.quadsConnF

    # Get number of elements
    nTria = triaConn.shape[0]
    nQuads = quadsConn.shape[0]

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
    # sharedBarInfo[ii,:] = [14 (start point ID),
    #                        52 (end point ID),
    #                        -190 (first element ID, negative because it is a quad),
    #                        36 (second element ID)]

    initSize = 10000
    sharedBarInfo = np.zeros((initSize, 4))

    # Create function to reallocate sharedBarInfo if necessary
    def reallocate_sharedBarInfo(sharedBarInfo):

        # Get old size
        oldSize = sharedBarInfo.shape[0]

        # Allocate new array
        newSharedBarInfo = np.zeros((oldSize + initSize, 4), dtype=int)

        # Copy values from the old array
        newSharedBarInfo[:oldSize, :] = sharedBarInfo

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
    # An edge that is connected to a single element will have None in its values
    #
    # We will pop an entry from the dictionary as soon as it is complete, to avoid getting
    # a big dictionary. In other words, as soon as we detect both elements that share an
    # edge, we transfer that info from edge2Elem to sharedBarInfo (which is preallocated),
    # and then pop the entry from edge2Elem (which can't be preallocated).

    # Initialize dictionary
    edge2Elem = {}

    # Initialize edge counter
    nEdge = 0

    # Loop over all trias
    for triaID in range(nTria):

        if np.mod(triaID, 200) == 0:
            print("Checking tria", triaID + 1, "of", nTria)

        # Get current element connectivity
        # The minus 1 is due to Fortran indexing present in triaConnF
        currConn = triaConnF[triaID, :] - 1

        # Determine the bars that define the element.
        # We sort the node order so that they will point to the same
        # place in the dictionary.
        bars = [
            tuple(np.sort([currConn[0], currConn[1]])),
            tuple(np.sort([currConn[1], currConn[2]])),
            tuple(np.sort([currConn[2], currConn[0]])),
        ]

        # Check each bar
        for bar in bars:

            # Check if the bar is not defined in the dictionary yet
            if bar not in list(edge2Elem.keys()):

                # Then add a new entry for this bar, containing
                # the first element we found
                edge2Elem[bar] = [triaID, None]

            else:

                # The entry already exists. This means we already found
                # the first element that shares this edge and we just need
                # to add the second one.
                edge2Elem[bar][1] = triaID

                # Now that the information is complete, we can transfer it to the pre-allocated array.
                # But we need to check if we should increase the array size.

                if nEdge == sharedBarInfo.shape[0]:
                    sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

                # Assign new data
                sharedBarInfo[nEdge, :] = [bar[0], bar[1], edge2Elem[bar][0], edge2Elem[bar][1]]

                # Increment number of edges detected
                nEdge = nEdge + 1

                # Pop entry from dictionary
                edge2Elem.pop(bar)

    # Loop over all quads
    for quadID in range(nQuads):

        if np.mod(quadID, 200) == 0:
            print("Checking quad", quadID + 1, "of", nQuads)

        # Get current element connectivity
        # The minus 1 is due to Fortran indexing present in quadsConnF
        currConn = quadsConnF[quadID, :] - 1

        # Determine the bars that define the element.
        # We sort the node order so that they will point to the same
        # place in the dictionary.
        # We need to use tuples so that the connectivity can be used as
        # a dictionary key. Lists are not accepted.
        bars = [
            tuple(np.sort([currConn[0], currConn[1]])),
            tuple(np.sort([currConn[1], currConn[2]])),
            tuple(np.sort([currConn[2], currConn[3]])),
            tuple(np.sort([currConn[3], currConn[0]])),
        ]

        # Check each bar
        for bar in bars:

            # Check if the bar is not defined in the dictionary yet
            if bar not in list(edge2Elem.keys()):

                # Then add a new entry for this bar, containing
                # the first element we found.
                # Remember that quads are flagged with negative IDS.
                edge2Elem[bar] = [-quadID, None]

            else:

                # The entry already exists. This means we already found
                # the first element that shares this edge and we just need
                # to add the second one.
                # Remember that quads are flagged with negative IDS.
                edge2Elem[bar][1] = -quadID

                # Now that the information is complete, we can transfer it to the pre-allocated array.
                # But we need to check if we should increase the array size.
                if nEdge == sharedBarInfo.shape[0]:
                    sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

                # Assign new data
                sharedBarInfo[nEdge, :] = [bar[0], bar[1], edge2Elem[bar][0], edge2Elem[bar][1]]

                # Increment number of edges detected
                nEdge = nEdge + 1

                # Pop entry from dictionary
                edge2Elem.pop(bar)

    # Now we transfer the remaining edges in edge2Elem to sharedBarInfo. These remaining bars
    # are at domain boundaries, so they are linked to a single element only
    for bar in list(edge2Elem.keys()):

        # Check if we should increase the array size.
        if nEdge == sharedBarInfo.shape[0]:
            sharedBarInfo = reallocate_sharedBarInfo(sharedBarInfo)

        # Copy bar to sharedBarInfo
        sharedBarInfo[nEdge, :] = [bar[0], bar[1], edge2Elem[bar][0], edge2Elem[bar][1]]

        # Increment number of edges detected
        nEdge = nEdge + 1

        # Pop entry from dictionary
        edge2Elem.pop(bar)

    # Now we crop sharedBarInfo to have the exact number of edges
    sharedBarInfo = sharedBarInfo[:nEdge, :]

    # CURVE DETECTION

    # Phew! Now sharedBarInfo has all the information we need! It tells which elements
    # share a given edge. Just to recap, here's the sharedBarInfo structure:
    #
    # For instance, assume that we have a bar that connects nodes 14 and 52,
    # And the elements that share this bar are quad 190 and tria 36. Then the entry
    # associated with this bar will be:
    #
    # sharedBarInfo[ii,:] = [14 (start point ID),
    #                        52 (end point ID),
    #                        -190 (first element ID, negative because it is a quad),
    #                        36 (second element ID)]

    # Initialize list of bars that satisfy the selection criteria
    selectedBarsConn = []

    # Loop over the bars to extract features

    for barID in range(sharedBarInfo.shape[0]):

        # Gather data from sharedBarInfo
        node1 = sharedBarInfo[barID, 0]
        node2 = sharedBarInfo[barID, 1]
        element1 = sharedBarInfo[barID, 2]
        element2 = sharedBarInfo[barID, 3]

        # Run function to detect a given feature
        featureIsPresent = detect_feature(node1, node2, element1, element2, coor, triaConnF, quadsConnF, feature)

        # Store bar if feature is present
        if featureIsPresent:
            selectedBarsConn.append([node1, node2])

    # Run sorting algorithm on the detected bars
    # Remember that selectedBarsConn now will be a list containing
    # connectivities of multiple disconnect curves
    selectedBarsConn, dummy_map = FEsort(selectedBarsConn)

    # GENERATE CURVE OBJECTS

    # Initialize list of curve objects
    featurecurves = []

    # Initialize curve counter
    curveID = -1

    # Now we initiliaze curve objects for each disconnect curve we detected
    for currCurveConn in selectedBarsConn:

        # Increment curve counter
        curveID = curveID + 1

        # Create name for this curve
        curveName = "extracted_" + feature + "_%03d" % curveID

        # Create new curve object
        newCurve = tsurf_component.TSurfCurve(coor, currCurveConn, curveName)

        # Initialize curve object and append it to the list
        TSurfGeometry.add_curve(newCurve)

    # Print log
    print("Number of extracted curves: ", len(selectedBarsConn))


# =================================================================

# =================================================================
# AUXILIARY SURFACE FUNCTIONS
# =================================================================

# =================================================================


def merge_surface_sections(sectionDict, selectedSections):

    """
    This function merges the connectivity data of all surface sections
    specified in selectedSections.

    Parameters
    ----------
    sectionDict: dictionary{sectionName,sDict} -> Dictionary that contains surface info. The
                  keys are section names, while the values are also dictionaries with the
                  following fields:'quadsConnF','triaConnF' for surface sections and 'barsConn'
                  for curve sections.
    """

    # Initialize global connectivities
    triaConnF = None
    quadsConnF = None

    # Loop over the selected surfaces to gather connectivities
    for sectionName in selectedSections:

        if sectionName not in list(sectionDict.keys()):

            # User wants a section that is not defined. Print a warning message:
            print("ERROR: Surface", sectionName, "is not defined. Check surface names in your CGNS file.")
            quit()

        elif "triaConnF" in list(sectionDict[sectionName].keys()):  # Then we have a surface section

            # Assign connectivities
            if triaConnF is None:
                # Start new connectivities if we have none
                triaConnF = sectionDict[sectionName]["triaConnF"]
                quadsConnF = sectionDict[sectionName]["quadsConnF"]
            else:
                # Append new connectivities
                triaConnF = np.vstack([triaConnF, sectionDict[sectionName]["triaConnF"]])
                quadsConnF = np.vstack([quadsConnF, sectionDict[sectionName]["quadsConnF"]])

        else:

            # The user provided a name that is not a surface section
            print(sectionName, "is not a surface section.")

    return triaConnF, quadsConnF


# =================================================================

# =================================================================
# AUXILIARY CURVE FUNCTIONS
# =================================================================

# =================================================================


def create_curve_from_points(coor, curveName, periodic=False, mergeTol=1e-7):

    """
    This function will generate a curve object from a list
    of points provided by the user.
    This function automatically generate the bar element connectivity
    assuming that the points in coor are ordered.

    Parameters
    ----------

    coor -> float[numNodes,3] : coordinates (X,Y,Z) of the curve nodes.

    curveName -> string : name for the new curve object.

    periodic -> Boolean : if True, a bar element will connect the first and
                          the last node defined in coor.

    mergeTol -> float : threshold used to merge nearby nodes.

    Returns
    -------

    curve -> TSurfCurve object : curve object defined by the given nodes.

    Ney Secco 2017-02
    """

    # Check if we have 3D data
    if coor.shape[1] != 3:
        raise ValueError("coor should be of shape [numNodes,3]")

    # Determine the number of nodes
    numNodes = coor.shape[0]

    # Create bar element connectivity
    barsConn = np.zeros((numNodes - 1, 2))
    barsConn[:, 0] = list(range(numNodes - 1))
    barsConn[:, 1] = barsConn[:, 0] + 1

    # Add an extra element if the curve is periodic
    if periodic:
        newBar = np.array([[barsConn[-1, -1], 0]])
        np.vstack([barsConn, newBar])

    # Create curve object
    curve = tsurf_component.TSurfCurve(coor, barsConn, curveName, mergeTol)

    # Return curve object
    return curve


# =================================================================


def merge_curves(curveDict, mergedCurveName, curvesToMerge=None):

    """
    This function will merge curves in the curveDict dictonary whose
    name is included in curveNames.
    If curvesToMerge=None, we will merge all curves in curveDict

    Ney Secco 2016-11
    """

    # Initialize set of coordinates and connectivities
    mergedCoor = np.zeros((0, 3))
    mergedBarsConn = np.zeros((0, 2), dtype="int32")

    # Merge all curves if user provided none
    if curvesToMerge == None:
        curvesToMerge = list(curveDict.keys())

    # Initialize index offset. The new curves cannot reference
    # nodes of the old ones.
    indexOffset = 0

    # Loop over every curve we want to merge
    for curveName in curvesToMerge:

        # Check if we have access to this curve
        if curveName in curveDict:

            # Get coor and connectivities of this curve.
            # Remember to apply the offset to the connectivities
            coor = curveDict[curveName].coor
            barsConn = curveDict[curveName].barsConn + indexOffset

            # Update offset for next curve
            indexOffset = indexOffset + coor.shape[0]

            # Merge data
            mergedCoor = np.vstack([mergedCoor, coor])
            mergedBarsConn = np.vstack([mergedBarsConn, barsConn])

    # Create new curve with the merged data
    mergedCurve = tsurf_component.TSurfCurve(mergedCoor, mergedBarsConn, mergedCurveName)

    # Return merged curve
    return mergedCurve


# =================================================================


def split_curves(curveDict, curvesToBeSplit=None, optionsDict={}, criteria="sharpness"):

    """
    This function will loop over all curves in curveDict and split
    them based on a given criteria. The available criteria are:
    criteria=['sharpness', 'curve']

    Returns
    =======
    breakList, list -> List of bar element IDs where the guideCurves split
                       the initial intersection. These are outputted and can
                       be used to always split the curves at the same bar
                       elements, which avoids an issue where the mesh may
                       become discontinuous.

    Ney Secco 2016-08
    John Jasa 2016-11
    """

    # Split all curves in the dictionary if the user provided None as input
    if curvesToBeSplit is None:
        curvesToBeSplit = list(curveDict.keys())

    """
    # Set the default curvesToBeSplit; useful if criteria = 'sharpness'
    options = {'curvesToBeSplit' : []}
    options.update(optionsDict)
    """

    # Loop over every curve to check for splits
    for curveName in list(curveDict.keys()):

        if curveName in curvesToBeSplit:

            # First we remove the curve component from the dictionary
            curve = curveDict.pop(curveName)

            # Now we run the split function for this single curve
            splitcurves = split_curve_single(curve, curveName, optionsDict=optionsDict, criteria=criteria)

            # Now we add the split curves to the original curves dictionary
            curveDict.update(splitcurves)


# =============================================================


def split_curve_single(curve, curveName, optionsDict={}, criteria="sharpness"):

    """
    This function receives a single curve object, splits it according to a criteria,
    and then return a new dictionary with split curves.

    ATTENTION: This function assumes that the FE data is sorted.

    Parameters
    ----------
    curve: Curve object -> Curve that will be split

    curveName: string -> Name of the original curve. This name will be used
               to name the split curves.

    optionsDict: dictionary -> Dictionary containing special options for each split
    criteria. The following fields are needed:
    For 'sharpness' criteria: -'angle': angle (in degrees) used to split the curve.
                               the default value is 60 deg.
    For 'curve' criteria: -'splittingCurves': list of curve objects. We will split
                           the current curve in the intersection points with other curves.
    For 'node' criteria: -'splittingNodes': list of integers. They represent the node IDs
                          where we will split the curve.

    criteria: string -> Criteria that will be used to split curves. The options
              available for now are: ['sharpness', 'curve', 'node']

    Returns
    -------
    splitCurveDict: dictionary[curve objects] -> Dictionary containing split curves.

    Ney Secco 2016-08
    """

    # PARAMETERS
    sharpAngle = 60 * np.pi / 180.0

    # READING INPUTS

    # Get coordinates and connectivities
    coor = curve.coor
    barsConn = curve.barsConn

    # Get the number of elements
    nElem = barsConn.shape[0]

    isPeriodic = False

    # DETECT KINKS

    # Initialize list of break elements. The break is considered to be
    # at the first point of the element
    breakList = []

    # The next step will depend on the criteria

    if criteria == "sharpness":

        # This criteria will split the curve if its sharpness (defined here
        # as the change in tangential direction) is beyond a given threshold.

        # Get the split angle from the options
        if "angle" in list(optionsDict.keys()):
            sharpAngle = optionsDict["angle"] * np.pi / 180.0
        else:
            sharpAngle = 60 * np.pi / 180.0

        # Get the tangent direction of the first bar element
        prevTan = coor[barsConn[0, 1], :] - coor[barsConn[0, 0], :]
        prevTan = prevTan / np.linalg.norm(prevTan)

        # Loop over the remaining bars to find sharp kinks
        for elemID in range(1, nElem):

            # Compute tangent direction of the current element
            currTan = coor[barsConn[elemID, 1], :] - coor[barsConn[elemID, 0], :]
            currTan = currTan / np.linalg.norm(currTan)

            # Compute change in direction between consecutive tangents
            angle = np.arccos(np.min([np.max([prevTan.dot(currTan), -1.0]), 1.0]))

            # Check if the angle is beyond a certain threshold
            if angle > sharpAngle:

                # Store the current element as a break position
                breakList.append(elemID)

            # Now the current tangent will become the previous tangent of
            # the next iteration
            prevTan = currTan.copy()

    # If the curve is periodic, we need to check for breaks between the first and last elements.
    if barsConn[0, 0] == barsConn[-1, 1]:

        # Compute tangent direction of the current element
        prevTan = coor[barsConn[-1, 1], :] - coor[barsConn[-1, 0], :]
        prevTan = prevTan / np.linalg.norm(prevTan)

        # Get the tangent direction of the first bar element
        currTan = coor[barsConn[0, 1], :] - coor[barsConn[0, 0], :]
        currTan = currTan / np.linalg.norm(currTan)

        # Compute change in direction between consecutive tangents
        angle = np.arccos(np.min([np.max([prevTan.dot(currTan), -1.0]), 1.0]))

        # Check if the angle is beyond a certain threshold
        if angle > sharpAngle:

            # We should break the curve at the initial point
            # so it will not be periodic anymore
            isPeriodic = False

        else:
            isPeriodic = True

    if criteria == "curve":
        if "splittingCurves" in list(optionsDict.keys()):
            for guideCurve in optionsDict["splittingCurves"]:
                split_index = closest_node(guideCurve, coor)

                elemID = np.where(barsConn[:, 0] == split_index)[0]

                # Store the current element as a break position
                breakList.append(elemID[0])

    elif criteria == "node":
        if "splittingNodes" in list(optionsDict.keys()):
            for node in optionsDict["splittingNodes"]:
                breakList.append(node)
        else:
            raise ValueError("node IDs are not provided. Cannot split curve.")

    # Make sure that the obtained breakList is in order so that we get
    # the correct split curves
    breakList = sorted(breakList)

    if 0 in breakList:
        breakList.remove(0)

    # CHECK IF BREAKS WERE DETECTED
    # We can stop this function earlier if we detect no break points
    if len(breakList) == 0:

        # In this case, we just create a dictionary with the original curve

        # Initialize dictionary that will contain the split curves
        splitcurvesDict = {}

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # Slice the original connectivity matrix
        splitBarsConn = barsConn[:, :]

        # Generate a name for this new curve
        splitCurveName = curveName + "_split"

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitcurvesDict[splitCurveName] = splitCurve

        # Return dictionary with the single curve
        return splitcurvesDict

    # CREATE SPLIT CURVE OBJECTS

    # Now that we have the list of break points (given in breakList), we can create multiple
    # curve objects.

    # Count the number of curves that should be between the first and last break point
    nInnercurves = len(breakList) - 1

    # Initialize dictionary that will contain the split curves
    splitcurvesDict = {}

    # Initialize dictionary that will hold split information
    curve.extra_data["splitCurves"] = {}

    # First we will define all curves that are between the first and the last
    # break point. The remaining part of the curve will be determined depending
    # if the original curve is periodic or not.
    for splitID in range(nInnercurves):

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # Slice the original connectivity matrix
        splitBarsConn = barsConn[breakList[splitID] : breakList[splitID + 1], :]

        # Generate a name for this new curve
        splitCurveName = curveName + "_" + "%02d" % (splitID + 1)

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitcurvesDict[splitCurveName] = splitCurve

        # Store split information in the parent curve
        curve.extra_data["splitCurves"][splitCurveName] = [breakList[splitID], breakList[splitID + 1]]

    # Now we need to create curves with elements that come before the first break point and after
    # the last break point.
    # We need to treat periodic curves differently. We use the same check to determine
    # the number of split curves. The periodic case will have minus one curve compared
    # to a non-periodic one.
    if isPeriodic:  # Curve is periodic, so we need to define one curve

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = np.vstack([barsConn[breakList[-1] :, :], barsConn[: breakList[0], :]])

        # Generate a name for this new curve
        splitCurveName = curveName + "_" + "%02d" % 0

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitcurvesDict[splitCurveName] = splitCurve

        # Store split information in the parent curve
        curve.extra_data["splitCurves"][splitCurveName] = [breakList[-1], breakList[0]]

    else:  # Curve is not periodic, so we need to define two curves

        # CURVE 0 : before the first break point

        # For now copy the original set of nodes (the unused ones will be removed)
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = barsConn[: breakList[0], :]

        # Generate a name for this new curve
        splitCurveName = curveName + "_" + "%02d" % 0

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitcurvesDict[splitCurveName] = splitCurve

        # Store split information in the parent curve
        curve.extra_data["splitCurves"][splitCurveName] = [0, breakList[0]]

        # CURVE 1 : after the first break point

        # For now copy the original set of nodes
        splitCoor = coor.copy()

        # We need to wrap around connectivities
        splitBarsConn = barsConn[breakList[-1] :, :]

        # Generate a name for this new curve
        splitCurveName = curveName + "_" + "%02d" % (nInnercurves + 1)

        # Create curve object
        splitCurve = tsurf_component.TSurfCurve(splitCoor, splitBarsConn, splitCurveName)

        # Append new curve object to the dictionary
        splitcurvesDict[splitCurveName] = splitCurve

        # Store split information in the parent curve
        curve.extra_data["splitCurves"][splitCurveName] = [breakList[-1], barsConn.shape[0]]

    # Return the dictionary of new curves and the list of breaking elements
    return splitcurvesDict


# =================================================================

# ===================================
# INTERSECTION FUNCTIONS
# ===================================


def compute_intersections(TSurfGeometryList, distTol=1e-7, comm=MPI.COMM_WORLD):

    """
    This function will just compute pair-wise intersections
    for all components given in TSurfGeometryList.

    distTol is a distance tolerance used to merge nearby nodes when
    generating the intersection finite element data.
    """

    # Get number of components
    numGeometries = len(TSurfGeometryList)

    # Stop if user gives only one component
    if numGeometries < 2:
        print("ERROR: Cannot compute intersections with just one component")
        quit()

    # Initialize dictionary of intersections
    Intersections = {}

    # Call intersection function for each pair
    for ii in range(numGeometries):
        for jj in range(ii + 1, numGeometries):

            # Compute new intersections for the current pair
            newIntersections = _compute_pair_intersection(TSurfGeometryList[ii], TSurfGeometryList[jj], distTol, comm)

            # Append new curve objects to the dictionary
            for curve in newIntersections:
                Intersections[curve.name] = curve

    # Print log
    print("Computed", len(Intersections), "intersection curves.")

    # Return all intersections
    return Intersections


# =================================================================


def _compute_pair_intersection(TSurfGeometryA, TSurfGeometryB, distTol):

    """
    This function finds intersection curves between components A and B
    and returns a list of curve objects corresponding to these intersections.
    """

    # Get communicator from the first object
    comm = TSurfGeometryA.comm

    ### Call Fortran code to find intersections

    # This will de done in two steps:
    # The first Fortran call will give us a tuple containing the number of nodes, bars, and parent triangles
    # These are the sizes of the allocatable arrays that were defined in the Fortran call.
    # The second Fortran call will use this information to actually retrieve the intersection data.
    # We have to use this two steps approach so that Python does not need direct access to the
    # allocatable arrays, since this does not work when we use Intel compilers.

    arraySizes = intersectionAPI.intersectionapi.computeintersection(
        TSurfGeometryA.coor.T,
        TSurfGeometryA.triaConnF.T,
        TSurfGeometryA.quadsConnF.T,
        TSurfGeometryB.coor.T,
        TSurfGeometryB.triaConnF.T,
        TSurfGeometryB.quadsConnF.T,
        distTol,
        comm.py2f(),
    )

    # Retrieve results from Fortran if we have an intersection
    if np.max(arraySizes[1:]) > 0:

        # Second Fortran call to retrieve data from the CGNS file.
        intersectionArrays = intersectionAPI.intersectionapi.retrievedata(*arraySizes)

        # Split results.
        # We need to do actual copies, otherwise data will be overwritten if we compute another intersection.
        # We subtract one to make indices consistent with the Python 0-based indices.
        # We also need to transpose it since Python and Fortran use different orderings to store matrices
        # in memory.
        coor = np.array(intersectionArrays[0]).T
        barsConn = np.array(intersectionArrays[1]).T - 1
        parentTria = np.array(intersectionArrays[2]).T - 1

    else:

        # Assign empty arrays
        coor = np.zeros((0, 3))
        barsConn = np.zeros((0, 2))
        parentTria = np.zeros((0, 2))

    # Release memory used by Fortran so we can run another intersection in the future
    intersectionAPI.intersectionapi.releasememory()

    # Initialize list to hold all intersection curves
    Intersections = []

    if len(coor) > 0:

        # Sort FE data. After this step, barsConn may become a list if we
        # have two disconnect intersection curves.
        barsConn, newMap = FEsort(barsConn.tolist())

        # Now barsConn is a list of curves. Each element of the list brings
        # FE data for a continuous curve. If the intersection generates
        # multiple disconnect curves, barsConn will be a list with
        # multiple elements.
        # For each curve (and therefore for each barsConn element), we
        # will initialize a Curve object to hold intersection FE data.

        # Set intersection counter
        intCounter = -1

        for (currConn, currMap) in zip(barsConn, newMap):

            # Increment counter
            intCounter = intCounter + 1

            # Gather names of parent components
            name1 = TSurfGeometryA.name
            name2 = TSurfGeometryB.name

            # Create a name for the curve
            curveName = "int_" + name1 + "_" + name2 + "_%03d" % intCounter

            # Slice the parent triangles array using the sorted mapping
            currParents = parentTria[currMap, :]

            # Create new curve object
            newCurve = tsurf_component.TSurfCurve(coor, currConn, curveName)

            # Store parent triangles as extra data
            newCurve.extra_data["parentTria"] = currParents
            newCurve.extra_data["parentGeoms"] = [name1, name2]  # This will also be done by the manager

            # Initialize curve object and append it to the list
            Intersections.append(newCurve)

    # Return intersection FE data
    return Intersections


# =================================================================


def _compute_pair_intersection_d(TSurfGeometryA, TSurfGeometryB, intCurve, coorAd, coorBd, distTol):

    """
    This function is the backward mode of the intersection computation.
    This can only be called when we already have the intersection curve.

    Parameters
    ----------

    TSurfGeometryA: TSurfGeometry object

    TSurfGeometryB: TSurfGeometry object

    intCurve: Intersection curve defined by both components (taken from component A)

    coorAd: float(nNodesA,3) -> Derivative seeds of the component A nodal coordinates

    coorBd: float(nNodesB,3) -> Derivative seeds of the component B nodal coordinates

    distTol: float -> Distance used to check if the bar elements were flipped. Used the
                      same distTol used to merge nearby nodes.

    Returns
    -------

    coorIntd: float(nNodesInt,3) -> Derivative seeds of the intersection curve nodal coordinates

    Ney Secco 2016-09
    """

    # Call Fortran code to find derivatives
    # The transposes are because Fortran stores data by rows while Python stores by columns.
    # The +1 are due to the Fortran indexing, that starts at 1
    coorIntd = intersectionAPI.intersectionapi.computeintersection_d(
        TSurfGeometryA.coor.T,
        coorAd.T,
        TSurfGeometryA.triaConnF.T,
        TSurfGeometryA.quadsConnF.T,
        TSurfGeometryB.coor.T,
        coorBd.T,
        TSurfGeometryB.triaConnF.T,
        TSurfGeometryB.quadsConnF.T,
        intCurve.coor.T,
        intCurve.barsConn.T + 1,
        intCurve.extra_data["parentTria"].T + 1,
        distTol,
    )

    # Return derivatives
    return coorIntd.T


# =================================================================


def _compute_pair_intersection_b(TSurfGeometryA, TSurfGeometryB, intCurve, coorIntb, distTol):

    """
    This function is the backward mode of the intersection computation.
    This can only be called when we already have the intersection curve.

    Parameters
    ----------

    TSurfGeometryA: TSurfGeometry object

    TSurfGeometryB: TSurfGeometry object

    intCurve: Intersection curve defined by both components (taken from component A)

    coorIntb: float(nNodesInt,3) -> Derivative seeds of the intersection curve nodal coordinates

    distTol: float -> Distance used to check if the bar elements were flipped. Used the
                      same distTol used to merge nearby nodes.

    Returns
    -------

    coorAb: float(nNodesA,3) -> Derivative seeds of the component A nodal coordinates

    coorBb: float(nNodesB,3) -> Derivative seeds of the component B nodal coordinates

    Ney Secco 2016-09
    """

    # Call Fortran code to find derivatives
    # The transposes are because Fortran stores data by rows while Python stores by columns.
    # The +1 are due to the Fortran indexing, that starts at 1
    coorAb, coorBb = intersectionAPI.intersectionapi.computeintersection_b(
        TSurfGeometryA.coor.T,
        TSurfGeometryA.triaConnF.T,
        TSurfGeometryA.quadsConnF.T,
        TSurfGeometryB.coor.T,
        TSurfGeometryB.triaConnF.T,
        TSurfGeometryB.quadsConnF.T,
        intCurve.coor.T,
        coorIntb.T,
        intCurve.barsConn.T + 1,
        intCurve.extra_data["parentTria"].T + 1,
        distTol,
    )

    # Return derivatives
    return coorAb.T, coorBb.T


# =================================================================

# ===================================
# GENERAL AUXILIARY FUNCTIONS
# ===================================


def formatStringArray(fortranArray):

    """
    This functions receives a numpy array of characters coming from Fortran
    and returns a list of formatted data.
    """

    # Get array size
    shape = fortranArray.shape

    # We need to tranpose the array considering the Fortran ordering
    fortranArray = fortranArray.T

    # Initialize new array
    pythonArray = []

    # Now we format and concatenate the elements
    for index in range(shape[0]):

        # format string
        newString = "".join(str(fortranArray[:, index], "utf-8")).strip().lower()

        # Add formatted string to the list
        pythonArray.append(newString)

    # Return formatted array
    return pythonArray


# =============================================================


def FEsort(barsConn):

    """
    This function can be used to sort connectivities coming from the CGNS file

    barsConn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

    newConn will be set to None if the sorting algorithm fails.
    """

    # We will solve this in an iterative process until the connectivities do not change
    # anymore

    # Initialize the newest connectivities as our current input.
    newConn = barsConn

    # Say that we still need to search
    keep_searching = True

    # Get the number of bars
    nBars = len(barsConn)

    # Initialize mapping that will link indices of the given FE data to the sorted one.
    # This mapping will tell which curve the bar element was assigned to.
    # For instance, newMap[5] gives the barIDs that belong to curve[5].
    newMap = [[i] for i in range(nBars)]

    while keep_searching:

        # If nothing happens, we will get out of the loop
        keep_searching = False

        # Update "old" connectivities and mapping
        oldConn = newConn[:]
        oldMap = newMap[:]

        # Initialize list of new connectivities and mapping with the first element of the old ones.
        # We will also pop this element from oldConn. Remember that in subsequent iterations, each
        # element of oldConn represents in fact a concatenation of bar elements.
        newConn = [oldConn.pop(0)]
        newMap = [oldMap.pop(0)]

        # We will keep sorting until all elements of oldConn are assigned and popped out
        # NOTE: An "Element" here may be a concatenation of several bar elements

        while len(oldConn) > 0:

            # Pop another element from oldConn
            oldElement = oldConn.pop(0)
            oldElemMap = oldMap.pop(0)

            # We do not know if this element could be linked beforehand
            linked_element = False

            # Now we check if we can fit this element into any new connectivity
            newElemCounter = 0
            for (newElement, newElemMap) in zip(newConn, newMap):

                if oldElement[0] == newElement[0]:

                    # We will flip the old element and place it at the beginning of the
                    # new connectivity
                    newConn[newElemCounter] = oldElement[::-1] + newElement[1:]

                    # Also assign ID the beginning of the map
                    newMap[newElemCounter] = oldElemMap[::-1] + newElemMap[:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[0] == newElement[-1]:

                    # Just append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement

                    # Also assign mapping
                    newMap[newElemCounter] = newElemMap[:] + oldElemMap

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[0]:

                    # Place the old element at the beginning
                    newConn[newElemCounter] = oldElement + newElement[1:]

                    # Also assign mapping
                    newMap[newElemCounter] = oldElemMap + newElemMap[:]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                elif oldElement[-1] == newElement[-1]:

                    # Flip and append the old element
                    newConn[newElemCounter] = newElement[:-1] + oldElement[::-1]

                    # Also assign mapping
                    newMap[newElemCounter] = newElemMap[:] + oldElemMap[::-1]

                    # We need to keep searching as we are still making changes
                    linked_element = True

                    # We only need to keep searching if we have updates
                    keep_searching = True

                    # Jump out of the for loop
                    break

                # Increment element counter
                newElemCounter = newElemCounter + 1

            if not linked_element:  # We have an element that is not connected to anyone

                # Define a new curve in newConn
                newConn.append(oldElement)

                # Define a new mapping in newMap
                newMap.append(oldElemMap)

    # Right now, newConn represent each line by an 1D array of point IDs.
    # e.g. [[2,3,4,5],[1,7,8]]
    # We need to convert this back to FE format, represented by 2D arrays.
    # e.g. [[[2,3,4],[3,4,5]],[[1,7],[7,8]]]

    # Initialize FE array for each curve
    newConnFE = [np.zeros((len(curve) - 1, 2), dtype=int) for curve in newConn]

    # Fill FE data for each curve
    for curveID in range(len(newConn)):

        # Get current curve
        curve = newConn[curveID]

        # Get FE list of the current curve
        FEcurve = newConnFE[curveID]

        # Assign FE connectivity
        for pointID in range(1, len(curve)):

            # Get point indices
            prevPoint = curve[pointID - 1]
            currPoint = curve[pointID]

            # Assign bar FE
            FEcurve[pointID - 1, 0] = prevPoint
            FEcurve[pointID - 1, 1] = currPoint

    # We just need to crop the last index of the mapping arrays because it has a repeated
    # number

    # Now we do a final check to remove degenerate bars (bars that begin and end at the same point)
    # Take every disconnect curve in newConnFE:
    for curveID in range(len(newConnFE)):

        # We convert the connectivity array to list so we can 'pop' elements
        curveFE = newConnFE[curveID].tolist()

        # Get the mapping as well
        currMap = newMap[curveID]

        # Define FE counter
        FEcounter = 0

        while FEcounter < len(curveFE):

            # Get current finite element
            FE = curveFE[FEcounter]

            # Check if the start and end points are the same
            if FE[0] == FE[1]:

                # Remove FE
                curveFE.pop(FEcounter)

                # Remove mapping (remove bar element from this curve)
                currMap.pop(FEcounter)

            else:

                # Increment counter
                FEcounter = FEcounter + 1

        # Convert connectivity back to numpy array
        newConnFE[curveID] = np.array(curveFE)

    # Now remove empty curves
    curveID = 0
    while curveID < len(newConnFE):

        # Check if current curve has no element
        if len(newConnFE[curveID].tolist()) == 0:

            # If this is the case, remove the curve and its mapping
            newConnFE.pop(curveID)
            newMap.pop(curveID)

        else:

            # Increment counter
            curveID = curveID + 1

    # Return the sorted array and the mapping
    return newConnFE, newMap


# =============================================================


def remove_unused_points(coor, triaConnF=np.zeros((0, 0)), quadsConnF=np.zeros((0, 0)), barsConn=np.zeros((0, 0))):

    """
    This function will removed unused points from coor and
    also update all connectivities.

    ATTENTION:
    All inputs must be arrays. The user does not need to provide all inputs, only
    the connectivities that matter to select used points.

    Returns
    -------
    This function returns cropCoor, which is the cropped set of nodes.
    This function also updates triaConn, quadsConn and barsConn.
    """

    # Get total number of points and elements
    nPoints = coor.shape[0]
    nTria = triaConnF.shape[0]
    nQuads = quadsConnF.shape[0]
    nBars = barsConn.shape[0]

    # Initialize mask to identify used points
    usedPtsMask = np.zeros(nPoints, dtype=int)

    # First we loop over all elements to create a mask that indicates used points
    for triaID in range(nTria):
        # Flag used points
        # Remember that triaConnF is in Fortran indexing (starting at 1)
        # So we need to subtract 1 to get correct values
        usedPtsMask[triaConnF[triaID, 0] - 1] = 1
        usedPtsMask[triaConnF[triaID, 1] - 1] = 1
        usedPtsMask[triaConnF[triaID, 2] - 1] = 1

    for quadID in range(nQuads):
        # Flag used points
        # Remember that quadsConnF is in Fortran indexing (starting at 1)
        # So we need to subtract 1 to get correct values
        usedPtsMask[quadsConnF[quadID, 0] - 1] = 1
        usedPtsMask[quadsConnF[quadID, 1] - 1] = 1
        usedPtsMask[quadsConnF[quadID, 2] - 1] = 1
        usedPtsMask[quadsConnF[quadID, 3] - 1] = 1

    for barID in range(nBars):
        # Flag used points
        usedPtsMask[barsConn[barID, 0]] = 1
        usedPtsMask[barsConn[barID, 1]] = 1

    # Now we can compute the number of points actually used
    nUsedPts = np.sum(usedPtsMask)

    # Initialize new coordinate array
    cropCoor = np.zeros((nUsedPts, 3))

    # Initialize counter to fill cropCoor
    cropPointID = -1

    # Now we fill the points of the cropped array
    for pointID in range(nPoints):

        # Check if the point is used
        if usedPtsMask[pointID] == 1:

            # Increment counter
            cropPointID = cropPointID + 1

            # Add point to the cropped array
            cropCoor[cropPointID, :] = coor[pointID, :]

            # Now we replace the value in the mask array so that we
            # can use it as a pointer from coor to cropCoor when we
            # update element connectivities.
            usedPtsMask[pointID] = cropPointID

        else:

            # Unused points receive a negative flag
            usedPtsMask[pointID] = -1

    # Now we need to update connectivities so that they point to the the correct
    # indices of the cropped array
    for triaID in range(nTria):
        # Use pointer to update connectivity
        # Remember that triaConnF is in Fortran indexing (starting at 1)
        # So we need to subtract 1 to get correct values.
        triaConnF[triaID, 0] = usedPtsMask[triaConnF[triaID, 0] - 1] + 1
        triaConnF[triaID, 1] = usedPtsMask[triaConnF[triaID, 1] - 1] + 1
        triaConnF[triaID, 2] = usedPtsMask[triaConnF[triaID, 2] - 1] + 1

    for quadID in range(nQuads):
        # Use pointer to update connectivity
        # Remember that quadsConnF is in Fortran indexing (starting at 1)
        # So we need to subtract 1 to get correct values
        quadsConnF[quadID, 0] = usedPtsMask[quadsConnF[quadID, 0] - 1] + 1
        quadsConnF[quadID, 1] = usedPtsMask[quadsConnF[quadID, 1] - 1] + 1
        quadsConnF[quadID, 2] = usedPtsMask[quadsConnF[quadID, 2] - 1] + 1
        quadsConnF[quadID, 3] = usedPtsMask[quadsConnF[quadID, 3] - 1] + 1

    for barID in range(nBars):
        # Use pointer to update connectivity
        barsConn[barID, 0] = usedPtsMask[barsConn[barID, 0]]
        barsConn[barID, 1] = usedPtsMask[barsConn[barID, 1]]

    # Return the new set of nodes and the mask used to map nodes
    return cropCoor, usedPtsMask


# =============================================================


def detect_feature(node1, node2, element1, element2, coor, triaConnF, quadsConnF, feature):

    """
    This function checks if the bar that connects node1 and node2, and
    is shared between element1 and element2 has a desired feature. This is
    used by extract_curves_from_surface function.

    Parameters
    ----------

    node1: integer -> node1 index in coor

    node2: integer -> node2 index in coor

    element1: integer -> element index in triaConn (if positive) or quadsConn
              (if negative).

    element2: integer -> element index in triaConn (if positive) or quadsConn
              (if negative). Could be None if bar is not linked to any other element.

    coor: float[nNodes,3] -> nodal coordinates

    triaConnF: integer[nTria,3] -> triangle connectivities. Remember that indices should
                                   start at 1, as in Fortran ordering

    quadsConnF: integer[nQuads,4] -> quad connectivities. Remember that indices should
                                   start at 1, as in Fortran ordering

    feature: string -> feature that should be detected. Available options are:
             ['sharpness','open_ends']

    Returns
    -------

    featureIsDetected: logical -> True if feature is detected in this edge. Otherwise it is False.

    Ney Secco 2016-08
    """

    # Check the criteria

    if feature == "sharpness":

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
        if element2 == None:
            featureIsDetected = False
            return featureIsDetected

        elif element2 > 0:  # We have a tria

            # Get element nodes
            node1Coor = coor[triaConnF[element2, 0] - 1, :]
            node2Coor = coor[triaConnF[element2, 1] - 1, :]
            node3Coor = coor[triaConnF[element2, 2] - 1, :]

            # Get inplane vectors so that the normal points outside
            v12 = node2Coor - node1Coor
            v22 = node3Coor - node1Coor

        elif element2 < 0:  # We have a quad

            # Get element nodes
            node1Coor = coor[quadsConnF[-element2, 0] - 1, :]
            node2Coor = coor[quadsConnF[-element2, 1] - 1, :]
            node3Coor = coor[quadsConnF[-element2, 2] - 1, :]
            node4Coor = coor[quadsConnF[-element2, 3] - 1, :]

            # Get inplane vectors so that the normal points outside
            v12 = node3Coor - node1Coor
            v22 = node4Coor - node2Coor

        # ELEMENT 1
        if element1 > 0:  # We have a tria

            # Get element nodes
            node1Coor = coor[triaConnF[element1, 0] - 1, :]
            node2Coor = coor[triaConnF[element1, 1] - 1, :]
            node3Coor = coor[triaConnF[element1, 2] - 1, :]

            # Get inplane vectors so that the normal points outside
            v11 = node2Coor - node1Coor
            v21 = node3Coor - node1Coor

        elif element1 < 0:  # We have a quad

            # Get element nodes
            node1Coor = coor[quadsConnF[-element1, 0] - 1, :]
            node2Coor = coor[quadsConnF[-element1, 1] - 1, :]
            node3Coor = coor[quadsConnF[-element1, 2] - 1, :]
            node4Coor = coor[quadsConnF[-element1, 3] - 1, :]

            # Get inplane vectors so that the normal points outside
            v11 = node3Coor - node1Coor
            v21 = node4Coor - node2Coor

        # Take the cross products to get normals and normalize the normals =P
        # Element 1
        n1 = np.cross(v11, v21)
        n1 = n1 / np.linalg.norm(n1)
        # Element 2
        n2 = np.cross(v12, v22)
        n2 = n2 / np.linalg.norm(n2)

        # Use dot product to find the angle between the normals
        angle = np.arccos(n1.dot(n2))

        # We have a "sharp" edge if this angle is beyond a certain threshold
        if angle > 60 * np.pi / 180:
            featureIsDetected = True
            return featureIsDetected
        else:
            featureIsDetected = False
            return featureIsDetected

    elif feature == "open_ends":

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

        print("ERROR: Feature", feature, "cannot be detected as it is not an option.")
        quit()


# ===================================
# CURVE MESHING FUNCTIONS
# ===================================

# Define hyperbolic tangent generator
def hypTanDist(Sp1, Sp2, N):

    # This is the hyperbolic tangential spacing developed by ICEMCFD
    # I got the equations from the ICEM manual
    # This bunching law is coarser at the middle and finer at the ends
    # of the interval, just like shown below:
    # || |   |    |    |   |  | ||

    # Sp1: initial spacing (within the [0,1] interval)
    # Sp2: final spacing (within the [0,1] interval)
    # N: number of nodes
    #
    # Ney Secco 2016-11

    # IMPORTS
    from numpy import arange, sqrt, sinh, tanh
    from scipy.optimize import broyden1

    # Check existence of result
    if (N - 1) * sqrt(Sp1 * Sp2) < 1:
        print("Hyperbolic distribution not possible")
        # exit()

    # First find b
    def findRootb(b):
        return sinh(b) - b / (N - 1) / sqrt(Sp1 * Sp2)

    b = broyden1(findRootb, 20, f_tol=1e-10)

    # # Manually created secant method for the above solve
    # b = 4.
    # step = 1.e-6
    # for i in range(1000):
    #     b_old = b
    #     f_prime = (findRootb(b) - findRootb(b-step)) / step
    #     b = b - findRootb(b) / f_prime
    #     if np.abs(b_old - b) < 1.e-10:
    #         break

    # Compute parameter A
    A = sqrt(Sp1 / Sp2)

    # Create nodes indices
    index = arange(N) + 1

    # Compute parameter R
    R = (index - 1) / (N - 1) - 1 / 2

    # Compute parameter U
    U = 1 + tanh(b * R) / tanh(b / 2)

    # Compute spacings
    S = U / (2 * A + (1 - A) * U)

    # Return spacings
    return S


####################################################
####################################################
####################################################

# Define tangent generator
def tanDist(Sp1, Sp2, N):

    # This is the tangential spacing developed by Ney Secco
    # This bunching law is coarser at the ends and finer at the middle
    # of the interval, just like shown below:
    # |    |   |  | || |  |   |    |

    # Sp1: initial spacing (within the [0,1] interval)
    # Sp2: final spacing (within the [0,1] interval)
    # N: number of nodes
    #
    # Ney Secco 2016-11

    # IMPORTS
    from numpy import tan, arange, pi

    try:
        from scipy.optimize import minimize
    except:
        print("ERROR: Scipy not available")

    # Convert number of nodes to number of cells
    N = N - 1

    # Define objective function
    def func(P):
        # Split variables
        a = P[0]
        e = P[1]
        c = P[2]
        # Find b
        b = e - c
        # Equations
        Eq1 = a * (tan(b + c) - tan(c)) - 1
        Eq2 = a * (tan(b / N + c) - tan(c)) - Sp1
        Eq3 = a * (tan(b + c) - tan(b * (1 - 1 / N) + c)) - Sp2
        # Cost function
        J = Eq1 ** 2 + Eq2 ** 2 + Eq3 ** 2
        # Return
        return J

    # Define bounds for the problem
    a_bounds = [(0, None)]
    e_bounds = [(0, pi / 2)]
    c_bounds = [(-pi / 2, 0)]
    bounds = a_bounds + e_bounds + c_bounds

    # Define initial guess
    a_start = 1.0
    e_start = pi / 4
    c_start = -pi / 4
    x_start = [a_start, e_start, c_start]

    # Optimize
    res = minimize(
        func, x_start, method="SLSQP", bounds=bounds, options={"maxiter": 1000, "disp": False, "ftol": 1e-12}
    )

    # Split variables
    a = res.x[0]
    e = res.x[1]
    c = res.x[2]

    # Find other parameters
    b = e - c
    d = -a * tan(c)

    # Generate spacing
    index = arange(N + 1)
    S = a * tan(b * index / N + c) + d

    # Return spacing
    return S


####################################################
####################################################
####################################################

# Define cubic generator
def cubicDist(Sp1, Sp2, N):

    # Sp1: initial spacing (within the [0,1] interval)
    # Sp2: final spacing (within the [0,1] interval)
    # N: number of nodes
    #
    # Ney Secco 2016-11

    # IMPORTS
    from numpy import array, arange
    from numpy.linalg import inv

    # Convert number of nodes to number of cells
    N = N - 1

    # Assemble linear system matrix
    A = array([[1, 1, 1], [(1 / N) ** 3, (1 / N) ** 2, (1 / N)], [(1 - 1 / N) ** 3, (1 - 1 / N) ** 2, (1 - 1 / N)]])
    b = array([1, Sp1, 1 - Sp2])

    # Solve the linear system
    x = inv(A).dot(b)
    a = x[0]
    b = x[1]
    c = x[2]

    # Generate spacing
    index = arange(N + 1) / N
    S = a * index ** 3 + b * index ** 2 + c * index

    # Return spacing
    return S


# ============================================================

# ============================================================
# ============================================================
# ============================================================


def translate(Geometry, x, y, z):
    """
    Translate coor based on given x, y, and z values.
    John Jasa 2016-08
    """

    Geometry.coor = Geometry.coor[:, :] + np.atleast_2d(np.array([x, y, z]))


def scale(Geometry, factor, point=None):
    """
    Scale coor about a defined point or the center.
    John Jasa 2016-08
    """

    # Get the points
    coor = Geometry.coor

    # Set the central point as the origin if no point provided
    if not point:
        point = [0.0, 0.0, 0.0]

    # Convert point to array so that we can subtract it from the coordinates
    point = np.atleast_2d(point)

    # Get the relative coordinates to the central point
    relCoor = coor.copy() - point

    # Scale the relative coordinates and add them back to the central point
    # to obtain the new coor
    relCoor *= factor
    coor[:, :] = relCoor + point

    Geometry.coor = coor


def rotate(Geometry, angle, axis, point=None):
    """
    Rotate coor about an axis at a defined angle (in degrees).
    Axis should be either 0 for x, 1 for y, or 2 for z.

    John Jasa 2016-08
    """

    # Get the points
    coor = Geometry.coor

    # Assign rotation origin if the user provided None
    if point is None:
        point = [0.0, 0.0, 0.0]

    # Convert point to array so that we can subtract it from the coordinates
    point = np.atleast_2d(point)

    # Define the rotation angle in radians
    angle = angle * np.pi / 180.0
    rotationMat = np.zeros((3, 3))

    # Depending on the axis selected, define the indices of the matrix
    # where we will have the correct rotation functions.
    if axis == 0:
        ind1 = 1
        ind2 = 2
    elif axis == 1:
        ind1 = 2
        ind2 = 0
    else:
        ind1 = 0
        ind2 = 1

    # Set up the rotation matrix based on the axis selected
    rotationMat[ind1, ind1] = np.cos(angle)
    rotationMat[ind1, ind2] = -np.sin(angle)
    rotationMat[axis, axis] = 1
    rotationMat[ind2, ind1] = np.sin(angle)
    rotationMat[ind2, ind2] = np.cos(angle)

    # Obtain the relative coordinates to the central point
    relCoor = coor - point

    # Compute the rotation around the selected axis of the relative coordinates
    relCoor[:, :] = (rotationMat.dot(relCoor.T)).T

    # Obtain and output the rotated coordinates
    coor = relCoor + point

    Geometry.coor = coor


# ============================================================
# ============================================================
# ============================================================

# NORMALIZATION FUNCTION

# ORIGINAL FUNCTION
def normalize(vec):

    """
    This is a simple function to normalize an array of vectors.
    Each row of vec will be treated as a vector that should be normalized.

    INPUTS
    vec: array[n x m] -> Array of n vectors of size m that should be normalized

    Returns
    -------
    normalVec: array[n x m] -> Array of n normalized vectors

    Ney Secco 2016-10
    """

    # STEP 1: compute sum of squares
    vecNorms = np.array([np.sum(vec ** 2, axis=1)]).T

    # STEP 2: compute norms
    vecNorms = np.sqrt(vecNorms)

    # STEP 3: normalize vec
    normalVec = vec / vecNorms

    return normalVec


# FORWARD MODE
def normalize_d(vec, vecd):

    """
    This is the forward differentiation of the normalize routine

    Ney Secco 2016-10
    """

    # STEP 1_d
    vecNorms = np.array([np.sum(vec ** 2, axis=1)]).T
    vecNormsd = np.array([np.sum(2 * vec * vecd, axis=1)]).T

    # STEP 2_d
    vecNorms = np.sqrt(vecNorms)
    vecNormsd = 0.5 * vecNormsd / vecNorms

    # STEP 3_d
    normalVec = vec / vecNorms
    normalVecd = vecd / vecNorms - vec * vecNormsd / vecNorms ** 2

    return normalVec, normalVecd


# REVERSE MODE
def normalize_b(vec, normalVecb):

    """
    This is the backward differentiation of the normalize routine

    Ney Secco 2016-10
    """

    # STEP 1
    vecNorms1 = np.array([np.sum(vec ** 2, axis=1)]).T

    # STEP 2
    vecNorms2 = np.sqrt(vecNorms1)

    # STEP 3
    normalVec = vec / vecNorms2

    # STEP 3_b
    vecNorms2b = np.array([np.sum(-vec / vecNorms2 ** 2 * normalVecb, axis=1)]).T
    vecb = normalVecb / vecNorms2

    # STEP 2_b
    vecNorms1b = vecNorms2b / 2.0 / vecNorms2

    # STEP 1_b
    vecb = vecb + 2.0 * vec * vecNorms1b

    return vecb


### HELPER FUNCTIONS ###


def closest_node(guideCurve, curve):
    """ Find closest node from a list of node coordinates. """
    curve = np.asarray(curve)

    # Get number of points in the seed curve
    nPoints = curve.shape[0]

    # Initialize arrays to call projection function
    dist2 = np.ones(nPoints) * 1e10
    xyzProj = np.zeros((nPoints, 3))
    tangents = np.zeros((nPoints, 3))
    elemIDs = np.zeros((nPoints), dtype="int32")

    # Call projection function to find the distance between every node of the
    # seed curve to the guide curve
    guideCurve.project(curve, dist2, xyzProj, tangents, elemIDs)

    # Find closest point
    ind = np.argmin(dist2)

    return ind


# ===================================
# INTERSECTION COMPOSITION FUNCTIONS
# ===================================
"""
The next set of functions can be used to expedite the intersection treatment of
common types of intersections. They already comprise a series of operations such as
splits and remeshes.
"""


def airfoil_intersection(
    manager,
    intCurveName,
    LECurve,
    UpperTECurve,
    LowerTECurve,
    vertDir,
    numTENodes,
    numSkinNodes,
    TESpacing,
    LESpacing,
    mergedCurveName,
    optionsBodyMesh,
    optionsWingMesh,
    wingGeomName,
    exportFeatures=False,
):

    """
    This script works for blunt traling edges

    intCurveName: string -> Name of the intersection curve that should be
                            treated as a wing-body intersection.
                            This intersection curve should already be in the
                            manager. manager.intersect() returns curve names.

    LECurve: TSurf curve object -> Curve that follows the wing leading edge.

    UpperTECurve: TSurf curve object -> Curve that follows the wing upper trailing edge.

    LowerTECurve: TSurf curve object -> Curve that follows the wing upper trailing edge.

    mergedCurveName: string -> Name that will be assigned to the intersection curve and
                               its corresponding collar mesh.

    optionsBodyMesh: Dictionary -> Options to grow the body mesh. The user does not need
                                     to specify boundary conditions. This is done automatically.
                                     Here is an example of dictionary:
    optionsBodyMesh = {
        'dStart' : 0.01*2/fact,
        'numLayers' : int(48/2*fact)+1, #Should be /4
        'extension' : 1.8,
        'epsE0' : 3.5,
        'theta' : -0.5,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 0,
        'nuArea' : 0.16,
        'numAreaPasses' : 20,
        'sigmaSplay' : 0.3,
        'cMax' : 5.0,
        'ratioGuess' : 10.0,
    }

    optionsWingMesh: Dictionary -> Options to grow the wing mesh. The user does not need
                                     to specify boundary conditions. This is done automatically.
                                     Here is an example of dictionary:

    optionsWingMesh = {
        'dStart' : 0.01*2/fact,
        'numLayers' : int(48/2*fact)+1,
        'extension' : 2.2,
        'epsE0' : 12.5,
        'theta' : -0.8,
        'alphaP0' : 0.25,
        'numSmoothingPasses' : 4,
        'nuArea' : 0.16,
        'numAreaPasses' : 0,
        'sigmaSplay' : 0.3,
        'cMax' : 5.0,
        'ratioGuess' : 10.0,
        'guideCurves' : ['curve_te_low','curve_te_upp'],
        'remesh' : True
    }

    wingGeomName: string -> Name of the geometry component that represents the wing.
                            We need this information because pySurf may store either
                            the wing or the body as the first marching side for the
                            collar mesh generation.

    """

    # Get the names of the components that share the intersection
    parentGeomNames = manager.intCurves[intCurveName].extra_data["parentGeoms"]

    # Check if the wing component is the first one that will be marched on
    if parentGeomNames[0] == wingGeomName:
        wingIsFirst = True
        bodyGeomName = parentGeomNames[1]
    else:
        wingIsFirst = False
        bodyGeomName = parentGeomNames[0]

    # SPLIT

    # Split curves based on TE and LE curves
    optionsDict = {"splittingCurves": [LECurve, UpperTECurve, LowerTECurve]}

    splitCurveNames = manager.split_intCurve(intCurveName, optionsDict, criteria="curve")

    # REMESH

    # Find the highest vertical coordinate of the entire intersection (vertical position)
    # We also use the same loop to find the curve with the minimum number of nodes,
    # since this is probably the trailing edge

    # Initialize variables with dummy values
    maxCoor = -99999
    minNumNodes = 99999

    # Loop over every curve after the split
    for curveName in splitCurveNames:

        # Get pointer to the curve object
        curve = manager.intCurves[curveName]

        # Check maximum vertical position
        curr_maxCoor = np.max(curve.coor[:, vertDir])

        # Compare against the maximum found so far
        maxCoor = max(maxCoor, curr_maxCoor)

        # Get the number of nodes
        curr_numNodes = curve.numNodes

        # Compare against the minimum found so far
        minNumNodes = min(minNumNodes, curr_numNodes)

    # Now we can identify and remesh each curve properly
    for curveName in splitCurveNames:

        # Get pointer to the curve object
        curve = manager.intCurves[curveName]

        # The trailing edge curve probably has less nodes than the other ones
        if curve.numNodes == minNumNodes:

            # This is the trailing edge curve.
            # Just apply an uniform spacing
            optionsDict = {"nNewNodes": int(numTENodes)}
            TE_curveName = manager.remesh_intCurve(curveName, optionsDict)

        else:

            # We have an upper or lower skin curve.
            # First let's identify if the curve is defined from
            # LE to TE or vice-versa
            curveCoor = curve.get_points()
            deltaX = curveCoor[-1, 0] - curveCoor[0, 0]

            if deltaX > 0:
                LE_to_TE = True
            else:
                LE_to_TE = False

            # Compute the highest vertical coordinate of the curve
            curr_maxCoor = np.max(curve.coor[:, vertDir])

            # Now we can determine if we have upper or lower skin
            if curr_maxCoor < maxCoor:

                # We are at the lower skin

                if LE_to_TE:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": TESpacing,
                        "finalSpacing": LESpacing,
                    }

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": LESpacing,
                        "finalSpacing": TESpacing,
                    }

                    LS_curveName = manager.remesh_intCurve(curveName, optionsDict)

            else:

                # We are at the upper skin

                if LE_to_TE:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": TESpacing,
                        "finalSpacing": LESpacing,
                    }

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

                else:

                    optionsDict = {
                        "nNewNodes": numSkinNodes,
                        "spacing": "hypTan",
                        "initialSpacing": LESpacing,
                        "finalSpacing": TESpacing,
                    }

                    US_curveName = manager.remesh_intCurve(curveName, optionsDict)

    # MERGE
    curveNames = [TE_curveName, LS_curveName, US_curveName]
    manager.merge_intCurves(curveNames, mergedCurveName)

    # REORDER
    # The first node should be the leading edge node
    manager.intCurves[mergedCurveName].shift_end_nodes(criteria="curve", curveObject=LECurve)

    # Flip the curve for marching if necessary
    # Here we verify if the first two points (which are on the leading edge)
    # are going up or down. They should always go down so that the mesh
    # marhcing direction points towards the wing tip.
    mergedCurveCoor = manager.intCurves[mergedCurveName].get_points()

    deltaCoor = mergedCurveCoor[1, :] - mergedCurveCoor[0, :]

    _, bodyNormal, _ = manager.geoms[bodyGeomName].project_on_surface(np.array([mergedCurveCoor[0, :]]))

    bodyMarchDir = np.cross(bodyNormal[0, :], deltaCoor)

    # Make sure that the marching direction points to the wing tip if the wing
    # is the first body to be marched on.
    # Similarly, we can also allow the marching direction to point to the symmetry plane
    # if the body is the first one to be marched on.
    if (bodyMarchDir[0] < 0) and (wingIsFirst):
        manager.intCurves[mergedCurveName].flip()
    elif (bodyMarchDir[0] > 0) and (not wingIsFirst):
        manager.intCurves[mergedCurveName].flip()

    # Export final curve
    if exportFeatures:
        manager.intCurves[mergedCurveName].export_tecplot(mergedCurveName)

    # MARCH SURFACE MESHES

    # Add approriate boundary conditions
    optionsBodyMesh["bc1"] = "continuous"
    optionsBodyMesh["bc2"] = "continuous"

    optionsWingMesh["bc1"] = "curve:" + LECurve.name
    optionsWingMesh["bc2"] = "curve:" + LECurve.name
    optionsWingMesh["guideCurves"] = [LowerTECurve.name, UpperTECurve.name]
    optionsWingMesh["remesh"] = True

    # Now we need to verify which side of the collar mesh is marched first
    # to assign the correct options.
    if wingIsFirst:
        options0 = optionsWingMesh
        options1 = optionsBodyMesh
    else:
        options0 = optionsBodyMesh
        options1 = optionsWingMesh

    # Extrude mesh
    meshNames = manager.march_intCurve_surfaceMesh(
        mergedCurveName, options0=options0, options1=options1, meshName=mergedCurveName
    )

    # Export meshes if requested by the user
    if exportFeatures:
        for meshName in meshNames:
            manager.meshGenerators[meshName].export_plot3d(meshName + ".xyz")
