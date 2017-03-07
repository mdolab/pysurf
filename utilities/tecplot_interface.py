# IMPORTS
import numpy as np

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

#==================================================================
#==================================================================

def readTecplotFEdata(fileName):

    # Written by Ney Secco
    # This script will read sections from a Tecplot FE file

    # IMPORTS
    from numpy import array

    # Initialize list that will hold node and connectivity data for each slice
    sectionName = []
    sectionData = []
    sectionConn = []

    # Initialize counters
    sectionID = -1

    # Read the lines of the file
    with open(fileName,'r') as fid:
        lines = fid.readlines()

    # Loop over the lines to gather data
    for line in lines:

        # Split data
        data = line.split()

        # Assing something to data if we got an empty list so that the
        # rest of the code does no break.
        if len(data) == 0:
            data = [None]

        # Detect if we reached a new section definition
        if data[0].lower() == 'zone':

            # Gather name of the new zone
            zoneName = line.split('"')[1]

            # Save zone name
            sectionName.append(zoneName)

            # Save data from previous section
            if sectionID >= 0: # Only if we have a previous section...
                sectionData.append(array(currData))
                sectionConn.append(array(currConn))

            # Increment section counter
            sectionID = sectionID + 1

            # Initialize lists that will hold the current section data
            currData = []
            currConn = []

        # Check if we have a data line
        if len(data) == 3:

            # Try to convert the elements of this line to numbers.
            # If this doesn't work, it means we do not have a data line.
            try:
                currData.append(map(float,data))
            except:
                pass

        # Check if we have a connectivity line
        if len(data) == 2:

            # Try to convert the elements of this line to numbers.
            # If this doesn't work, it means we do not have a connectivity line.
            try:
                currConn.append(map(int,data))
            except:
                pass

    # Append the last section
    sectionData.append(array(currData))
    sectionConn.append(array(currConn))

    # RETURNS
    return sectionName, sectionData, sectionConn

#==================================================================
#==================================================================

def writeTecplotFEdata(coor,barsConn,curveName,fileName):

    # This script will write sections from a Tecplot FE file
    # Written by John Hwang. Adapted by Ney Secco.

    # Add extension to the file name
    fileName = fileName + '.plt'

    # Open file
    fileID = open(fileName, 'w')

    # Add title
    fileID.write('Title = \"TSurf curve FE data\" \n')

    # Add variable names
    fileID.write('Variables = ')
    var_names = ['X', 'Y', 'Z']
    for name in var_names:
        fileID.write('\"Coordinate' + name + '\" ')
    fileID.write('\n')
    
    # Gather number of nodes and finite elements
    nNodes = coor.shape[0]
    nBars = barsConn.shape[0]

    # Write curve data
    fileID.write('Zone T= \"'+curveName+'\"\n')
    fileID.write('Nodes=' + str(nNodes) + ', Elements=' + str(nBars) + ', ZONETYPE=FELineSeg\n')
    fileID.write('DATAPACKING=POINT\n')
    
    # Write nodal coordinates
    np.savetxt(fileID, coor)

    # Write connectivities
    np.savetxt(fileID, barsConn+1, fmt='%i')

    # Close output file
    fileID.close()

    # Print log
    print 'Curve '+curveName+' saved to file '+fileName

#==================================================================
#==================================================================

def writeTecplotSurfaceFEData(coor,triaConn,quadsConn,surfName,fileName):

    '''
    This method will export the triangulated surface data into a tecplot format.
    The will write all elements as quad data. The triangle elements will be exported
    as degenerated quads (quads with a repeated node).
    
    Ney Secco 2017-02
    '''

    # Add extension to the file name
    fileName = fileName + '.plt'

    # Open file
    fileID = open(fileName, 'w')

    # Add title
    fileID.write('Title = \"TSurf Surface FE data\" \n')

    # Add variable names
    fileID.write('Variables = ')
    var_names = ['X', 'Y', 'Z']
    for name in var_names:
        fileID.write('\"Coordinate' + name + '\" ')
    fileID.write('\n')
    
    # Gather number of nodes and finite elements
    nNodes = coor.shape[1]
    nTria = triaConn.shape[1]
    nQuads = quadsConn.shape[1]

    # Write curve data
    fileID.write('Zone T= \"'+surfName+'\"\n')
    fileID.write('Nodes=' + str(nNodes) + ', Elements=' + str(nTria+nQuads) + ', ZONETYPE=FEQUADRILATERAL\n')
    fileID.write('DATAPACKING=POINT\n')
    
    # Write nodal coordinates
    np.savetxt(fileID, coor.T)

    # Extend tria connectivities to include degenerate node
    triaConnExt = np.vstack([triaConn, triaConn[-1,:]])

    # Write tria connectivities
    np.savetxt(fileID, triaConnExt.T, fmt='%i')

    # Write quad connectivities
    np.savetxt(fileID, quadsConn.T, fmt='%i')

    # Close output file
    fileID.close()

    # Print log
    print 'Surface '+surfName+' saved to file '+fileName

#==================================================================
#==================================================================

def convert(oldData, oldConn):

    # Written by Gaetan Kenway
    # This script will organize Tecplots FE data
    # The data output will be ordered, an will go from LE,
    # over the top skin to TE, then follow the lower skin back to LE.

    import numpy
    import pygeo
    import sys
    import os
    import re

    X = numpy.zeros((len(oldData), 3))
    X = oldData[:, 0:3]

    uniquePts, link = pygeo.geo_utils.pointReduce(X, nodeTol=1e-12)
    nUnique = len(uniquePts)

    # Create the mask for the unique data:
    mask = numpy.zeros(nUnique, 'intc')
    for i in range(len(link)):
        mask[link[i]] = i

    # De-duplicate the data
    data = oldData[mask, :]

    # Update the conn for the new data...remember conn is 1 based
    conn = numpy.zeros_like(oldConn)
    for i in range(len(oldConn)):
        conn[i, 0] = link[oldConn[i, 0]-1]
        conn[i, 1] = link[oldConn[i, 1]-1]

    # Create the node-to-elem structure for easier searching
    nodeToElem = -1*numpy.ones((nUnique, 2))

    for i in range(len(conn)):
        for jj in range(2):
            n = conn[i, jj]
            if nodeToElem[n ,0] == -1:
                nodeToElem[n, 0] = i
            elif nodeToElem[n, 1] == -1:
                nodeToElem[n, 1] = i
            else:
                pass
                # This should probably be an error...

    # Find the minimum x-coordinate
    iMin = numpy.argmin(data[:, 0])

    # This is the "reorder mask" the ordering that will make everything
    # I-ordered
    rMask = -1*numpy.ones(nUnique, 'intc')

    # First point is min x.
    rMask[0] = iMin

    # Now we need to determine which way to go: Get the two nodes
    # connected to iMin

    elem = nodeToElem[iMin, 0]
    nodes = conn[elem]
    if nodes[0] == iMin:
        n1 = nodes[1]
    else:
        n1 = nodes[0]

    elem = nodeToElem[iMin, 1]
    nodes = conn[elem]
    if nodes[0] == iMin:
        n2 = nodes[1]
    else:
        n2 = nodes[0]

    # Now Go with the node that is higher in 'z'
    if data[n1, 2] > data[iMin, 2]:
        nextNode = n1
    elif data[n2, 2] > data[iMin, 2]:
        nextNode = n2
    else:
        print 'Something wierd happened'
        print data[iMin, 2], data[n1, 2], data[n2, 2]
        sys.exit(0)

    # So now we know the first and second nodes of rMask
    rMask[1] = nextNode

    # We have two nodes done
    i = 2

    while i < nUnique:

        # Find the index of the i-th value in rMask. We can guarantee that
        # i-1 and i-2 entries exist. 

        curNode = rMask[i-1]
        oldNode = rMask[i-2]

        for jj in range(2):
            elem = nodeToElem[curNode, jj]
            nodes = conn[elem]
            for kk in range(2):
                if nodes[kk] != curNode and nodes[kk] != oldNode:
                    nextNode = nodes[kk]

        rMask[i] = nextNode
        i += 1

    # Finally we can reorder the data for a final time:
    data = data[rMask, :]

    return data
