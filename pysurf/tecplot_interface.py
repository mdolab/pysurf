# IMPORTS
import numpy as np


def write_tecplot_scatter(filename, title, variable_names, data_points):

    # Open the data file
    fid = open(filename, "w")

    # Write the title
    fid.write("title = " + title + "\n")

    # Write tha variable names
    varnames_commas = ",".join(variable_names)  # Merge names in a single string separated by commas
    fid.write("variables = " + varnames_commas + ",\n")  # Write to the file

    # Write data points
    if type(data_points) is list:  # Check if user provided a list
        for point in data_points:
            str_points = [str(x) for x in point]  # Covert each entry to string
            str_points = " ".join(str_points)  # Merge all entries in a single string separated by whitespace
            fid.write(str_points + "\n")  # Write to file
    else:  # The user probably provided a numpy array
        for index in range(data_points.shape[0]):
            str_points = [str(x) for x in data_points[index, :]]
            str_points = " ".join(str_points)  # Merge all entries in a single string separated by whitespace
            fid.write(str_points + "\n")  # Write to file

    # Close file
    fid.close()


# ==================================================================
# ==================================================================


def readTecplotFEdata(fileName):

    # Written by Ney Secco
    # This script will read sections from a Tecplot FE file
    # FOR NOW THIS JUST READS BAR FE FILES FOR THE CURVE FUNCTIONS

    # IMPORTS
    from numpy import array

    # Initialize list that will hold node and connectivity data for each slice
    sectionName = []
    sectionData = []
    sectionConn = []

    # Initialize counters
    sectionID = -1

    # Read the lines of the file
    with open(fileName, "r") as fid:
        lines = fid.readlines()

    # Loop over the lines to gather data
    for line in lines:

        # Split data
        data = line.split()

        # Adding something to data if we got an empty list so that the
        # rest of the code does not break
        if len(data) == 0:
            data = [None]

        # Detect if we reached a new section definition
        if data[0].lower() == "zone":

            # Gather name of the new zone
            zoneName = line.split('"')[1]

            # Save zone name
            sectionName.append(zoneName)

            # Save data from previous section
            if sectionID >= 0:  # Only if we have a previous section...
                sectionData.append(array(currData))  # noqa: F821
                sectionConn.append(array(currConn))  # noqa: F821

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
                currData.append(list(map(float, data)))
            except Exception:
                pass

        # Check if we have a connectivity line
        if len(data) == 2:

            # Try to convert the elements of this line to numbers.
            # If this doesn't work, it means we do not have a connectivity line.
            try:
                currConn.append(list(map(int, data)))
            except Exception:
                pass

    # Append the last section
    sectionData.append(array(currData))
    sectionConn.append(array(currConn))

    # RETURNS
    return sectionName, sectionData, sectionConn


def readTecplotCurves(fileName):

    """
    This function will read a curve definition from a Tecplot FE data file
    and assign coor and barsConn.

    Written by Ney Secco 2016-11
    """

    # Read information from file
    sectionName, sectionData, sectionConn = readTecplotFEdata(fileName)

    # Create curves for every section
    curves = []
    for secID in range(len(sectionName)):

        # Gather data
        # The -1 is to adjust connectivities to Python indexing, which starts at zero.
        coor = np.array(sectionData[secID])
        barsConn = np.array(sectionConn[secID], dtype="int32") - 1
        curveName = sectionName[secID]

        # Create curve object
        currCurve = Curve(coor, barsConn, curveName)

        # Append new curve to the dictionary
        curves.append(currCurve)

    # Return curves
    return curves


class Curve:
    def __init__(self, coor, barsConn, curveName):
        self.coor = coor
        self.barsConn = barsConn
        self.name = curveName

    def export_tecplot(self, fileName="curve"):

        writeTecplotFEdata(self.coor, self.barsConn, self.name, fileName)


# ==================================================================
# ==================================================================


def readTecplotFEdataSurf(fileName):

    # Written by Ney Secco
    # This script will read sections from a Tecplot FE file with quads written by pysurf.
    # The Tecplot file should have a single zone

    # Initialize arrays
    coor = []
    triaConn = []
    quadsConn = []

    # Read the lines of the file
    with open(fileName, "r") as fid:
        lines = fid.readlines()

    # Loop over the lines to gather data
    for line in lines:

        # Split data
        data = line.split()

        # Assing something to data if we got an empty list so that the
        # rest of the code does not break.
        if len(data) == 0:
            data = [None]

        # Check if we have a coordinate line
        if len(data) == 3:

            # Try to convert the elements of this line to numbers.
            # If this doesn't work, it means we do not have a data line.
            try:
                coor.append(list(map(float, data)))
            except Exception:
                pass

        # Check if we have a quad connectivity line
        if len(data) == 4:

            # Try to convert the elements of this line to numbers.
            # If this doesn't work, it means we do not have a connectivity line.
            try:
                # We also do set() to eliminate equal entries, so we can detect triangles as degenerate quads
                currConn = list(set(map(int, data)))

                # See if we have a tria or a quad
                if len(currConn) == 3:
                    triaConn.append(currConn)
                elif len(currConn) == 4:
                    quadsConn.append(currConn)
            except Exception:
                pass

    # Transform everything to numpy arrays
    # The -1 is to make it consistent with python indexing
    coor = np.array(coor)
    triaConn = np.array(triaConn) - 1
    quadsConn = np.array(quadsConn) - 1

    # RETURNS
    return coor, triaConn, quadsConn


# ==================================================================
# ==================================================================


def writeTecplotFEdata(coor, barsConn, curveName, fileName):

    # This script will write sections from a Tecplot FE file
    # Written by John Hwang. Adapted by Ney Secco.

    # Add extension to the file name
    fileName = fileName + ".plt"

    # Open file
    fileID = open(fileName, "w")

    # Add title
    fileID.write('Title = "TSurf curve FE data" \n')

    # Add variable names
    fileID.write("Variables = ")
    var_names = ["X", "Y", "Z"]
    for name in var_names:
        fileID.write('"Coordinate' + name + '" ')
    fileID.write("\n")

    # Gather number of nodes and finite elements
    nNodes = coor.shape[0]
    nBars = barsConn.shape[0]

    # Write curve data
    fileID.write('Zone T= "' + curveName + '"\n')
    fileID.write("Nodes=" + str(nNodes) + ", Elements=" + str(nBars) + ", ZONETYPE=FELineSeg\n")
    fileID.write("DATAPACKING=POINT\n")

    # Write nodal coordinates
    np.savetxt(fileID, coor)

    # Write connectivities
    np.savetxt(fileID, barsConn + 1, fmt="%i")

    # Close output file
    fileID.close()

    # Print log
    print("Curve " + curveName + " saved to file " + fileName)


# ==================================================================
# ==================================================================


def writeTecplotSurfaceFEData(coor, triaConn, quadsConn, surfName, fileName):

    """
    This method will export the triangulated surface data into a tecplot format.
    The will write all elements as quad data. The triangle elements will be exported
    as degenerated quads (quads with a repeated node).

    Ney Secco 2017-02
    """

    # Add extension to the file name
    fileName = fileName + ".plt"

    # Open file
    fileID = open(fileName, "w")

    # Add title
    fileID.write('Title = "TSurf Surface FE data" \n')

    # Add variable names
    fileID.write("Variables = ")
    var_names = ["X", "Y", "Z"]
    for name in var_names:
        fileID.write('"Coordinate' + name + '" ')
    fileID.write("\n")

    # Gather number of nodes and finite elements
    nNodes = coor.shape[0]
    nTria = triaConn.shape[0]
    nQuads = quadsConn.shape[0]

    # Write curve data
    fileID.write('Zone T= "' + surfName + '"\n')
    fileID.write("Nodes=" + str(nNodes) + ", Elements=" + str(nTria + nQuads) + ", ZONETYPE=FEQUADRILATERAL\n")
    fileID.write("DATAPACKING=POINT\n")

    # Write nodal coordinates
    np.savetxt(fileID, coor)

    # Extend tria connectivities to include degenerate node
    triaConnExt = np.hstack([triaConn, np.array([triaConn[:, -1]]).T])

    # Write tria connectivities
    np.savetxt(fileID, triaConnExt + 1, fmt="%i")

    # Write quad connectivities
    np.savetxt(fileID, quadsConn + 1, fmt="%i")

    # Close output file
    fileID.close()

    # Print log
    print("Surface " + surfName + " saved to file " + fileName)
