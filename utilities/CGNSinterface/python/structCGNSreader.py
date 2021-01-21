"""
This module contains basic function to read structured
CGNS files

Ney Secco 2017-02
"""

# IMPORTS

import numpy as np
from . import structCGNS_API
from collections import OrderedDict

# ===========================================
# READING FUNCTIONS


def readFile(fileName):

    """
    This function reads a CGNS file and outputs a mesh object
    with the nodal coordinates of all blocks.

    Parameters
    ----------

    fileName: string -> Name of the CGNS file to be read

    Returns
    -------

    Blocks: dictionary -> Dictionary whose keys are block names, and
    whose fields are the nodal coordinates. If the function read a
    3D (Volume) mesh, then the fields are arrays of size [imax,jmax,kmax,3].
    If we have a 2D (Surface) mesh instead, then the fields are arrays
    of size [imax,jmax,3].
    """

    # Print log
    print("")
    print("Reading structured CGNS file:", fileName)

    # Open the CGNS file, and receive the file handle
    cg = structCGNS_API.utils.openfile(fileName, 0)

    # Determine the dimension of the grid
    cellDim = structCGNS_API.utils.getgriddimension(cg)

    # Get the number of blocks
    numBlocks = structCGNS_API.utils.getnblocks(cg)

    # Initialize dictionary to store block information
    Blocks = OrderedDict()

    # Loop over the blocks to gather information.
    # We need to start at 1 due to CGNS ordering.
    for blockID in range(1, numBlocks + 1):

        # Get blcok information from the CGNS file
        blockName, dims, nboco, nb2b = structCGNS_API.utils.getblockinfo(cg, blockID)

        # Remove whitespaces from the block name
        blockName = blockName.strip()

        # Get block coordinates and assign them to the dictionary.
        # We need to call the corresponding version depending on the
        # mesh dimension
        if cellDim == 2:
            Blocks[blockName] = structCGNS_API.utils.getcoordinates2d(cg, blockID, dims[0], dims[1])

        elif cellDim == 3:
            Blocks[blockName] = structCGNS_API.utils.getcoordinates3d(cg, blockID, dims[0], dims[1], dims[2])

    # Print log
    print("Mesh dimension:", cellDim)
    print("Number of blocks:", len(Blocks))

    # Return the nodal coordinates separated by blocks
    return Blocks, cellDim
