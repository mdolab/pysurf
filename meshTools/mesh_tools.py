'''
This script contains general functions to operate on
Mesh objects.
'''

# IMPORTS
from __future__ import division
import numpy as np
import pysurf
from collections import OrderedDict

#=====================================
# MERGING FUNCTIONS

def merge_meshes(mergedMeshName, meshList):

    '''
    This function will combine all blocks from all meshes in
    meshList into a single mesh object.

    We rename all blocks using a new standard to avoid any duplicate issues.

    Ney Secco 2017-02
    '''

    # Initialize domain counter
    domainID = 1

    # Initialize dictionary that will gather all blocks
    coor = OrderedDict()

    for mesh in meshList:

        # Loop over all the blocks since we have to rename them anyway
        for blockName in mesh.coor:

            # Generate name of the new block
            mergedBlockName = 'domain.%05d'%domainID

            # Copy coordinates
            coor[mergedBlockName] = mesh.coor[blockName].copy()

            # Increment domain counter
            domainID = domainID + 1

    # Initialize new mesh object
    mergedMesh = pysurf.SurfaceMesh(mergedMeshName, coor)

    # Return the merged mesh
    return mergedMesh
