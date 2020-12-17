# plot3d_interface
# Written by:
# Ney Rafael Secco
# neysecco@umich.edu
# Jan 2016


import numpy as np

# DEFINING CLASSES

class Block():

    def __init__(self, X, Y, Z):

        # X, Y, and Z should be 3D arrays

        # Initialize
        self.shape = X.shape
        self.X = X
        self.Y = Y
        self.Z = Z

    def update(self, X, Y, Z):

        # Update coordinates and size

        # Initialize
        self.shape = X.shape
        self.X = X
        self.Y = Y
        self.Z = Z

class Grid():

    def __init__(self):
        self.blocks = []
        self.num_blocks = 0

    def add_block(self, X, Y, Z):
        # Create an block object with the given parameters
        dummy_block = Block(X, Y, Z)
        # append this block to the block list
        self.blocks.append(dummy_block)
        # Increase number of blocks
        self.num_blocks = self.num_blocks + 1

    def remove_block(self, blockID):
        # Remove this block to the block list
        self.blocks.pop(blockID)
        # Decrease number of blocks
        self.num_blocks = self.num_blocks - 1

    def convert_2dto3d(self,zSpan=1,num_i=2):

        '''
        If the size of the blocks are 1 x num_j x num_k (which is the case if you use the read_plot3d
        function for dimension = 2) then this function will convert these blocks to num_i x num_j x num_k
        by creating multiple identical slices of the 2D coordinates. The expansion will be done in
        the Z coordinate.
        '''

        # EXECUTION
        for iblock in range(self.num_blocks):

            # Expand the coordinate matrices
            X3d = np.array([np.copy(self.blocks[iblock].X[0,:,:]) for dummy in range(num_i)])
            Y3d = np.array([np.copy(self.blocks[iblock].Y[0,:,:]) for dummy in range(num_i)])

            # Create the span-wise coordinates
            Z3d = np.array([z*np.ones((self.blocks[iblock].shape[1:])) for z in np.linspace(0,zSpan,num_i)])

            # Update coordinates
            self.blocks[iblock].update(X3d, Y3d, Z3d)

    def flip_YZ(self):

        '''
        This function will flip the Y and Z coordinates
        '''

        for iblock in range(self.num_blocks):

            # Gather old coordinates
            X = self.blocks[iblock].X
            Y = self.blocks[iblock].Y
            Z = self.blocks[iblock].Z

            # Update values with flipped order
            self.blocks[iblock].update(X,Z,Y)

    def remove_curves(self):

        '''
        This method will remove unidimensional blocks.
        This is useful when preparing Plot3D inputs for pyHyp

        Ney Secco 2016-10
        '''

        # Initialize block counter
        blockID = 0

        # Loop until all blocks are checked
        while blockID < self.num_blocks:

            # Get shape of the current block
            shape = self.blocks[blockID].shape

            # Check if block is unidimensional by checking how many "1"
            # we have in its shape
            numDim = shape.count(1)

            if numDim > 1: # Block is unidimensional if we have more than 1 "1"
                # Remove block
                self.remove_block(blockID)
            else:
                # Check next block
                blockID = blockID + 1

class Curve():

    def __init__(self, coor, barsConn, curveName):
        self.coor = coor
        self.barsConn = barsConn
        self.name = curveName

    def export_tecplot(self, fileName='curve'):

        from . import tecplot_interface as ti

        ti.writeTecplotFEdata(self.coor, self.barsConn, self.name, fileName)

#===========================================

# DEFINE FUNCTIONS

#===========================================

def read_plot3d(fileName,dimension):

    '''
    The function will return a Grid object containg the coordinates.
    dimension should be 2 or 3 depending on the dimension of the plot3d file.
    '''

    # EXECUTION

    # Open file
    fid = open(fileName,'r')

    # Read all numbers and oppend to the same list
    data_list = []
    for line in fid:  #Line is a string
        # Split the string on whitespace, return a list of numbers (as strings)
        numbers_str = line.split()
        # Convert numbers to floats and append them to the list
        for x in numbers_str:
            data_list.append(float(x))

    # Close file
    fid.close()

    # Read the number of blocks
    num_blocks = int(data_list[0])

    # Read dimensions of each block. Each entry will have I,J,K sizes for each block
    shape_block = [list(map(int, data_list[dimension*(block_id)+1:dimension*(block_id+1)+1])) for block_id in range(num_blocks)]

    # Define index position from where we will start reading coordinates
    # NOTE: I define it as a list so I can change it inside the read_coord function
    reading_index = [dimension*(num_blocks)+1]

    # Define coordinate reader function
    def read_coord(num_i,num_j,num_k):
        # Initialize array with coordinates
        coord = np.zeros((num_i, num_j, num_k))
        # k is the outermost variable that should change
        for ind_k in range(num_k):
            # j should be the second variable to change
            for ind_j in range(num_j):
                # The last one that should change is i
                for ind_i in range(num_i):
                    # Read coordinate
                    coord[ind_i,ind_j,ind_k] = data_list[reading_index[0]]
                    # Shift index
                    reading_index[0] = reading_index[0] + 1
        # Return coordinates
        return coord

    # Create Grid object
    newGrid = Grid()

    # Read coordinates of each block using auxiliary function
    for iblock in range(num_blocks):

        # Read coordinates
        if dimension == 2:
            X = read_coord(shape_block[iblock][0],shape_block[iblock][1],1)
            Y = read_coord(shape_block[iblock][0],shape_block[iblock][1],1)
            Z = np.zeros(Y.shape)
        elif dimension == 3:
            X = read_coord(shape_block[iblock][0],shape_block[iblock][1],shape_block[iblock][2])
            Y = read_coord(shape_block[iblock][0],shape_block[iblock][1],shape_block[iblock][2])
            Z = read_coord(shape_block[iblock][0],shape_block[iblock][1],shape_block[iblock][2])

        # Store new block
        newGrid.add_block(X,Y,Z)

    # RETURNS
    return newGrid

#===========================================
#===========================================
#===========================================

def export_plot3d(grid, fileName, saveNumpy=False):

    '''
    grid should be a Grid object
    '''

    # Open file
    fid = open(fileName,'w')

    # Write number of blocks
    fid.write('%d\n'% grid.num_blocks)

    # Write dimensions of each block
    for iblock in range(grid.num_blocks):
        # Get the number of points
        block_shape = grid.blocks[iblock].X.shape
        # Write 3 integers corresponding to each of the block dimensions
        fid.write('%d %d %d\n'% (block_shape[0], block_shape[1], block_shape[2]))

    # Define coordinate writer function
    def write_coord(Coord):
        # k is the outermost variable that should change
        for ind_k in range(Coord.shape[2]):
            # j should be the second variable to change
            for ind_j in range(Coord.shape[1]):
                # The last one that should change is i
                for ind_i in range(Coord.shape[0]):
                    fid.write('%20.15g\n'% Coord[ind_i, ind_j, ind_k])

    # Initialize empty array of shape (0, 3). We will add coordinates of each
    # block to this array.
    if saveNumpy:
        fullCoords = np.zeros((0, 3))

    # Write coordinates of each block
    for iblock in range(grid.num_blocks):

        # Get current block
        block = grid.blocks[iblock]

        # Write coordinates
        write_coord(block.X)
        write_coord(block.Y)
        write_coord(block.Z)

        # Line break to start another block
        fid.write('\n')

        # Save this current block's coordinates in an array then update the
        # fullCoords with all information up to and including the current block.
        if saveNumpy:
            blockCoords = np.array([block.X.flatten(order='F'), block.Y.flatten(order='F'), block.Z.flatten(order='F')]).T
            fullCoords = np.vstack((fullCoords, blockCoords))

    # Actually save the file that contains all of the coordinate info.
    if saveNumpy:
        np.save(fileName[:-4], fullCoords)

    # Close file
    fid.close()

#===========================================
#===========================================
#===========================================

def merge_plot3d(files, flips, outputFile='merged.xyz'):

    '''
    This subroutine will merge two plot3d files into one. It may also
    flip the orientation of some meshes so that the final block remains
    structured.

    INPUTS:

    files: list of strings -> Name of one plot3d file to be merged.

    flips: list of boolean[3] -> If true, it will flip the node ordering along the
                                 the corresponding dimension for files. You should
                                 provide one list for every grid.

    outputFile: string -> Name of the output file (without extension).

    EXAMPLE:
    merge_plot3d(['a.xyz','b.xyz'], [[0,1,0],[0,0,0]])
    -- This will flip the j coordinates in a.xyz and then merge with b.xyz

    ATTENTION: I still need to implement a smart way to check shared faces
    to automatically concatenate. For now I will concatenate along the j
    direction.

    Ney Secco 2016-10
    '''

    # Print log
    print('Merging Plot3D files')

    # Initialize list of grids and number of blocks
    Grids = []
    numBlocks = []

    # Initialize merged grid
    mergedGrid = Grid()

    # Read grids
    for ii in range(len(files)):
        # Read input files and append to grid list
        Grids.append(read_plot3d(files[ii],3))

        # Print log
        print('Block dimensions from '+files[ii])
        for jj in range(Grids[ii].num_blocks):
            print(Grids[ii].blocks[jj].shape)

    # Flip coordinates along each dimension if necessary
    for problem in zip(Grids,flips):

        # Split Grid and flip from the tuple given by zip
        curGrid = problem[0]
        curFlip = problem[1]

        for dim in range(3):

            # Generate flip command as a string
            if dim == 0:
                flipCommand = '[::-1,:,:]'
            elif dim == 1:
                flipCommand = '[:,::-1,:]'
            elif dim == 2:
                flipCommand = '[:,:,::-1]'

            if curFlip[dim]: # Check if user requested to flip current dimension

                for ii in range(curGrid.num_blocks):

                    # Get pointers to each block coordinate
                    X = curGrid.blocks[ii].X
                    Y = curGrid.blocks[ii].Y
                    Z = curGrid.blocks[ii].Z

                    # Flip each coordinate
                    X[:,:,:] = eval('X'+flipCommand)
                    Y[:,:,:] = eval('Y'+flipCommand)
                    Z[:,:,:] = eval('Z'+flipCommand)

        # Now add blocks to the merged grid
        for ii in range(curGrid.num_blocks):

            # Get pointers to each block coordinate
            X = curGrid.blocks[ii].X
            Y = curGrid.blocks[ii].Y
            Z = curGrid.blocks[ii].Z

            # Add new block
            mergedGrid.add_block(X, Y, Z)

    # Now we can finally export the new grid
    export_plot3d(mergedGrid, outputFile)

#===========================================
#===========================================
#===========================================

def read_tecplot_curves(fileName):

    '''
    This function will read a curve definition from a Tecplot FE data file
    and assign coor and barsConn.

    Written by Ney Secco 2016-11
    '''

    # IMPORTS
    from . import tecplot_interface as ti

    # Read information from file
    sectionName, sectionData, sectionConn = ti.readTecplotFEdata(fileName)

    # Create curves for every section
    curves = []
    for secID in range(len(sectionName)):

        # Gather data
        # The -1 is to adjust connectivities to Python indexing, which starts at zero.
        coor = np.array(sectionData[secID])
        barsConn = np.array(sectionConn[secID],dtype='int32') - 1
        curveName = sectionName[secID]

        # Create curve object
        currCurve = Curve(coor, barsConn, curveName)

        # Append new curve to the dictionary
        curves.append(currCurve)

    # Return curves
    return curves
