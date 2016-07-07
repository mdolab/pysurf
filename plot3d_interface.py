# plot3d_interface
# Written by:
# Ney Rafael Secco
# neysecco@umich.edu
# Jan 2016

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

    def convert_2dto3d(self,zSpan=1,num_i=2):

        '''
        If the size of the blocks are 1 x num_j x num_k (which is the case if you use the read_plot3d
        function for dimension = 2) then this function will convert these blocks to num_i x num_j x num_k
        by creating multiple identical slices of the 2D coordinates. The expansion will be done in
        the Z coordinate.
        '''

        # IMPORTS
        from numpy import array, copy, ones, linspace

        # EXECUTION
        for iblock in range(self.num_blocks):

            # Expand the coordinate matrices
            X3d = array([copy(self.blocks[iblock].X[0,:,:]) for dummy in range(num_i)])
            Y3d = array([copy(self.blocks[iblock].Y[0,:,:]) for dummy in range(num_i)])

            # Create the span-wise coordinates
            Z3d = array([z*ones((self.blocks[iblock].shape[1:])) for z in linspace(0,zSpan,num_i)])

            # Update coordinates
            self.blocks[iblock].update(X3d, Y3d, Z3d)

    def flip_YZ(self):

        '''
        This function will flip the Y and Z coordinates
        '''

        for iblock in range(self.num_blocks):
            # Gather old coordiantes
            X = self.blocks[iblock].X
            Y = self.blocks[iblock].Y
            Z = self.blocks[iblock].Z
            # Update values with flipped order
            self.blocks[iblock].update(X,Z,Y)

#===========================================

# DEFINE FUNCTIONS

#===========================================

def read_plot3d(fileName,dimension):

    '''
    First you should initialize a grid object with the file name, then
    call the read_plot3d method!
    This function will always return 3D shaped arrays. In the case of 2d grids, the
    Z coordinate will be all zeroes.
    The function will return a Grid object containg the coordinates.
    dimension should be 2 or 3 depending on the dimension of the plot3d file.
    '''
    
    # IMPORTS
    from numpy import zeros
    
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
    shape_block = [map(int, data_list[dimension*(block_id)+1:dimension*(block_id+1)+1]) for block_id in range(num_blocks)]
    
    # Define index position from where we will start reading coordinates
    # NOTE: I define it as a list so I can change it inside the read_coord function
    reading_index = [dimension*(num_blocks)+1]

    # Define coordinate reader function
    def read_coord(num_i,num_j,num_k):
        # Initialize array with coordinates
        coord = zeros((num_i, num_j, num_k))
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
            X = read_coord(1,shape_block[iblock][0],shape_block[iblock][1])
            Y = read_coord(1,shape_block[iblock][0],shape_block[iblock][1])
            Z = zeros(Y.shape)
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

def export_plot3d(grid, fileName):

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

    # Close file
    fid.close()


