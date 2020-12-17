
import numpy as np
from mpi4py import MPI
import os
import pysurf
from collections import OrderedDict

class SurfaceMesh(object):

    #=================================================
    # BASE METHODS

    def __init__(self, meshName, arg):
        '''
        This function initializes a multi-block surface mesh object.

        INPUTS:

        meshName: string with name of the new mesh object

        arg: could be a dictionary or a string

        In case it is a dictionary, the keys should be block names of
        the structured surface mesh, and the values should be 3D arrays
        of size [imax,jmax,3] containing the nodal coordinates of the
        corresponding blocks.

        In case it is a string, it should specify the file name that
        contains the surface coordinates. If the file extension is
        ".cgns" we will use the CGNS reader. Otherwise, we will treat it
        as a formatted Plot3D file.

        CLASS ATTRIBUTES:

        self.coor : dictionary -> the keys are block names of
        the structured surface mesh, and the values are 3D arrays
        of size [imax,jmax,3] containing the nodal coordinates of the
        corresponding blocks.

        self.coord : dictionary -> the keys are block names of
        the structured surface mesh, and the values are 3D arrays
        of size [imax,jmax,3] containing the forward AD seeds of the
        corresponding nodal coordinates.

        self.coorb : dictionary -> the keys are block names of
        the structured surface mesh, and the values are 3D arrays
        of size [imax,jmax,3] containing the forward AD seeds of the
        corresponding nodal coordinates.

        self.numPts : integer -> Total number of points in this mesh
        '''
        

        # Assing name
        self.name = meshName

        # Figure out what we need to do based on the argument type

        if isinstance(arg, dict):

            # If the user gave a dictionary, we can assign it directly to
            # the mesh, but we will check the shape of the arrays

            for block in arg:

                dim = arg[block].shape[-1]

                if dim != 3:
                    raise ValueError('The coordinate array should be [imax,jmax,3]')

            # We can assign the data if it passes the test
            self.coor = OrderedDict(arg)

        elif isinstance(arg,str):

            # We received a filename
            fileName = arg

            # The user probably specified a file name. Let's detect if it
            # is a CGNS file, otherwise we will treat it as a Plot3D file
            extension = os.path.splitext(fileName)[1]            

            if extension == '.cgns':

                # Use the CGNS reader to get all coordinates
                self.coor, cellDim = pysurf.structCGNSreader.readFile(fileName)

                # Stop everything if we get a volume mesh rather than a surface mesh
                if cellDim != 2:
                    raise ValueError('Got a volume mesh instead of a surface mesh.')

            else: # Assume we have a plot3d file

                # Get the blocks from the plot3d file
                grid = pysurf.plot3d_interface.read_plot3d(fileName,3)
            
                # Initialize block dictionary
                self.coor = OrderedDict()

                # Rearrange block information so that if follows the standard ordering
                for blockID in range(grid.num_blocks):

                    # Generate block name
                    blockName = 'domain.%05d'%(blockID+1)

                    # Initialize 3D array to hold coordinates
                    blockCoor = np.zeros([grid.blocks[blockID].X.shape[0],
                                          grid.blocks[blockID].X.shape[1],
                                          3])

                    # Assign coordinates
                    blockCoor[:,:,0] = grid.blocks[blockID].X[:,:,0]
                    blockCoor[:,:,1] = grid.blocks[blockID].Y[:,:,0]
                    blockCoor[:,:,2] = grid.blocks[blockID].Z[:,:,0]
                    
                    # Add dictionary entry
                    self.coor[blockName] = blockCoor
        
        # Now we need to initialize arrays to hold derivative seeds.
        # Use the same loop to get the total number of points in this mesh
        self.coord = {}
        self.coorb = {}
        numPts = 0

        for block in self.coor:
            blockShape = self.coor[block].shape

            self.coord[block] = np.zeros(blockShape)
            self.coorb[block] = np.zeros(blockShape)

            numPts = numPts + blockShape[0]*blockShape[1]

        # Assign number of points
        self.numPts = numPts

    def __str__(self):

        '''
        This method gives the "print" behavior for a mesh object.
        If you call print mesh, this method will print relevant
        mesh info.
        '''

        # Blank space
        print('')

        # Print mesh name
        print('Printing log of mesh:',self.name)

        # Print the number of blocks
        print('Number of blocks:',len(self.coor))

        # Print block dimensions
        print('Block dimensions:')
        for block in self.coor:
            dims = self.coor[block].shape
            print(block,':',dims)

        # Blank space
        print('')

        # Return a blank string
        return ''

    #=================================================
    # NODAL COORDINATES METHODS

    def get_points(self):

        '''
        This function returns a flattened array with all nodal coordinates
        '''

        # Initialize array
        flatArray = np.array([[0,0,0]])

        # Append each point of every block
        for block in self.coor:

            # We need the order='F' so that we add all j=1 points first,
            # and then proceed to j=2, ...
            # Note that this does not change the actual memory layout
            # but just the reading order.
            currArray = self.coor[block].reshape((-1,3),order='F')

            flatArray = np.vstack([flatArray,currArray])

        # Remove the first dummy point
        flatArray = flatArray[1:,:]

        # Return the coordinate array
        return flatArray

    def set_points(self, flatArray):

        '''
        This function receives a flattened array with all nodal coordinates, and
        updates all block coordinates with it.
        '''

        # Define offset to slice the flattened array
        offset = 0

        # Append each point of every block
        for block in self.coor:

            # Get block dimensions
            dims = self.coor[block].shape

            # Compute how many points we have in this block
            numPoints = dims[0]*dims[1]

            # Assign the new points
            self.coor[block][:,:,:] = flatArray[offset:offset+numPoints,:].reshape((dims[0],dims[1],3),order='F')

            # Update the offset
            offset = offset + numPoints

    #=================================================
    # FORWARD DERIVATIVE SEEDS METHODS

    def get_forwardADSeeds(self):

        '''
        This function returns a flattened array with all nodal coordinate derivatives
        '''

        # Initialize array
        flatArray = np.array([[0,0,0]])

        # Append each point of every block
        for block in self.coor:

            currArray = self.coord[block].reshape((-1,3),order='F')

            flatArray = np.vstack([flatArray,currArray])

        # Remove the first dummy point
        flatArray = flatArray[1:,:]

        # Return the coordinate array
        return flatArray

    def set_forwardADSeeds(self, flatArray):

        '''
        This function receives a flattened array with all nodal coordinate derivatives, and
        updates all block coordinates with it.
        '''

        # Define offset to slice the flattened array
        offset = 0

        # Append each point of every block
        for block in self.coor:

            # Get block dimensions
            dims = self.coor[block].shape

            # Compute how many points we have in this block
            numPoints = dims[0]*dims[1]

            # Assign the new points
            self.coord[block] = flatArray[offset:offset+numPoints,:].reshape((dims[0],dims[1],3),order='F')

            # Update the offset
            offset = offset + numPoints

    #=================================================
    # REVERSE DERIVATIVE SEEDS METHODS

    def get_reverseADSeeds(self, clean=True):

        '''
        This function returns a flattened array with all nodal coordinate derivatives
        '''

        # Initialize array
        flatArray = np.array([[0,0,0]])

        # Append each point of every block
        for block in self.coor:

            currArray = self.coorb[block].reshape((-1,3),order='F')

            flatArray = np.vstack([flatArray,currArray])

            # Clean seeds if necessary
            if clean:
                self.coorb[block][:,:,:] = 0.0

        # Remove the first dummy point
        flatArray = flatArray[1:,:]

        # Return the coordinate array
        return flatArray

    def set_reverseADSeeds(self, flatArray):

        '''
        This function receives a flattened array with all nodal coordinate derivatives, and
        updates all block coordinates with it.
        '''

        # Define offset to slice the flattened array
        offset = 0

        # Append each point of every block
        for block in self.coor:

            # Get block dimensions
            dims = self.coor[block].shape

            # Compute how many points we have in this block
            numPoints = dims[0]*dims[1]

            # Assign the derivatives
            self.coorb[block] = flatArray[offset:offset+numPoints,:].reshape((dims[0],dims[1],3),order='F')

            # Update the offset
            offset = offset + numPoints

    def accumulate_reverseADSeeds(self, flatArray):

        '''
        This function receives a flattened array with all nodal coordinate derivatives, and
        updates all block coordinates with it.
        '''

        # Define offset to slice the flattened array
        offset = 0

        # Append each point of every block
        for block in self.coor:

            # Get block dimensions
            dims = self.coor[block].shape

            # Compute how many points we have in this block
            numPoints = dims[0]*dims[1]

            # Accumulate the new seeds
            self.coorb[block] = self.coorb[block] + flatArray[offset:offset+numPoints,:].reshape((dims[0],dims[1],3),order='F')

            # Update the offset
            offset = offset + numPoints

    def clean_reverseADSeeds(self):

        '''
        This function erases all reverse derivative seeds.
        '''

        # Loop over all blocks
        for block in self.coor:

            # Erase seeds
            self.coorb[block][:,:,:] = 0.0

    #=================================================
    # GENERAL DERIVATIVE SEEDS METHODS

    def set_randomADSeeds(self, mode='both', fixedSeed=True):

        # See if we should use a fixed seed for the RNG
        if fixedSeed:
            np.random.seed(123)

        # Set forward AD seeds
        if mode=='forward' or mode=='both':

            for block in self.coor:

                coord = np.array(np.random.rand(self.coor[block].shape[0],
                                                self.coor[block].shape[1],
                                                self.coor[block].shape[2]))
                coord = coord/np.sqrt(np.sum(coord**2))
                self.coord[block] = coord

        # Set reverse AD seeds
        if mode=='reverse' or mode=='both':

            for block in self.coor:

                coorb = np.array(np.random.rand(self.coor[block].shape[0],
                                                self.coor[block].shape[1],
                                                self.coor[block].shape[2]))
                coorb = coorb/np.sqrt(np.sum(coorb**2))
                self.coorb[block] = coorb


        # Return the randomly generated seeds

        if mode == 'forward':
            return self.get_forwardADSeeds()

        elif mode == 'reverse':
            return self.get_reverseADSeeds(clean=False)

        elif mode == 'both':
            return self.get_forwardADSeeds(), self.get_reverseADSeeds(clean=False)

    #=================================================
    # VISUALIZATION METHODS

    def export_plot3d(self, fileName):

        # Initialize Plot3D grid object
        myGrid = pysurf.plot3d_interface.Grid()

        # Add each block to the grid
        for block in self.coor:
            myGrid.add_block(self.coor[block][:,:,0:1],
                             self.coor[block][:,:,1:2],
                             self.coor[block][:,:,2:3])

        # Export grid
        pysurf.plot3d_interface.export_plot3d(myGrid, fileName)
