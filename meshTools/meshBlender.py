# GENERAL IMPORT
import numpy as np

class collarBlender():

    def __init__(self, geom, bDec=-5.0, bExp=2.0):

        '''
        This function blends the specified collar mesh
        '''

        # Store references
        self.geom = geom
        self.hasRefCollar = False
        self.bDec = bDec
        self.bExp = bExp

    def project_collar(self, collarMeshObj):

        '''
        This function will find the parametric coordinates of the
        collar mesh nodes on the triangulated surface.
        This will update self.projDict.
        This also recomputes the blending weights.
        '''

        # Flatten the collar mesh nodes
        xyz = collarMeshObj.get_points()

        # Project the collar mesh nodes
        # The parametric coordinates are stored in projDict
        xyzProj, normProj, projDict = self.geom.project_on_surface(xyz)
        
        # Save projection information
        self.projDict = projDict

        # Compute and store blending weights
        self.weights = build_weight_array(collarMeshObj, self.bDec, self.bExp)

        # Now we have a reference mesh
        self.hasRefCollar = True

    def get_warped_collar(self):

        '''
        This will use the parametric coordinates obtained in
        self.project_collar to obtain physical coordinates on the
        current situation of self.geom.
        '''

        # Get important results of the projection
        elemType = self.projDict['elementType']
        elemID = self.projDict['elementID']
        uvw = self.projDict['uvw']

        # We need to adjust the elemID because ADT gives negative ID values
        # if the sampling node is not exactly on the element
        elemID = abs(elemID+1)-1

        # Initialize array of warped nodes
        warpedNodes = np.zeros((len(elemType),3))

        # Compute the physical coordinates of each collar node
        for nodeID in range(len(elemType)):

            currElemType = elemType[nodeID]
            currElemID = elemID[nodeID]
            u = uvw[nodeID,0]
            v = uvw[nodeID,1]

            if currElemType == 1: # We have a triangle

                node1 = self.geom.triaConnF[currElemID,0]-1
                node2 = self.geom.triaConnF[currElemID,1]-1
                node3 = self.geom.triaConnF[currElemID,2]-1
                X1 = self.geom.coor[node1,:]
                X2 = self.geom.coor[node2,:]
                X3 = self.geom.coor[node3,:]
                warpedNodes[nodeID,:] = get_tria_proj(X1,X2,X3,u,v)

            elif currElemType == 2: # We have a quad

                node1 = self.geom.quadsConnF[currElemID,0]-1
                node2 = self.geom.quadsConnF[currElemID,1]-1
                node3 = self.geom.quadsConnF[currElemID,2]-1
                node4 = self.geom.quadsConnF[currElemID,3]-1
                X1 = self.geom.coor[node1,:]
                X2 = self.geom.coor[node2,:]
                X3 = self.geom.coor[node3,:]
                X4 = self.geom.coor[node4,:]
                warpedNodes[nodeID,:] = get_quad_proj(X1,X2,X3,X4,u,v)

            if currElemID < 0:
                warpedNodes[nodeID,:] = 0.0

        # Return warped nodes
        return warpedNodes

    def blend(self, collarMeshObj):

        # Set the collar mesh as current reference if this blended still has None:
        if self.hasRefCollar is False:

            # Do the first projection to get parametric coordinates
            # for the collar mesh nodes on the triangulated surface.
            self.project_collar(collarMeshObj)

        else: # We already have a reference collar mesh to do the blending

            # Get regenerated mesh nodes
            regenNodes = collarMeshObj.get_points()

            # Get warped nodes (that have no new intersection information)
            warpedNodes = self.get_warped_collar()

            # Get weights
            weights = self.weights

            # Compute blended coordinates
            blendedNodes = weights*regenNodes + (1.0-weights)*warpedNodes

            #Project the blended nodes to the triangulated surface
            projBlendedNodes, normProj, projDict = self.geom.project_on_surface(blendedNodes)

            # Run quality check on regenerated mesh
            minQuality_regen, maxQuality_regen, quality_regen = check_collar_quality(collarMeshObj)

            # Assign blended nodes to the collar
            collarMeshObj.set_points(projBlendedNodes)

            # Run quality check on blended mesh (note that we already replaced the nodes)
            minQuality, maxQuality, quality = check_collar_quality(collarMeshObj)

            # Run quality checks for the blended mesh
            badQuality = False
            if np.min(quality/quality_regen) < 0.5:
                badQuality = True

            if badQuality:

                # Go back to regenerated mesh if quality is too bad
                collarMeshObj.set_points(regenNodes)

                # Update the reference mesh
                self.project_collar(collarMeshObj)

# AUXILIARY FUNCTIONS

def get_tria_proj(X1, X2, X3, u, v):

    '''
    This function finds the point in space defined by the parametric
    coordinates u and v on the triangle with nodes X1, X2, and X3.
    '''

    # Get nodal weights
    w1 = 1.0 - u - v
    w2 = u
    w3 = v

    # Compute physical node
    Xf = w1*X1 + w2*X2 + w3*X3

    # Return physical node
    return Xf

def get_quad_proj(X1, X2, X3, X4, u, v):

    '''
    This function finds the point in space defined by the parametric
    coordinates u and v on the quad with nodes X1, X2, X3, and X4.
    '''

    # Get nodal weights
    w1 = (1.0-u)*(1.0-v)
    w2 = u*(1.0-v)
    w3 = u*v
    w4 = (1.0-u)*v

    # Compute physical node
    Xf = w1*X1 + w2*X2 + w3*X3 + w4*X4

    # Return physical node
    return Xf

def build_weight_array(collarMesh, bDec=-5.0, bExp=2.0):

    '''
    This function will create the weight array used to blend
    the warped mesh and the regenerated mesh.
    We need to call this once we set a new regenerated mesh
    as reference. We do this to avoid weight dependency to the
    nodal coordinates, facilitating the differentiation process.
    Here we assume that the intersection curve is at jlow, and the overset
    boundary is at jhigh.
    '''

    # Get name of the collar block (assuming collar mesh always have one block)
    blockName = collarMesh.coor.keys()[0]

    # Get collar mesh coordinates in block form
    coor = collarMesh.coor[blockName]

    # Get mesh size
    ni,nj = coor[:,:,0].shape

    # Initialize weight matrix
    # weights near 1.0 preserves the regenerated mesh
    # weights near 0.0 preserves the warped mesh
    weights = np.ones((ni,nj),dtype=float)

    # Now we need to compute the variations in distance for each node
    # along the j direction
    dcoor = coor[:,1:,:] - coor[:,:-1,:]
    dist = np.sqrt(np.sum(dcoor**2, axis=2))

    # Accumulate distances along the j direction
    for jj in range(1,dist.shape[1]):
        dist[:,jj] = dist[:,jj-1] + dist[:,jj]

    # Now we normalize the distances by the maximum distance in the i direction
    dist = dist/dist[:,-1][:,None]

    # Now we use the blending function to compute weights.
    # We do not compute weights for j=0 because it is already 1.0
    weights[:,1:] = np.exp(bDec*dist**bExp)

    # Flatten the weights matrix
    # We need the order='F' so that we add all j=1 points first,
    # and then proceed to j=2, ...
    # Note that this does not change the actual memory layout
    # but just the reading order.
    weights_flat = weights.reshape((-1,1),order='F')

    # Return the weight matrix
    return weights_flat

def check_collar_quality(collarMesh):

    '''
    This checks the quality of the collar mesh surface cells.
    Only the first block is checked since it is the only expected one.
    '''

    # Get name of the collar block (assuming collar mesh always have one block)
    blockName = collarMesh.coor.keys()[0]

    # Get collar mesh coordinates in block form
    coor = collarMesh.coor[blockName]

    # Run quality check
    minQuality, maxQuality, quality = mesh_quality(coor)

    # Return quality metrics
    return minQuality, maxQuality, quality

def mesh_quality(X):

    '''
    X[ni,nj,3] : mesh nodes
    '''

    # Get mesh size
    ni = X.shape[0]
    nj = X.shape[1]

    # We need to compute Jacobians on every node.
    # First let's compute derivatives.
    # BUUT even before that... we allocate arrays for the derivatives
    dXdu = np.zeros(X.shape)
    dXdv = np.zeros(X.shape)
    normal = np.zeros(X.shape)

    # U derivatives
    # Compute derivatives for the center nodes
    dXdu[1:-1,:,:] = (X[2:,:,:] - X[:-2,:,:])/2.0/ni
    # Compute derivatives for the boundaries
    dXdu[0,:,:] = (X[1,:,:] - X[0,:,:])/ni
    dXdu[-1,:,:] = (X[-1,:,:] - X[-2,:,:])/ni

    # V derivatives
    # Compute derivatives for the center nodes
    dXdv[:,1:-1,:] = (X[:,2:,:] - X[:,:-2,:])/2.0/nj
    # Compute derivatives for the boundaries
    dXdv[:,0,:] = (X[:,1,:] - X[:,0,:])/nj
    dXdv[:,-1,:] = (X[:,-1,:] - X[:,-2,:])/nj

    # Normals
    for ii in range(ni):
        for jj in range(nj):

            # The derivatives roughly represent the tangent directions,
            # so we can use them to compute normals
            currNormal = np.cross(dXdu[ii,jj,:], dXdv[ii,jj,:])
            normal[ii,jj,:] = currNormal/np.linalg.norm(currNormal)

    # Now we can compute the determinants of the Jacobians for every node
    # First we allocate the arrays of determinants
    dets = np.zeros((X.shape[0],X.shape[1]))

    # Loop over every node
    for ii in range(ni):
        for jj in range(nj):

            # Assemble Jacobian matrix
            Jac = np.array([dXdu[ii,jj,:], dXdv[ii,jj,:], normal[ii,jj,:]])

            # Compute and store determinant
            dets[ii,jj] = np.linalg.det(Jac)

    # Compute quality of every cell.
    # The quality is the ratio between the maximum and minimum determinants
    # of all nodes of a cell

    # Allocate array of qualities
    quality = np.zeros((X.shape[0]-1,X.shape[1]-1))

    # Loop over every CELL
    for ii in range(ni-1):
        for jj in range(nj-1):

            # Get minimum determinant of the surrounding nodes
            minDet = min(dets[ii,  jj],
                         dets[ii+1,jj],
                         dets[ii,  jj+1],
                         dets[ii+1,jj+1])

            # Get maximum determinant of the surrounding nodes
            maxDet = max(dets[ii,  jj],
                         dets[ii+1,jj],
                         dets[ii,  jj+1],
                         dets[ii+1,jj+1])

            # Compute quality
            quality[ii,jj] = minDet/maxDet

    # Compute overall quality
    minQuality = np.min(quality)
    maxQuality = np.max(quality)

    # RETURNS
    return minQuality, maxQuality, quality
