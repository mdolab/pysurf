module Intersection

use Utilities
use precision
implicit none

contains


subroutine computeBBox(coor, BBox)

  ! This subroutine computes the bounding box of all points given in coor.
  ! Ney Secco - 2016-08
  !
  ! INPUTS:
  !
  ! coor: real(3, nNodes) -> coordinates (x,y,z) of all nodes used in the
  !       triangulated surface definition.
  !
  ! OUTPUTS:
  !
  ! BBox: real(3, 2) -> coordinates defining the bounding box
  !                     [xmin, xmax]
  !                     [ymin, ymax]
  !                     [zmin, zmax]


  implicit none

  ! INPUTS
  real(kind=realType), dimension(:,:), intent(in) :: coor

  ! OUTPUTS
  real(kind=realType), dimension(3,2), intent(out) :: BBox

  ! EXECUTION

  ! Get bounding values
  BBox(:, 1) = minval(coor, 2)
  BBox(:, 2) = maxval(coor, 2)

end subroutine computeBBox

!============================================================

subroutine computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

  ! This subroutine gives a rectangle (BBoxAB) that correponds to the intersection
  ! of two bounding boxes (BBoxA, and BBoxB).
  ! If there is no intersection, we will return overlap = .false., but BBoxAB will
  ! still have meaningless value.
  ! Ney Secco - 2016-08
  !
  ! INPUTS
  !
  ! BBoxA: real(3,2) -> bounding box coordinates of component A
  !
  ! BBoxA: real(3,2) -> bounding box coordinates of component B
  !
  ! OUTPUTS
  !
  ! BBoxAB: real(3,2) -> coordinates of the A-B intersection bounding box.
  !         ATTENTION: This variable may still have meaningless values if
  !         overlap = .false.
  !
  ! overlap: logical -> logical indicating if there is an intersection
  !          between boxes A and B. If overlap = .false., disregard all
  !          values given in BBoxAB.

  implicit none

  ! INPUTS
  real(kind=realType), dimension(3,2), intent(in) :: BBoxA, BBoxB

  ! OUTPUTS
  real(kind=realType), dimension(3,2), intent(out) :: BBoxAB
  logical, intent(out) :: overlap

  ! EXECUTION

  ! Check overlaps along each dimension
  ! We have an actual BBox intersection if there are overlaps in all dimensions

  ! X overlap
  call lineIntersectionInterval(BBoxA(1,1), BBoxA(1,2), BBoxB(1,1), BBoxB(1,2), &
                                BBoxAB(1,1), BBoxAB(1,2), overlap)

  if (overlap) then
     ! Y overlap
     call lineIntersectionInterval(BBoxA(2,1), BBoxA(2,2), BBoxB(2,1), BBoxB(2,2), &
                                   BBoxAB(2,1), BBoxAB(2,2), overlap)

     if (overlap) then
        ! Z overlap
        call lineIntersectionInterval(BBoxA(3,1), BBoxA(3,2), BBoxB(3,1), BBoxB(3,2), &
                                      BBoxAB(3,1), BBoxAB(3,2), overlap)

     end if

  end if

  ! Remember that overlap may be set to false in the Z overlap test

end subroutine computeBBoxIntersection

!============================================================

subroutine lineIntersectionInterval(xminA, xmaxA, xminB, xmaxB, xminAB, xmaxAB, overlap)

  ! This subroutine finds a 1D overlap interval between two lines.
  ! This step is used three times (one for each dimension) when
  ! computing bounding boxes intersection.
  ! If the bounding boxes do not overlap we return overlap = .false.
  ! Ney Secco - 2016-08

  implicit none

  ! INPUTS
  real(kind=realType), intent(in) :: xminA, xmaxA, xminB, xmaxB

  ! OUTPUTS
  real(kind=realType), intent(out) :: xminAB, xmaxAB
  logical, intent(out) :: overlap

  ! EXECUTION

  ! Get bounding interval
  xminAB = max(xminA, xminB)
  xmaxAB = min(xmaxA, xmaxB)

  ! Check if we actually have an overlap
  if (xminAB .gt. xmaxAB) then
     overlap = .false.
  else
     overlap = .true.
  end if

end subroutine lineIntersectionInterval

!============================================================

subroutine filterElements(coor, triaConn, quadsConn, BBox, &
                          innerTriaID, innerQuadsID)

  ! This subroutine finds elements that are inside (or partially inside)
  ! the given bounding box.
  ! Ney Secco - 2016-08
  !
  ! INPUTS:
  !
  ! coor: real(3,nNodes) -> Nodal coordinates of the local grid.
  !
  ! triaConn: int(3,nTria) -> Local connectivity of the triangles.
  !
  ! quadsConn: int(4,nQuads) -> Idem for the quadrilaterals.
  !
  ! BBox: real(3,2) -> Corners of the bounding box.
  !
  ! OUTPUTS:
  !
  ! innerTriaID: integer(nInnerTria) -> Indices of triangle elements that
  !              are inside the bounding box.
  !
  ! innerQuadsID: integer(nInnerQuads) -> Indices of quad elements that
  !               are inside the bounding box.

  implicit none

  ! INPUTS
  real(kind=realType), dimension(:,:), intent(in) :: coor
  integer(kind=intType), dimension(:,:), intent(in) :: triaConn, quadsConn
  real(kind=realType), dimension(3,2), intent(in) :: BBox

  ! OUTPUTS
  integer(kind=intType), dimension(:), allocatable, intent(out) :: innerTriaID, innerQuadsID

  ! WORKING
  integer(kind=intType), dimension(:), allocatable :: extInnerTriaID, extInnerQuadsID
  logical, dimension(:), allocatable :: nodeLocation
  integer(kind=intType), dimension(:, :), allocatable :: nodeOnTop
  integer(kind=intType) :: nNodes, nTria, nQuads
  integer(kind=intType) :: elemID, nodeID
  integer(kind=intType) :: numInnerTria, numInnerQuads, i
  integer(kind=intType), dimension(3) :: onTop
  real(kind=realType) :: nodeX, nodeY, nodeZ
  logical :: nodeIsInside

  ! EXECUTION

  ! Get problem size
  nNodes = size(coor,2)
  nTria = size(triaConn,2)
  nQuads = size(quadsConn,2)

  ! Allocate array that states whether a node is inside or outside of the BBox
  ! nodeLocation(nodeID) = .false. -> node is outside the BBox
  ! nodeLocation(nodeID) = .true.  -> node is inside the BBox
  allocate(nodeLocation(nNodes))
  allocate(nodeOnTop(nNodes, 3))

  ! Initialize arrays
  nodeLocation = .false.
  nodeOnTop = 0

  ! Now check each node
  nodeLoop: do nodeID = 1,nNodes

    ! Get current node coordinates
    nodeX = coor(1,nodeID)
    nodeY = coor(2,nodeID)
    nodeZ = coor(3,nodeID)

    ! Check each dimension
    if ((nodeX .ge. BBox(1,1)) .and. (nodeX .le. BBox(1,2))) then
      if ((nodeY .ge. BBox(2,1)) .and. (nodeY .le. BBox(2,2))) then
         if ((nodeZ .ge. BBox(3,1)) .and. (nodeZ .le. BBox(3,2))) then
            nodeLocation(nodeID) = .true.
         end if
      end if
    end if

    ! Check to see if all of the nodes are on the same side of the bounding box
    ! Check X dimension
    if (nodeX .ge. BBox(1,2)) then
      nodeOnTop(nodeID, 1) = 1
    end if

    ! Check Y dimension
    if (nodeY .ge. BBox(2,2)) then
      nodeOnTop(nodeID, 2) = 1
    end if

    ! Check Z dimension
    if (nodeZ .ge. BBox(3,2)) then
      nodeOnTop(nodeID, 3) = 1
    end if

  end do nodeLoop

  ! Allocate extended arrays to store indices of the interior elements. In the worst
  ! case, all elements would be inside, so we would need to store nTria and nQuads
  ! indices in total.
  allocate(extInnerTriaID(nTria), extInnerQuadsID(nQuads))

  ! Now we check all elements.
  ! If one node is inside the bounding box, then we flag the element as "inside"

  ! Initialize interior element counters
  numInnerTria = 0
  numInnerQuads = 0

  ! First we check triangles
  triaLoop: do elemID = 1,nTria

    onTop(:) = 0

    ! Loop over all the nodes of the element
    do nodeID = 1,3

      ! Get location flag of the current node
      nodeIsInside = nodeLocation(triaConn(nodeID, elemID))

      ! Flag the element if any node is inside
      if (nodeIsInside) then

        ! We do not need to check other nodes, so
        ! we jump out of the do loop
        exit

      end if

      ! Get a running sum of which nodes are on one side of the bounding box
      onTop = onTop + nodeOnTop(triaConn(nodeID, elemID), :)

    end do

    ! Check to see if all nodes are on the same side of the bounding box
    ! If any node is not on the same side as all the others, flag the element
    ! to check for intersections
    if (nodeIsInside .neqv. .true.) then !DEBUG
      do i=1,3
        if ((onTop(i) .ne. 0) .and. (onTop(i) .ne. 3)) then
          nodeIsInside = .true.
          exit
        end if
      end do
    end if

    if (nodeIsInside) then
      ! Increment number of inside elements
      numInnerTria = numInnerTria + 1

      ! Store element index in the extended array
      extInnerTriaID(numInnerTria) = elemID
    end if

  end do triaLoop

  ! Now we check quads
  quadsLoop: do elemID = 1,nQuads

    onTop = 0

    ! Loop over all the nodes of the element
    do nodeID = 1,4

      ! Get location flag of the current node
      nodeIsInside = nodeLocation(quadsConn(nodeID, elemID))

      ! Flag the element if any node is inside
      if (nodeIsInside) then

        ! We do not need to check other nodes, so
        ! we jump out of the do loop
        exit

      end if

      ! Get a running sum of which nodes are on one side of the bounding box
      onTop = onTop + nodeOnTop(quadsConn(nodeID, elemID), :)
    end do

    ! Check to see if all nodes are on the same side of the bounding box
    ! If any node is not on the same side as all the others, flag the element
    ! to check for intersections
    if (nodeIsInside .neqv. .true.) then !DEBUG
      do i=1,3
        if ((onTop(i) .ne. 0) .and. (onTop(i) .ne. 4)) then
          nodeIsInside = .true.
          exit
        end if
      end do
    end if

    if (nodeIsInside) then

     ! Increment number of inside elements
     numInnerQuads = numInnerQuads + 1

     ! Store element index in the extended array
     extInnerQuadsID(numInnerQuads) = elemID
    end if

  end do quadsLoop

  ! Now that we know exactly how many elements are inside, we can
  ! allocate outputs arrays with the proper size
  allocate(innerTriaID(numInnerTria), innerQuadsID(numInnerQuads))

  ! Transfer values from the extended arrays to the output arrays
  innerTriaID(:) = extInnerTriaID(1:numInnerTria)
  innerQuadsID(:) = extInnerQuadsID(1:numInnerQuads)

  ! We can finally deallocate the extended arrays
  deallocate(extInnerTriaID, extInnerQuadsID)

end subroutine filterElements

!============================================================

subroutine condenseBarFEs(distTol, coor, barsConn, newCoor)

  ! This subroutine receives a list of bar FE which may have repeated points and
  ! Then condenses (merges) points that are closer than a given tolerance.
  ! This helps us getting a continuous line from the intersection algorithm results.
  ! Ney Secco 2016-08
  !
  ! INPUTS
  !
  ! distTol: real -> Distance tolerance used to determine if two nodes are effectively the same.
  !
  ! coor: real(3,nNodes) -> Nodal coordinates. May have repeated coordinates.
  !
  ! OUTPUTS
  !
  ! newCoor: real(3,nUniqueNodes) -> Nodal coordinates with unique nodes only.
  !
  ! INPUTS/OUTPUTS
  !
  ! barsConn: int(2,nElements) -> Bar element connectivity array. Its values will be updated
  !           to use newCoor indices.


  implicit none

  ! INPUTS
  real(kind=realType), intent(in) :: distTol
  real(kind=realType), dimension(:,:), intent(in) :: coor

  ! OUTPUTS
  real(kind=realType), dimension(:,:), allocatable, intent(out) :: newCoor

  ! INPUTS/OUTPUTS
  integer(kind=intType), dimension(:,:), intent(inout) :: barsConn

  ! WORKING
  integer(kind=intType) :: nNodes, nElem, nUniqueNodes, nCopies
  integer(kind=intType) :: currNodeID, prevNodeID, link, elemID
  real(kind=realType), dimension(3) :: currCoor, prevCoor
  real(kind=realType) :: dist
  integer(kind=intType), dimension(size(coor,2)) :: linkOld2New

  ! EXECUTION

  ! Get problem size
  nNodes = size(coor,2)
  nElem = size(barsConn,2)

  ! Initialize number of unique nodes found so far.
  ! The first node is a unique one =P!
  nUniqueNodes = 1

  ! As the first node is unique, its old-to-new link should point to itself
  linkOld2New(1) = 1

  ! Loop over the nodes to find the unique ones
  do currNodeID = 2,nNodes

     ! Get coordinates of current node
     currCoor = coor(:,currNodeID)

     ! Now loop over the previous nodes to find if it is repeated
     do prevNodeID = 1,currNodeID-1

        ! Get coordinates of the previous node
        prevCoor = coor(:,prevNodeID)

        ! Compute distance between nodes
        call norm(currCoor-prevCoor, dist)

        ! Check if the distance is below the merging tolerance
        if (dist .le. distTol) then

           ! Update link array.
           ! The array linkOld2New will contain newCoor indices that correspond to each
           ! coor index.
           ! So the current node should use the same link as the previous node, as they will
           ! point to the same index of the new coordinate array.
           linkOld2New(currNodeID) = linkOld2New(prevNodeID)

           ! We can jump out of the prevNode do loop and check the next node.
           exit

        end if

     end do

     ! Check if we did not find any copy. In this case, we need to initialize a new
     ! unique node
     if (prevNodeID .eq. currNodeID) then

        ! Increase the number of unique nodes found so far
        nUniqueNodes = nUniqueNodes + 1

        ! Create new link
        linkOld2New(currNodeID) = nUniqueNodes

     end if

  end do

  ! Now that we know the number of unique nodes, we can allocate memory for the new
  ! coordinate array.
  allocate(newCoor(3,nUniqueNodes))
  newCoor = 0.0

  ! Initialize number of nodes copied so far
  nCopies = 0

  ! We loop once again over the nodes so we can copy the unique values
  do currNodeID = 1,nNodes

     ! Get coordinates of current node
     currCoor = coor(:,currNodeID)

     ! Get index of the current node in the new coordinate array
     link = linkOld2New(currNodeID)

     ! Check if the new link is already used
     if (link .gt. nCopies) then

        ! Increment number of copies done so far
        nCopies = nCopies + 1

        ! Copy coordinates
        newCoor(:,nCopies) = currCoor

     end if

  end do

  ! Now the last step is updating the bars connectivity.
  ! Loop over the elements
  do elemID = 1,nElem

     ! Update connectivities
     barsConn(1,elemID) = linkOld2New(barsConn(1,elemID))
     barsConn(2,elemID) = linkOld2New(barsConn(2,elemID))

  end do

end subroutine condenseBarFEs

!============================================================

subroutine getAllTrias(triaConn, quadsConn, innerTriaID, innerQuadsID, allTriaConn)

  implicit none

  ! INPUTS
  integer(kind=intType), dimension(:,:), intent(in) :: triaConn, quadsConn
  integer(kind=intType), dimension(:), intent(in) :: innerTriaID, innerQuadsID

  ! OUTPUTS
  integer(kind=intType), dimension(:,:), allocatable, intent(out) :: allTriaConn

  ! WORKING
  integer(kind=intType) :: nInnerTria, nInnerQuads, ii, node1, node2, node3, node4

  ! EXECUTION

  ! Get number of interior elements
  nInnerTria = size(innerTriaID)
  nInnerQuads = size(innerQuadsID)

  ! Allocate connectivity arrays for triangles split from quads
  allocate(allTriaConn(3, nInnerTria + nInnerQuads*2))

  ! Copy connectivities of the interior triangles
  do ii = 1,nInnerTria

     allTriaConn(:,ii) = triaConn(:, innerTriaID(ii))

  end do

  ! Loop over every interior element
  do ii = 1,nInnerQuads

     ! Get nodes of the current quad
     node1 = quadsConn(1, innerQuadsID(ii))
     node2 = quadsConn(2, innerQuadsID(ii))
     node3 = quadsConn(3, innerQuadsID(ii))
     node4 = quadsConn(4, innerQuadsID(ii))

     ! Create two triangle elements
     allTriaConn(:, nInnerTria + 2*ii-1) = [node1, node2, node3]
     allTriaConn(:, nInnerTria + 2*ii)   = [node3, node4, node1]

  end do

end subroutine getAllTrias

!============================================================

subroutine triTriIntersect(V0, V1, V2, U0, U1, U2, intersect, vecStart, vecEnd)

  ! This subroutine computes the line vector of intersection between two triangles
  ! John Jasa - 2016-08
  !
  ! Adapted from Moller's 1997 paper "A Fast Triangle-Triangle Intersection Test"
  !
  ! INPUTS:
  !
  ! V0: real(3) -> coordinates (x,y,z) of the first node of triangle V
  ! V0: real(3) -> coordinates (x,y,z) of the second node of triangle V
  ! V0: real(3) -> coordinates (x,y,z) of the third node of triangle V
  !
  ! U0: real(3) -> coordinates (x,y,z) of the first node of triangle U
  ! U0: real(3) -> coordinates (x,y,z) of the second node of triangle U
  ! U0: real(3) -> coordinates (x,y,z) of the third node of triangle U
  !
  ! OUTPUTS:
  !
  ! intersect: integer -> 1 if the triangles intersect, 0 if they do not
  ! vecStart: real(3) -> coordinates (x,y,z) of the start point of the intersection line
  ! vecEnd: real(3) -> coordinates (x,y,z) of the end point of the intersection line

  ! INPUTS
  real(kind=realType), dimension(3), intent(in) :: V0, V1, V2, U0, U1, U2

  ! OUTPUTS
  integer(kind=intType), intent(out) :: intersect
  real(kind=realType), dimension(3), intent(out) :: vecStart, vecEnd

  ! WORKING
  real(kind=realType), dimension(3) :: E1, E2, N1, N2, Dir
  real(kind=realType) :: d1, du0, du1, du2, du0du1, du0du2, epsilon
  real(kind=realType) :: d2, dv0, dv1, dv2, dv0dv1, dv0dv2, maxD, bb, cc
  real(kind=realType) :: up0, up1, up2, vp0, vp1, vp2
  real(kind=realType) :: isect1(2), isect2(2)
  real(kind=realType), dimension(3) :: isectpointA1, isectpointA2, isectpointB1, isectpointB2
  integer(kind=intType) :: index
  integer(kind=intType) :: coplanar, smallest1, smallest2

  ! Initialize intersect and coplanar values so the program does not stop prematurely
  intersect = 2
  coplanar = 0

  ! Compute plane of triangle (V0, V1, V2)
  E1 = V1 - V0
  E2 = V2 - V0
  call cross_product(E1, E2, N1)
  call dot(N1, V0, d1)
  d1 = -d1

  ! Get distances from U points to plane defined by V points
  call dot(N1, U0, du0)
  du0 = du0 + d1
  call dot(N1, U1, du1)
  du1 = du1 + d1
  call dot(N1, U2, du2)
  du2 = du2 + d1

  ! Compute the signed distance product to see which side of the plane each point is on
  du0du1 = du0 * du1
  du0du2 = du0 * du2

  ! If all the points of one triangle are on the same side of the other triangle,
  ! there is no intersection
  if (du0du1 > 0.0 .and. du0du2 > 0.0) then
    intersect = 0
    return
  end if

  ! Compute plane of triangle (U0, U1, U2)
  E1 = U1 - U0
  E2 = U2 - U0
  call cross_product(E1, E2, N2)
  call dot(N2, U0, d2)
  d2 = -d2

  ! Get distances from V points to plane defined by U points
  call dot(N2, V0, dv0)
  dv0 = dv0 + d2
  call dot(N2, V1, dv1)
  dv1 = dv1 + d2
  call dot(N2, V2, dv2)
  dv2 = dv2 + d2

  ! Compute the signed distance product to see which side of the plane each point is on
  dv0dv1 = dv0 * dv1
  dv0dv2 = dv0 * dv2

  ! If all the points of one triangle are on the same side of the other triangle,
  ! there is no intersection
  if (dv0dv1 > 0.0 .and. dv0dv2 > 0.0) then
    intersect = 0
    return
  end if

  ! Compute the direction of the intersection line
  call cross_product(N1, N2, Dir)

  ! Compute and index the largest component of Dir
  maxD = abs(Dir(1))
  index = 1
  bb = abs(Dir(2))
  cc = abs(Dir(3))
  if (bb > maxD) then
    maxD = bb
    index = 2
  end if
  if (cc > maxD) then
    index = 3
  end if

  ! Simplified projection onto L
  vp0 = V0(index)
  vp1 = V1(index)
  vp2 = V2(index)

  up0 = U0(index)
  up1 = U1(index)
  up2 = U2(index)

  ! isect1 and isect2 are the projected intersections of triangles 1 and 2,
  ! which are defined by the V and U points respectively.
  ! These isect values show where each triangle edge intersects the
  ! two triangle's intersection line

  ! Compute the intersection interval for the V points
  call compute_intervals_isectline(V0, V1, V2, vp0, vp1, vp2, dv0, dv1, dv2, &
   dv0dv1, dv0dv2, isect1(1), isect1(2), isectpointA1, isectpointA2, coplanar, intersect)

  if (intersect .eq. 0) return

  ! Compute the intersection interval for the U points
  call compute_intervals_isectline(U0, U1, U2, up0, up1, up2, du0, du1, du2, &
    du0du1, du0du2, isect2(1), isect2(2), isectpointB1, isectpointB2, coplanar, intersect)

  if (intersect .eq. 0) return

  ! Sort the projected intersections so that the first index contains the
  ! smallest value. Also index which case has the smallest max so we
  ! can compute the actual intersection line
  call sort(isect1(1), isect1(2), smallest1)
  call sort(isect2(1), isect2(2), smallest2)

  ! If there is no interval where the isects overlap, there there is no intersection
  if ((isect1(2) .lt. isect2(1)) .or. (isect2(2) .lt. isect1(1))) then
    intersect = 0

  ! If the intersection is only a point, we do not treat it as an intersection
  else if (max(isect1(1), isect2(1)) .eq. min(isect1(2), isect2(2))) then
    intersect = 0

  ! Only continue if the triangles are not coplanar.
  ! Choose the triangle edges that intersect the intersection line in the
  ! most restrictive interval.
  else if (coplanar .ne. 1) then

    if (isect2(1) .lt. isect1(1)) then

      if (smallest1 .eq. 0) then
        vecStart = isectpointA1
      else
        vecStart = isectpointA2
      end if

      if (isect2(2) .lt. isect1(2)) then

        if (smallest2 .eq. 0) then
          vecEnd = isectpointB2
        else
          vecEnd = isectpointB1
        end if

      else

        if (smallest1 == 0) then
          vecEnd = isectpointA2
        else
          vecEnd = isectpointA1
        end if

      end if

    else

      if (smallest2 .eq. 0) then
        vecStart = isectpointB1
      else
        vecStart = isectpointB2
      end if

      if (isect2(2) .gt. isect1(2)) then

        if (smallest1 .eq. 0) then
          vecEnd = isectpointA2
        else
          vecEnd = isectpointA1
        end if

      else

        if (smallest2 .eq. 0) then
          vecEnd = isectpointB2
        else
          vecEnd = isectpointB1
        end if

      end if

    end if

    intersect = 1

  end if

end subroutine triTriIntersect

subroutine compute_intervals_isectline(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, &
  			                           D0D1, D0D2, isect0, isect1, isectpoint0, isectpoint1, coplanar, intersect)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(in) :: VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2
  real(kind=realType), dimension(3), intent(in) :: VERT0, VERT1, VERT2
  real(kind=realType), intent(out) :: isectpoint0(3), isectpoint1(3), isect0, isect1
  integer(kind=intType), intent(out) :: coplanar, intersect

  ! Test if d0 and d1 are on the same side
  if (D0D1 .gt. 0.0) then
    call intersect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, &
      isect0, isect1, isectpoint0, isectpoint1)

  ! Test if d0 and d2 are on the same side
  else if (D0D2 .gt. 0.0) then
    call intersect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, &
      isect0, isect1, isectpoint0, isectpoint1)

  ! Test if d1 and d2 are on the same side or if d0 is not equal to 0
  else if ((D1*D2 .gt. 0.0) .or. (D0 .ne. 0.0)) then
    call intersect2(VERT0, VERT1, VERT2, VV0, VV1, VV2, D0, D1, D2, &
      isect0, isect1, isectpoint0, isectpoint1)

  ! Called only if d0 is 0
  else if (D1 .ne. 0.0) then
    call intersect2(VERT1, VERT0, VERT2, VV1, VV0, VV2, D1, D0, D2, &
      isect0, isect1, isectpoint0, isectpoint1)

  ! Called only if d1 and d0 is 0
  else if (D2 .ne. 0.0) then
    call intersect2(VERT2, VERT0, VERT1, VV2, VV0, VV1, D2, D0, D1, &
      isect0, isect1, isectpoint0, isectpoint1)

  else
    ! For now, we are not interested in coplanar triangle intersections
    !intersect = coplanarTriTri(N1, V0, V1, V2, U0, U1, U2)
    intersect = 0
    coplanar = 1

  end if

end subroutine compute_intervals_isectline

subroutine intersect2(VTX0, VTX1, VTX2, VV0, VV1, VV2, &
	     D0, D1, D2, isect0, isect1, isectpoint0, isectpoint1)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(in) :: VV0, VV1, VV2, D0, D1, D2
  real(kind=realType), dimension(3), intent(in) :: VTX0, VTX1, VTX2
  real(kind=realType), intent(out) :: isect0, isect1, isectpoint0(3), isectpoint1(3)

  real(kind=realType) :: tmp
  real(kind=realType), dimension(3) :: diff

  tmp = D0 / (D0 - D1)
  isect0 = VV0 + (VV1 - VV0) * tmp
  diff = VTX1 - VTX0
  diff = diff * tmp
  isectpoint0 = diff + VTX0

  tmp = D0 / (D0 - D2)
  isect1 = VV0 + (VV2 - VV0) * tmp
  diff = VTX2 - VTX0
  diff = diff * tmp
  isectpoint1 = diff + VTX0

end subroutine intersect2

function coplanarTriTri(N, V0, V1, V2, U0, U1, U2) result(intersect)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), dimension(3), intent(in) :: N, V0, V1, V2, U0, U1, U2
  integer(kind=intType) :: intersect
  real(kind=realType), dimension(3) :: A
  integer(kind=intType) :: i0, i1

  ! First project onto an axis-aligned plane, that maximizes the area
  ! of the triangles, compute indices: i0,i1.

  A(1) = abs(N(1))
  A(2) = abs(N(2))
  A(3) = abs(N(3))

  if (A(1) .gt. A(2)) then
    if (A(1) .gt. A(3)) then
      i0 = 2
      i1 = 3
    else
      i0 = 1
      i1 = 2
    end if
  else
    if (A(3) .gt. A(2)) then
      i0 = 1
      i1 = 2
    else
      i0 = 1
      i1 = 3
    end if
  end if

  ! Test all edges of triangle 1 against the edges of triangle 2
  intersect = edgeAgainstTriEdges(V0, V1, U0, U1, U2, i0, i1)
  if (intersect .eq. 1) return
  intersect = edgeAgainstTriEdges(V1, V2, U0, U1, U2, i0, i1)
  if (intersect .eq. 1) return
  intersect = edgeAgainstTriEdges(V2, V0, U0, U1, U2, i0, i1)
  if (intersect .eq. 1) return

  ! Finally, test if tri1 is totally contained in tri2 or vice versa
  intersect = pointInTri(V0, U0, U1, U2, i0, i1)
  if (intersect .eq. 1) return
  intersect = pointInTri(U0, V0, V1, V2, i0, i1)
  if (intersect .eq. 1) return

end function coplanarTriTri

function pointInTri(V0, U0, U1, U2, i0, i1) result(intersect)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), dimension(3), intent(in) :: V0, U0, U1, U2
  integer(kind=intType) :: intersect
  real(kind=realType) :: a, b, c, d0, d1, d2
  integer(kind=intType) :: i0, i1

  intersect = 0
  a = U1(i1) - U0(i1)
  b = -(U1(i0) - U0(i0))
  c = -a * U0(i0) - b * U0(i1)
  d0 = a * V0(i0) + b * V0(i1) + c

  a = U2(i1) - U1(i1)
  b = -(U2(i0) - U1(i0))
  c = -a * U1(i0) - b * U1(i1)
  d1 = a * V0(i0) + b * V0(i1) + c

  a = U0(i1) - U2(i1)
  b = -(U0(i0) - U2(i0))
  c = -a * U2(i0) - b * U2(i1)
  d2 = a * V0(i0) + b * V0(i1) + c

  if (d0 * d1 .gt. 0.) then
    if (d0 * d2 .gt. 0.) intersect = 1
    return
  end if

end function pointInTri


function edgeAgainstTriEdges(V0, V1, U0, U1, U2, i0, i1) result(intersect)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), dimension(3), intent(in) :: V0, V1, U0, U1, U2
  integer(kind=intType), intent(in) :: i0, i1
  integer(kind=intType) :: intersect

  real(kind=realType) :: Ax, Ay

  Ax = V1(i0) - V0(i0)
  Ay = V1(i1) - V0(i1)

  intersect = EDGE_EDGE_TEST(V0, U0, U1, Ax, Ay, i0, i1)
  if (intersect .eq. 1) return
  intersect = EDGE_EDGE_TEST(V0, U1, U2, Ax, Ay, i0, i1)
  if (intersect .eq. 1) return
  intersect = EDGE_EDGE_TEST(V0, U2, U0, Ax, Ay, i0, i1)
  if (intersect .eq. 1) return

end function edgeAgainstTriEdges

function EDGE_EDGE_TEST(V0, U0, U1, Ax, Ay, i0, i1) result(intersect)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), dimension(3), intent(in) :: V0, U0, U1
  real(kind=realType), intent(in) :: Ax, Ay
  integer(kind=intType), intent(in) :: i0, i1
  real(kind=realType) :: Bx, By, Cx, Cy, f, d, e

  integer(kind=intType) :: intersect

  Bx = U0(i0) - U1(i0)
  By = U0(i1) - U1(i1)
  Cx = V0(i0) - U0(i0)
  Cy = V0(i1) - U0(i1)
  f = Ay * Bx - Ax * By
  d = By * Cx - Bx * Cy

  intersect = 0
  if (((f .gt. 0) .and. (d .ge. 0) .and. (d .le. f)) .or. &
    ((f .lt. 0) .and. (d .le. 0) .and. (d .ge. f))) then
    e = Ax * Cy - Ay * Cx
    if (f>0) then
      if((e .ge. 0) .and. (e .le. f)) intersect = 1
      return
    else
      if ((e .le. 0) .and. (e .ge. f)) intersect = 1
      return
    end if
  end if

end function EDGE_EDGE_TEST

subroutine cross_product(A, B, C)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(in) :: A(3), B(3)
  real(kind=realType), intent(out) :: C(3)

  C(1) = A(2) * B(3) - A(3) * B(2)
  C(2) = A(3) * B(1) - A(1) * B(3)
  C(3) = A(1) * B(2) - A(2) * B(1)

end subroutine cross_product

subroutine sort(a, b, smallest)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(inout) :: a, b
  integer(kind=intType), intent(out) :: smallest
  real(kind=realType) :: c

  if (a .gt. b) then
    c = a
    a = b
    b = c
    smallest = 1
  else
    smallest = 0
  end if

end subroutine sort

end module Intersection
