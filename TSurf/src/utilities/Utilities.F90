module Utilities

use precision
implicit none

contains

!============================================================

subroutine condenseBarNodes_main(nNodes, nElem, distTol, &
                                 coor, barsConn, nUniqueNodes)

  ! This subroutine receives a list of bar FE which may have repeated points and
  ! Then condenses (merges) points that are closer than a given tolerance.
  ! This helps us getting a continuous line from the intersection algorithm results.
  ! Ney Secco 2016-08
  !
  ! INPUTS
  !
  ! distTol: real -> Distance tolerance used to determine if two nodes are effectively the same.
  !
  ! OUTPUTS
  !
  ! nUniqueNodes: integer -> Number of unique nodes.
  !
  ! INPUTS/OUTPUTS
  !
  ! coor: real(3,nNodes) -> Nodal coordinates. May have repeated coordinates at the beginning.
  !                         At the end, it will have zeros if the total number of nodes is reduced.
  !                         Then you should use nUniqueNodes to crop the matrix accordingly.
  !
  ! barsConn: int(2,nElements) -> Bar element connectivity array. Its values will be updated
  !           to use newCoor indices.


  implicit none

  ! INPUTS
  integer(kind=intType), intent(in) :: nNodes, nElem
  real(kind=realType), intent(in) :: distTol

  ! INPUTS/OUTPUTS
  real(kind=realType), dimension(3,nNodes), intent(inout) :: coor
  integer(kind=intType), dimension(2,nElem), intent(inout) :: barsConn

  ! OUTPUTS
  integer(kind=intType), intent(out) :: nUniqueNodes

  ! WORKING
  integer(kind=intType) :: nCopies
  integer(kind=intType) :: currNodeID, prevNodeID, link, elemID
  real(kind=realType), dimension(3) :: currCoor, prevCoor
  real(kind=realType) :: dist
  integer(kind=intType), dimension(size(coor,2)) :: linkOld2New

  ! EXECUTION

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
        dist = norm2(currCoor - prevCoor)

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

  ! Initialize number of nodes copied so far
  nCopies = 0

  ! We loop once again over the nodes so we can copy the unique values
  do currNodeID = 1,nNodes

     ! Get index of the current node in the new coordinate array
     link = linkOld2New(currNodeID)

     ! Check if the new link is already used
     if (link .gt. nCopies) then

        ! Get coordinates of current node
        currCoor = coor(:,currNodeID)

        ! Increment number of copies done so far
        nCopies = nCopies + 1

        ! Copy coordinates
        coor(:,nCopies) = currCoor

     end if

  end do

  ! Now zero out the unused point in coor
  coor(:,nCopies+1:nNodes) = 0.0

  ! Now the last step is updating the bars connectivity.
  ! Loop over the elements
  do elemID = 1,nElem

     ! Update connectivities
     barsConn(1,elemID) = linkOld2New(barsConn(1,elemID))
     barsConn(2,elemID) = linkOld2New(barsConn(2,elemID))

  end do

end subroutine condenseBarNodes_main

!============================================================

end module Utilities
