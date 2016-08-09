module intersectionAPI

use precision
implicit none

! Here we define output variables of unknown shape as attributes, so we
! can read them from Python

! Intersection FE data
real(kind=realType), dimension(:,:), allocatable :: coor ! Intersection nodes
integer(kind=intType), dimension(:,:), allocatable :: barsConn ! Intersection bar connectivities

contains

subroutine computeIntersection(nNodesA, nTriaA, nQuadsA, &
                               nNodesB, nTriaB, nQuadsB, &
                               coorA, triaConnA, quadsConnA, &
                               coorB, triaConnB, quadsConnB, &
                               distTol)

  ! This function computes the intersection curve between two
  ! triangulated surface components A and B.
  !
  ! INPUTS
  !
  ! distTol: real -> distance tolerance to merge nearby nodes in
  !          the intersection curve.
  !
  ! OUTPUTS
  ! This subroutine has no explicit outputs. It updates the variables coor, and barsConn,
  ! which should be called from Python as attributes of the intersectionAPI module.

  use Intersection
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nNodesA, nTriaA, nQuadsA
  integer(kind=intType), intent(in) :: nNodesB, nTriaB, nQuadsB
  real(kind=realType), dimension(3,nNodesA), intent(in) :: coorA
  integer(kind=intType), dimension(3,nTriaA), intent(in) :: triaConnA
  integer(kind=intType), dimension(4,nQuadsA), intent(in) :: quadsConnA
  real(kind=realType), dimension(3,nNodesB), intent(in) :: coorB
  integer(kind=intType), dimension(3,nTriaB), intent(in) :: triaConnB
  integer(kind=intType), dimension(4,nQuadsB), intent(in) :: quadsConnB
  real(kind=realType), intent(in) :: distTol
  !f2py intent(in) nNodesA, nTriaA, nQuadsA
  !f2py intent(in) nNodesB, nTriaB, nQuadsB
  !f2py intent(in) coorA, triaConnA, quadsConnA
  !f2py intent(in) coorB, triaConnB, quadsConnB
  !f2py intent(in) distTol

  ! Working variables
  real(kind=realType), dimension(3,2) :: BBoxA, BBoxB, BBoxAB
  logical :: overlap
  integer(kind=intType), dimension(:), allocatable :: innerTriaID_A, innerQuadsID_A
  integer(kind=intType), dimension(:), allocatable :: innerTriaID_B, innerQuadsID_B
  integer(kind=intType), dimension(:,:), allocatable :: allTriaConnA, allTriaConnB
  integer(kind=intType) :: nInnerTriaA, nInnerQuadsA, nInnerTriaB, nInnerQuadsB
  integer(kind=intType) :: ii, jj
  integer(kind=intType) :: intersect
  integer(kind=intType) :: arraySize, nAllocations, nBarsConn
  integer(kind=intType), dimension(:,:), allocatable :: extBarsConn, extTempInt
  real(kind=realType), dimension(:,:), allocatable :: extCoor, extTempReal
  real(kind=realType), dimension(3) :: node1A, node2A, node3A
  real(kind=realType), dimension(3) :: node1B, node2B, node3B
  real(kind=realType), dimension(3) :: vecStart, vecEnd

  ! EXECUTION

  ! Compute bounding boxes for each component
  call computeBBox(coorA, BBoxA)
  call computeBBox(coorB, BBoxB)

  ! Compute bounding boxes intersection
  call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

  ! We can stop if there is no bounding box intersection
  if (.not. overlap) then
     print *,'Components do not overlap.'
     allocate(barsConn(0, 0), coor(0, 0))
     return
  else
     print *,'Components overlap.'
  end if

  ! Filter elements that are inside the intersected bounding box
  call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
                      innerTriaID_A, innerQuadsID_A)
  call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
                      innerTriaID_B, innerQuadsID_B)

  ! Get number of inner elements
  nInnerTriaA = size(innerTriaID_A)
  nInnerQuadsA = size(innerQuadsID_A)
  nInnerTriaB = size(innerTriaID_B)
  nInnerQuadsB = size(innerQuadsID_B)

  ! Print log
  print *,'Number of interior elements in A:'
  print *,nInnerTriaA + nInnerQuadsA,'of',nTriaA + nQuadsA
  print *,'Number of interior elements in B:'
  print *,nInnerTriaB + nInnerQuadsB,'of',nTriaB + nQuadsB

  ! Split quads into triangles
  call getAllTrias(triaConnA, quadsConnA, innerTriaID_A, innerQuadsID_A, allTriaConnA)
  call getAllTrias(triaConnB, quadsConnB, innerTriaID_B, innerQuadsID_B, allTriaConnB)

  ! Update number of interior triangles
  nInnerTriaA = size(allTriaConnA,2)
  nInnerTriaB = size(allTriaConnB,2)

  ! Initialize extended connectivity arrays for the intersection curves.
  ! These arrays will be cropped to get the actual outputs.
  ! For now we will allocate arrays of size arraySize. We will increase them
  ! if necessary
  arraySize = 1000
  allocate(extCoor(3,arraySize), extBarsConn(2,arraySize))

  ! Initialize values
  extCoor = 0.0
  extBarsConn = 0

  ! For now, we had just one allocations
  nAllocations = 1

  ! Initialize the number of intersection connectivities known so far.
  nBarsConn = 0

  ! Compute all triangle-triangle intersections.
  ! We will call every pair between triangles in A and B
  do ii = 1,nInnerTriaA

     ! Get nodes of the element in A
     node1A = coorA(:, allTriaConnA(1,ii))
     node2A = coorA(:, allTriaConnA(2,ii))
     node3A = coorA(:, allTriaConnA(3,ii))

     ! Now loop over the triangles of the other component
     do jj = 1,nInnerTriaB

        ! Get nodes of the element in B
        node1B = coorB(:, allTriaConnB(1,jj))
        node2B = coorB(:, allTriaConnB(2,jj))
        node3B = coorB(:, allTriaConnB(3,jj))

        ! Call triangle-triangle intersection routine
        call triTriIntersect(node1A, node2A, node3A, &
                             node1B, node2B, node3B, &
                             intersect, vecStart, vecEnd)

        ! Check if the triangles actually intersect
        if (intersect .eq. 1) then

           ! Increase the number of intersections elements known so far
           nBarsConn = nBarsConn + 1

           ! Check if we already extrapolated the memory allocated so far
           if ((nBarsConn .gt. size(extBarsConn,2)) .or. (2*nBarsConn .gt. size(extCoor,2))) then

              ! We need to allocate more memory
              nAllocations = nAllocations + 1

              ! REALLOCATING extCoor

              ! Create temporary array and initialize
              allocate(extTempReal(3,nAllocations*arraySize))
              extTempReal = 0.0

              ! Tranfer data from original array
              extTempReal(:,:(nAllocations-1)*arraySize) = extCoor

              ! Now move the new allocation back to extCoor.
              ! extTemp is deallocated in this process.
              call move_alloc(extTempReal, extCoor)

              ! REALLOCATING extBarsConn

              ! Create temporary array
              allocate(extTempInt(2,nAllocations*arraySize))
              extTempInt = 0

              ! Tranfer data from original array
              extTempInt(:,:(nAllocations-1)*arraySize) = extBarsConn

              ! Now move the new allocation back to extBarsConn.
              ! extTemp is deallocated in this process.
              call move_alloc(extTempInt, extBarsConn)

           end if

           ! Assign new nodes
           extCoor(:,2*nBarsConn-1) = vecStart
           extCoor(:,2*nBarsConn)   = vecEnd

           ! Assign new connectivity
           extBarsConn(:,nBarsConn) = [2*nBarsConn-1, 2*nBarsConn ]

        end if

     end do

  end do

  ! Crop the extended connectivity array
  allocate(barsConn(2,nBarsConn))
  barsConn(:,:) = extBarsConn(:,:nBarsConn)

  ! Merge close nodes to get continuous FE data. The condensed coordinates (coor)
  ! will be returned to Python.
  call condenseBarFEs(distTol, extCoor, barsConn, coor)

end subroutine computeIntersection

subroutine testTri(V0, V1, V2, U0, U1, U2, intersect, vecStart, vecEnd)

  ! This function computes the intersection curve between two
  ! triangulated surface components A and B.

  use Intersection
  implicit none

  real(kind=realType), dimension(3), intent(in) :: V0, V1, V2, U0, U1, U2
  integer(kind=intType), intent(out) :: intersect
  real(kind=realType), dimension(3), intent(out) :: vecStart, vecEnd

  call triTriIntersect(V0, V1, V2, U0, U1, U2, intersect, vecStart, vecEnd)

end subroutine testTri

end module
