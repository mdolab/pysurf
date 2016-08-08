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

  ! EXECUTION

  ! Compute bounding boxes for each component
  ! call computeBBox(coorA, BBoxA)
  ! call computeBBox(coorB, BBoxB)

  ! Compute bounding boxes intersection
  ! call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

  ! We can stop if there is no bounding box intersection
  if (.not. overlap) then
     print *,'Components do not overlap.'
     return
  else
     print *,'Components overlap.'
  end if

  ! Filter elements that are inside the intersected bounding box
  ! call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
  !                     innerTriaID_A, innerQuadsID_A)
  ! call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
                      ! innerTriaID_B, innerQuadsID_B)

  ! Print log
  print *,'Number of interior elements in A:'
  print *,size(innerTriaID_A) + size(innerQuadsID_A),'of',nTriaA + nQuadsA
  print *,'Number of interior elements in B:'
  print *,size(innerTriaID_B) + size(innerQuadsID_B),'of',nTriaB + nQuadsB

  ! ADD CODE TO COMPUTE TRIANGLE INTERSECTIONS HERE

  ! Merge close nodes to get continuous FE data. The condensed coordinates (coor)
  ! will be returned to Python.
  ! call condenseBarFEs(distTol, intCoor, barsConn, coor)

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
