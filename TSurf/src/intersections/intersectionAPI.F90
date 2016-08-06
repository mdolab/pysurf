module intersectionAPI

use precision
implicit none

contains

subroutine computeIntersection(nNodesA, nTriaA, nQuadsA, &
                               nNodesB, nTriaB, nQuadsB, &
                               coorA, triaConnA, quadsConnA, &
                               coorB, triaConnB, quadsConnB)

  ! This function computes the intersection curve between two
  ! triangulated surface components A and B.

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
  !f2py intent(in) nNodesA, nTriaA, nQuadsA
  !f2py intent(in) nNodesB, nTriaB, nQuadsB
  !f2py intent(in) coorA, triaConnA, quadsConnA
  !f2py intent(in) coorB, triaConnB, quadsConnB

  ! Working variables
  real(kind=realType), dimension(3,2) :: BBoxA, BBoxB, BBoxAB
  logical :: overlap
  integer(kind=intType), dimension(:), allocatable :: innerTriaID_A, innerQuadsID_A 
  integer(kind=intType), dimension(:), allocatable :: innerTriaID_B, innerQuadsID_B

  ! EXECUTION

  ! Compute bounding boxes for each component
  call computeBBox(coorA, BBoxA)
  call computeBBox(coorB, BBoxB)

  ! Compute bounding boxes intersection
  call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

  ! We can stop if there is no bounding box intersection
  if (.not. overlap) then
     print *,'Components do not overlap.'
     return
  else
     print *,'Components overlap.'
  end if

  ! Filter elements that are inside the intersected bounding box
  call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
                      innerTriaID_A, innerQuadsID_A)
  call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
                      innerTriaID_B, innerQuadsID_B)

  ! Print log
  print *,'Number of interior elements in A:'
  print *,size(innerTriaID_A) + size(innerQuadsID_A),'of',nTriaA + nQuadsA
  print *,'Number of interior elements in B:'
  print *,size(innerTriaID_B) + size(innerQuadsID_B),'of',nTriaB + nQuadsB

end subroutine computeIntersection

end module
