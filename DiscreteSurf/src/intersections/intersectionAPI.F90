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
  real(kind=realType), dimension(3,2) :: BBoxA, BBoxB

  ! EXECUTION

  ! Compute bounding boxes for each component
  call computeBoundingBox(coorA, BBoxA)
  call computeBoundingBox(coorB, BBoxB)

end subroutine computeIntersection

end module
