module Intersection

use precision
implicit none

contains


subroutine computeBoundingBox(coor, BBox)

  real(kind=realType), dimension(:,:), intent(in) :: coor

  real(kind=realType), dimension(3,2), intent(out) :: BBox

  BBox(:, 1) = minval(coor, 2)
  BBox(:, 2) = maxval(coor, 2)

  print *,BBox

end subroutine computeBoundingBox

subroutine computeIntervals(coor, BBox)

  real(kind=realType), dimension(:,:), intent(in) :: coor

  real(kind=realType), dimension(3,2), intent(out) :: BBox

  BBox(:, 1) = minval(coor, 2)
  BBox(:, 2) = maxval(coor, 2)

  print *,BBox

end subroutine computeIntervals





end module Intersection
