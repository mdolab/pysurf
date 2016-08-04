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

subroutine triTriIntersect(V0, V1, V2, U0, U1, U2, intersect)

  real(kind=realType), dimension(3), intent(in) :: V0, V1, V2, U0, U1, U2
  integer(kind=intType), intent(out) :: intersect

  real(kind=realType), dimension(3) :: E1, E2, N1
  real(kind=realType) :: d1, du0, du1, du2, du0du1, du0du2


  E1 = V1 - V0
  E2 = V2 - V0
  cross_product(E1, E2, N1)
  d1 = -dot_product(N1, V0)

  du0 = dot_product(N1, U0) + d1
  du1 = dot_product(N1, U1) + d1
  du2 = dot_product(N1, U2) + d1

  if (USE_EPSILON_TEST) then
    if (abs(du0)<EPSILON) du0 = 0.0
    if (abs(du1)<EPSILON) du1 = 0.0
    if (abs(du2)<EPSILON) du2 = 0.0
  end if

  du0du1 = du0 * du1
  du0du2 = du0 * du2

  if (du0du1 > 0.0 .and. du0du2 > 0.0) then
    intersect = 0
    return
  end if


end subroutine triTriIntersect

subroutine cross_product(A, B, C)

  implicit none

  real(kind=realType), intent(in) :: A(3), B(3)
  real(kind=realType), intent(out) :: C(3)

  C(1) = A(2) * B(3) - A(3) * B(2)
  C(2) = A(3) * B(1) - A(1) * B(3)
  C(3) = A(1) * B(2) - A(2) * B(1)

end subroutine cross_product





end module Intersection
