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

  real(kind=realType), dimension(3) :: E1, E2, N1, N2, Dir
  real(kind=realType) :: d1, du0, du1, du2, du0du1, du0du2, epsilon
  real(kind=realType) :: d2, dv0, dv1, dv2, dv0dv1, dv0dv2, maxD, bb, cc
  real(kind=realType) :: up0, up1, up2, vp0, vp1, vp2
  real(kind=realType) :: a, b, c, d, e, f, x0, x1, y0, y1

  integer(kind=intType) :: index

  epsilon = 1.e-5

  ! Compute plane of triangle (V0, V1, V2)
  E1 = V1 - V0
  E2 = V2 - V0
  call cross_product(E1, E2, N1)
  d1 = -dot_product(N1, V0)

  du0 = dot_product(N1, U0) + d1
  du1 = dot_product(N1, U1) + d1
  du2 = dot_product(N1, U2) + d1

  if (abs(du0)<EPSILON) du0 = 0.0
  if (abs(du1)<EPSILON) du1 = 0.0
  if (abs(du2)<EPSILON) du2 = 0.0

  du0du1 = du0 * du1
  du0du2 = du0 * du2

  if (du0du1 > 0.0 .and. du0du2 > 0.0) then
    intersect = 0
    return
  end if

  ! Compute plane of triangle (U0, U1, U2)
  E1 = U1 - U0
  E2 = U2 - U0
  call cross_product(E1, E2, N2)
  d2 = -dot_product(N2, U0)

  dv0 = dot_product(N2, V0) + d2
  dv1 = dot_product(N2, V1) + d2
  dv2 = dot_product(N2, V2) + d2

  if (abs(dv0)<EPSILON) dv0 = 0.0
  if (abs(dv1)<EPSILON) dv1 = 0.0
  if (abs(dv2)<EPSILON) dv2 = 0.0

  dv0dv1 = dv0 * dv1
  dv0dv2 = dv0 * dv2

  if (dv0dv1 > 0.0 .and. dv0dv2 > 0.0) then
    intersect = 0
    return
  end if

  ! Compute the direction of the intersection line
  call cross_product(N1, N2, Dir)

  ! Compute and index the largest component of D
  maxD = abs(Dir(1))
  index = 0
  bb = abs(Dir(2))
  cc = abs(Dir(3))
  if (bb > maxD) then
    maxD = bb
    index = 2
  end if
  if (cc > maxD) then
    maxD = cc
    index = 3
  end if

  ! Simplified projection onto L
  vp0 = V0(index)
  vp1 = V1(index)
  vp2 = V2(index)

  up0 = U0(index)
  up1 = U1(index)
  up2 = U2(index)

  call newcomputeIntervals(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, a, b, c, x0, x1, N1, V0, V1, V2, U0, U1, U2)

  call newcomputeIntervals(up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f, y0, y1, N1, V0, V1, V2, U0, U1, U2)

  ! xx = x0 * x1
  ! yy = y0 * y1
  ! xxyy = xx * yy
  !
  ! tmp = a * xxyy
  ! isect1(1) = tmp + b * x1 * yy
  ! isect1(2) = tmp + c * x0 * yy
  !
  ! tmp = d * xxyy
  ! isect2(1) = tmp + e * xx * y1
  ! isect2(2) = tmp + f * xx * y0
  !
  ! call sort(isect1(1), isect1(2))
  ! call sort(isect2(1), isect2(2))
  !
  ! if(isect1(2) .lt. isect2(1) .or. isect2(2) .lt. isect1(1)) return 0
  ! return 1

  intersect = 1


end subroutine triTriIntersect

subroutine newcomputeIntervals(VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2, A, B, C, X0, X1, N1, V0, V1, V2, U0, U1, U2)

  implicit none

  real(kind=realType), intent(in) :: VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2
  real(kind=realType), dimension(3), intent(in) :: N1, V0, V1, V2, U0, U1, U2
  real(kind=realType), intent(inout) :: A, B, C, X0, X1
  integer(kind=intType) :: ans

  if (D0D1 .gt. 0.0) then
    A = VV2
    B = (VV0 - VV2) * D2
    C = (VV1 - VV2) * D2
    X0 = D2 - D0
    X1 = D2 - D1

  else if (D0D2 .gt. 0.0) then
    A = VV1
    B = (VV0 - VV1) * D1
    C = (VV2 - VV1) * D1
    X0 = D1 - D0
    X1 = D1 - D2

  else if ((D1*D2 .gt. 0.0) .or. (D0 .ne. 0.0)) then
    A = VV0
    B = (VV1 - VV0) * D0
    C = (VV2 - VV0) * D0
    X0 = D0 - D1
    X1 = D0 - D2

  else if (D1 .ne. 0.0) then
    A = VV1
    B = (VV0 - VV1) * D1
    C = (VV2 - VV1) * D1
    X0 = D1 - D0
    X1 = D1 - D2

  else if (D2 .ne. 0.0) then
    A = VV2
    B = (VV0 - VV2) * D2
    C = (VV1 - VV2) * D2
    X0 = D2 - D0
    X1 = D2 - D1

  else
    print *,coplanarTriTri(N1, V0, V1, V2, U0, U1, U2)

  end if

end subroutine newcomputeIntervals

function coplanarTriTri(N, V0, V1, V2, U0, U1, U2) result(ans)

  implicit none

  real(kind=realType), dimension(3), intent(in) :: N, V0, V1, V2, U0, U1, U2
  integer(kind=intType) :: ans
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

    ! ! Test all edges of triangle 1 against the edges of triangle 2
    ! call edgeAgainstTriEdges(V0, V1, U0, U1, U2)
    ! call edgeAgainstTriEdges(V1, V2, U0, U1, U2)
    ! call edgeAgainstTriEdges(V2, V0, U0, U1, U2)
    !
    ! ! Finally, test if tri1 is totally contained in tri2 or vice versa
    ! call pointInTri(V0, U0, U1, U2)
    ! call pointInTri(U0, V0, V1, V2)

    ans = 0

end function coplanarTriTri

! subroutine edgeAgainstTriEdges(V0,V1,U0,U1,U2)
!
!   implicit none
!
!   real(kind=realType), dimension(3), intent(in) :: V0, V1, U0, U1, U2
!
!   float Ax,Ay,Bx,By,Cx,Cy,e,d,f;
!
!   Ax = V1(i0) - V0(i0)
!   Ay = V1(i1) - V0(i1)
!
!   EDGE_EDGE_TEST(V0,U0,U1)
!   EDGE_EDGE_TEST(V0,U1,U2)
!   EDGE_EDGE_TEST(V0,U2,U0)
!
! end subroutine edgeAgainstTriEdges
!
! function EDGE_EDGE_TEST(V0, U0, U1)
!
!   implicit none
!
!   real(kind=realType), dimension(3), intent(in) :: V0, U0, U1
!
!   Bx = U0(i0) - U1(i0)
!   By = U0(i1) - U1(i1)
!   Cx = V0(i0) - U0(i0)
!   Cy = V0(i1) - U0(i1)
!   f = Ay * Bx - Ax * By
!   d = By * Cx - Bx * Cy
!
!   if ((f .gt. 0 .and. d .ge. 0 .and. d .leq. f) .or. (f .lt. 0 .and. d .leq. 0 .and. d .geq. f))
!
!     e = Ax * Cy - Ay * Cx
!     if (f>0) then
!       if(e .geq. 0 .and. e .leq. f) return 1
!     else
!       if (e .leq. 0 .and. e .geq. f) return 1
!     end if
!   end if
!
!
! end function EDGE_EDGE_TEST
!
! function POINT_IN_TRI(V0,U0,U1,U2)
!
!   implicit none
!
!   real(kind=realType), dimension(3), intent(in) :: V0, U0, U1, U2
!
!   a = U1(i1) - U0(i1)
!   b = -(U1(i0) - U0(i0))
!   c = -a * U0(i0) - b * U0(i1)
!   d0 = a * V0(i0) + b * V0(i1) + c
!
!   a = U2(i1) - U1(i1)
!   b = -(U2(i0) - U1(i0))
!   c = -a * U1(i0) - b * U1(i1)
!   d1 = a * V0(i0) + b * V0(i1) + c
!
!   a = U0(i1) - U2(i1)
!   b = -(U0(i0) - U2(i0))
!   c = -a * U2(i0) - b * U2(i1)
!   d2 = a * V0(i0) + b * V0(i1) + c
!
!   if(d0 * d1 .gt. 0.0)
!     if(d0 * d2 .gt. 0.0) return 1
!   end if
!
! end function POINT_IN_TRI

subroutine cross_product(A, B, C)

  implicit none

  real(kind=realType), intent(in) :: A(3), B(3)
  real(kind=realType), intent(out) :: C(3)

  C(1) = A(2) * B(3) - A(3) * B(2)
  C(2) = A(3) * B(1) - A(1) * B(3)
  C(3) = A(1) * B(2) - A(2) * B(1)

end subroutine cross_product


end module Intersection
