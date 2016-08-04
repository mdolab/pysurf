program test

  use precision
  use Intersection
  implicit none

  real(kind=realType), dimension(:,:), allocatable :: coor

  real(kind=realType), dimension(3,2) :: BBox

  integer(kind=intType) :: intersect

  ! allocate(coor(3, 5))
  ! coor = reshape([-10, 3, 6, &
  !                 3, 144, 6, &
  !                 4, 3, 2, &
  !                 6, 7, 14, &
  !                 3, 2, 4] ,&
  !                 [3, 5])
  ! print *,coor
  ! call computeBoundingBox(coor, BBox)
  ! print *,BBox

  call triTriIntersect([0., 0., 0.], [1., 0., 0.], [0., 1., 0.],&
                  [.5, 0., 1.], [.5, 0., -1.], [.5, 1., 0.], intersect)

  print *,intersect

end program
