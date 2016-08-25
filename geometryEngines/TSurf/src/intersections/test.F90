program test

  use precision
  use Intersection
  implicit none

  real(kind=realType), dimension(:,:), allocatable :: coor

  real(kind=realType), dimension(3,2) :: BBox
  real(kind=realType) :: n, vecStart(3), vecEnd(3)

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

  ! ! 1: orthogonal triangles that do intersect
  ! call triTriIntersect([-2., 0., 0.], [2., 0., 0.], [0., 2., 0.],&
  !                 [-2., 1., -1.], [2., 1., -1.], [0., 1., 2.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 1.5: orthogonal triangles that do intersect
  ! call triTriIntersect([-2., 0., 0.]+[2., 0., 0.], [2., 0., 0.]+[2., 0., 0.], [0., 2., 0.]+[2., 0., 0.],&
  !                 [-2., 1., -1.]+[2., 0., 0.], [2., 1., -1.]+[2., 0., 0.], [0., 1., 2.]+[2., 0., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 1.75: orthogonal triangles that do intersect
  ! call triTriIntersect([-1., 0., 0.]+[2., 0., 0.], [1., 0., 0.]+[2., 0., 0.], [0., 1., 0.]+[2., 0., 0.],&
  !                 [-2., 1., -1.]+[2., 0., 0.], [2., 1., -1.]+[2., 0., 0.], [0., 1., 2.]+[2., 0., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''

  ! ! 2: orthogonal triangles that do not intersect
  ! call triTriIntersect([-1., 1., 10.], [1., 1., 10.], [0., -1., 10.],&
  !                 [-1., 0., 1.], [1., 0., 1.], [0., 0., -1.], intersect)
  ! if (0 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 3: coplanar triangles that do intersect
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-1., 0., 0.], [1., 0., 0.], [0., 1., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 4: coplanar triangles that do not intersect
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-1., -5., 0.], [1., -5., 0.], [0., -4., 0.], intersect)
  ! if (0 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 5: triangles that share an edge and do intersect
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-1., -1., 0.], [2., 2., 0.], [0., -4., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 6: triangles that share a point and do not intersect
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-1., -1., 0.], [-2., -2., 0.], [-10., -4., 0.], intersect)
  ! if (0 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 7: coplanar; triangle 1 is fully within triangle 2
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-2., 2., 0.], [2., 2., 0.], [0., -2., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 8: coplanar; triangle 2 is fully within triangle 1
  ! call triTriIntersect([-2., 2., 0.], [2., 2., 0.], [0., -2., 0.],&
  !                 [-1., 1., 0.], [1., 1., 0.], [0., -1., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''
  !
  ! ! 9: orthogonal triangles that do intersect
  ! call triTriIntersect([-1., 1., 0.], [1., 1., 0.], [0., -1., 0.],&
  !                 [-10., 0., 10.], [10., 0., 10.], [0., 0., -10.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''

  ! ! 10: noncoplanar triangles that share an edge
  ! call triTriIntersect([0., 0., 0.], [1., 0., 0.], [.5, 1., 1.],&
  !                 [0., 0., 0.], [1., 0., 0.], [.5, -1., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''

  ! ! 11: noncoplanar triangles that share an edge
  ! call triTriIntersect([-1.5, -2., 0.], [3., 4., 0.], [.5, 1., 1.],&
  !                 [0., 0., 0.], [6., 8., 0.], [.5, -1., 0.], intersect)
  ! if (1 .eq. intersect) then
  !   print *, 'good!'
  ! else
  !   print *, 'bad!'
  ! end if
  ! print *,''

  n = 0.
  ! 12: noncoplanar triangles that share an edge
  call triTriIntersect([-1.5, -2., 0.]+[0., n, 0.], [3., 4., 0.]+[0., n, 0.], [.5, 1., 1.]+[0., n, 0.],&
                  [0., 0., 0.]+[0., n, 0.], [6., 8., 0.]+[0., n, 0.], [.5, -1., 0.]+[0., n, 0.], intersect, vecStart, vecEnd)
  print *, vecStart
  print *, vecEnd

  if (1 .eq. intersect) then
    print *, 'good!'
  else
    print *, 'bad!'
  end if
  print *,''


end program
