module CGNSf2py

use precision
use CGNSinterface

!real(kind=realType), dimension(:,:), allocatable :: coor
!integer(kind=intType), dimension(:,:), allocatable :: triaConn, quadsConn

real(kind=realType), dimension(3,3) :: coor2
integer(kind=intType), dimension(3,3) :: triaConn2, quadsConn2

contains

subroutine CGNSf2py_routine(cgns_file, comm)

! INPUTS
character(32), intent(in) :: cgns_file
integer(kind=intType), intent(in) :: comm

real(kind=realType), dimension(:,:), allocatable :: coor
integer(kind=intType), dimension(:,:), allocatable :: triaConn, quadsConn

!f2py intent(in) cgns_file, comm

call readCGNS(cgns_file, comm, coor, triaConn, quadsConn)

end subroutine CGNSf2py_routine

end module CGNSf2py
