module CGNSf2py

use precision
use CGNSinterface

! OUTPUTS
real(kind=realType), dimension(:, :), allocatable :: coor
integer(kind=intType), dimension(:, :), allocatable :: triaConn, quadsConn


contains

subroutine CGNSf2py_routine(cgns_file, comm)

  implicit none

  !f2py intent(in) cgns_file, comm

  ! INPUTS
  character(32), intent(in) :: cgns_file
  integer(kind=intType), intent(in) :: comm

  call readCGNS(cgns_file, comm, coor, triaConn, quadsConn)

end subroutine CGNSf2py_routine

end module
