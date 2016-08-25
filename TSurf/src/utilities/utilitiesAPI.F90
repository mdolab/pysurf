module utilitiesAPI

use utilities
use precision
implicit none

! Here we define functions that can be generally used to handle FE data

contains

!=============================================

subroutine condenseBarNodes(nNodes, nBars, distTol, &
                            coor, barsConn, nUniqueNodes)

  ! This subroutine receives a list of bar FE which may have repeated points and
  ! Then condenses (merges) points that are closer than a given tolerance.
  ! This helps us getting a continuous line from the intersection algorithm results.
  ! Ney Secco 2016-08
  ! This is just an interface that allows us to call condenseBarFEs
  ! from Python when necessary.
  ! Remember to crop coorOut using nNodesOut in Python, because coorOut will
  ! have the same shape as coorIn.
  ! Ney Secco 2016-08

  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nNodes, nBars
  real(kind=realType), intent(in) :: distTol

  ! Input/Output variables
  real(kind=realType), dimension(3,nNodes), intent(inout) :: coor
  integer(kind=intType), dimension(2,nBars), intent(inout) :: barsConn

  ! Output variables
  integer(kind=intType), intent(out) :: nUniqueNodes

  !f2py intent(in) nNodes, nBars
  !f2py intent(in) distTol
  !f2py intent(inout) coor, barsConn
  !f2py intent(out) nUniqueNodes

  ! EXECUTION

  ! Just call the main routine that will do all the job
  call condenseBarNodes_main(nNodes, nBars, distTol, &
                             coor, barsConn, nUniqueNodes)

end subroutine condenseBarNodes

!=============================================

end module
