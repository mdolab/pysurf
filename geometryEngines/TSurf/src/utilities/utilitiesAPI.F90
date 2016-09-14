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
  !
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

subroutine remesh(nNodes, nNewNodes, coor, barsConn, method, spacing, newCoor, newBarsConn)

  ! Remesh the given curve based on user-set options to obtain a better spacing
  ! John Jasa 2016-09

  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nNodes, nNewNodes
  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  integer(kind=intType), dimension(2,nNodes-1), intent(in) :: barsConn

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoor
  integer(kind=intType), dimension(2,nNewNodes-1), intent(out) :: newBarsConn

  !!@!!!!! FIX THIS RIGHT HERE

  !f2py intent(in) nNodes, nNewNodes
  !f2py intent(in) coor, barsConn
  !f2py intent(out) newCoor, newBarsConn

  ! EXECUTION

  ! Just call the main routine that will do the actual job
  call remesh_main(nNodes, nNewNodes, coor, barsConn, method, spacing, newCoor, newBarsConn)


end subroutine remesh

!=============================================

end module
