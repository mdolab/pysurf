module utilitiesAPI

use utilities
use precision
implicit none

! Here we define functions that can be generally used to handle FE data

contains

!=============================================

subroutine condenseBarNodes(nNodes, nBars, distTol, &
                            coor, barsConn, nUniqueNodes, linkOld2New)

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
  integer(kind=intType), dimension(nNodes), intent(out) :: linkOld2New

  !f2py intent(in) nNodes, nBars
  !f2py intent(in) distTol
  !f2py intent(inout) coor, barsConn
  !f2py intent(out) nUniqueNodes, linkOld2New

  ! EXECUTION

  ! Just call the main routine that will do all the job
  call condenseBarNodes_main(nNodes, nBars, distTol, &
                             coor, barsConn, nUniqueNodes, linkOld2New)

end subroutine condenseBarNodes

subroutine remesh(nNodes, nElem, nNewNodes, coor, barsConn, method, spacing,&
  initialSpacing, finalSpacing, newCoor, newBarsConn)

  ! Remesh the given curve based on user-set options to obtain a better spacing
  ! John Jasa 2016-09

  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nElem, nNewNodes
  integer(kind=intType), intent(inout) :: nNodes

  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  integer(kind=intType), dimension(2,nElem), intent(in) :: barsConn
  real(kind=realType), intent(in) :: initialSpacing, finalSpacing

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoor
  integer(kind=intType), dimension(2,nNewNodes-1), intent(out) :: newBarsConn

  !f2py intent(in) nNodes, nNewNodes
  !f2py intent(in) coor, barsConn
  !f2py intent(out) newCoor, newBarsConn

  ! EXECUTION

  ! Just call the main routine that will do the actual job
  call remesh_main(nNodes, nElem, nNewNodes, coor, barsConn, method, spacing,&
    initialSpacing, finalSpacing, newCoor, newBarsConn)

end subroutine remesh

subroutine remesh_b(nNodes, nElem, nNewNodes, nNewElems, coor, newCoorb, barsConn,&
  method, spacing, initialSpacing, finalSpacing, newCoor, newBarsConn, coorb)

  ! Remesh the given curve based on user-set options to obtain a better spacing
  ! John Jasa 2016-09

  use utilities_b, only: remesh_main_b
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nElem, nNewNodes, nNewElems
  integer(kind=intType), intent(inout) :: nNodes
  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  real(kind=realType), dimension(3,nNewNodes), intent(in) :: newCoorb
  integer(kind=intType), dimension(2,nElem), intent(in) :: barsConn
  real(kind=realType), intent(in) :: initialSpacing, finalSpacing

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoor
  real(kind=realType), dimension(3,nNodes), intent(out) :: coorb
  integer(kind=intType), dimension(2,nNewElems), intent(out) :: newBarsConn

  integer(kind=intType) :: i

  !f2py intent(in) nNodes, nNewNodes
  !f2py intent(in) coor, barsConn, newCoorb
  !f2py intent(out) newCoor, newBarsConn, coorb

  ! EXECUTION

  ! Just call the main routine that will do the actual job
  call remesh_main_b(nNodes, nElem, nNewNodes, coor, coorb, barsConn, method,&
  spacing, initialSpacing, finalSpacing, newCoor, newCoorb, newBarsConn)

end subroutine remesh_b

subroutine remesh_d(nNodes, nElem, nNewNodes, nNewElems, coor, coord, barsConn,&
  method, spacing, initialSpacing, finalSpacing, newCoor, newCoord, newBarsConn)

  ! Remesh the given curve based on user-set options to obtain a better spacing
  ! John Jasa 2016-09

  use utilities_d, only: remesh_main_d
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nElem, nNewNodes, nNewElems
  integer(kind=intType), intent(inout) :: nNodes
  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  real(kind=realType), dimension(3,nNodes), intent(in) :: coord
  integer(kind=intType), dimension(2,nElem), intent(in) :: barsConn
  real(kind=realType), intent(in) :: initialSpacing, finalSpacing

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoor
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoord
  integer(kind=intType), dimension(2,nNewElems), intent(out) :: newBarsConn

  !f2py intent(in) nNodes, nNewNodes
  !f2py intent(in) coor, barsConn, coord
  !f2py intent(out) newCoor, newBarsConn, newCoord

  ! EXECUTION

  ! Just call the main routine that will do the actual job
  call remesh_main_d(nNodes, nElem, nNewNodes, coor, coord, barsConn, method,&
  spacing, initialSpacing, finalSpacing, newCoor, newCoord, newBarsConn)

end subroutine remesh_d

!=============================================

end module
