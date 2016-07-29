module CGNSapi

use precision
use CGNSinterface

! OUTPUTS
real(kind=realType), dimension(:, :), allocatable :: coor
integer(kind=intType), dimension(:, :), allocatable :: triaConn, quadsConn, barsConn
integer(kind=intType), dimension(:), allocatable :: surfTriaPtr, surfQuadsPtr, curveBarsPtr
character*32, dimension(:), allocatable :: surfNames, curveNames

! coor: real(3,numNodes) -> X,Y,Z coordinates of all nodes
! triaConn: real(3,numTria) -> Triangles connectivity
! quadsConn: integer(4,numQuads) -> Quads connectivity
! barsConn: integer(2,numBars) -> bars connectivity
! surfTriaPtr: integer(numSections) -> Pointer indicating index of triaConn where a surface section begins
! surfQuadsPtr: integer(numSections) -> Pointer indicating index of quadsConn where a surface section begins
! curveBarsPtr: integer(numCurves) -> Pointer indicating index of barsConn where a curve section begins

contains

subroutine readCGNS(cgns_file, comm)

  implicit none

  !f2py intent(in) cgns_file, comm

  ! INPUTS
  character(128), intent(in) :: cgns_file
  integer(kind=intType), intent(in) :: comm

  ! Working variables
  integer(kind=intType) :: index

  call readCGNSmain(cgns_file, comm, coor, triaConn, quadsConn, barsConn, &
                    surfTriaPtr, surfQuadsPtr, curveBarsPtr, &
                    surfNames, curveNames)

  ! Print log
  print *,'The following sections were found'
  print *,'-> Surfaces:'
  do index = 1,size(surfNames)
     print *,surfNames(index)
  end do
  print *,'-> Curves:'
  do index = 1,size(curveNames)
     print *,curveNames(index)
  end do

end subroutine readCGNS

end module CGNSapi
