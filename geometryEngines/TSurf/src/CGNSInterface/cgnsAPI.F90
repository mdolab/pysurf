module CGNSapi

use precision
use CGNSinterface

! OUTPUTS
real(kind=realType), dimension(:, :), allocatable, save :: coor
integer(kind=intType), dimension(:, :), allocatable, save :: triaConn, quadsConn, barsConn
integer(kind=intType), dimension(:), allocatable, save :: surfTriaPtr, surfQuadsPtr, curveBarsPtr
character*32, dimension(:), allocatable, save :: surfNames, curveNames

! coor: real(3,numNodes) -> X,Y,Z coordinates of all nodes
! triaConn: real(3,numTria) -> Triangles connectivity
! quadsConn: integer(4,numQuads) -> Quads connectivity
! barsConn: integer(2,numBars) -> bars connectivity
! surfTriaPtr: integer(numSections) -> Pointer indicating index of triaConn where a surface section begins
! surfQuadsPtr: integer(numSections) -> Pointer indicating index of quadsConn where a surface section begins
! curveBarsPtr: integer(numCurves) -> Pointer indicating index of barsConn where a curve section begins

contains

subroutine readCGNS(cgns_file, comm, &
                    numCoor, numTriaConn, numQuadsConn, numBarsConn, &
                    numSurfTriaPtr, numSurfQuadsPtr, numCurveBarsPtr, &
                    numSurfNames, numCurveNames)

  implicit none

  !f2py intent(in) cgns_file, comm
  !f2py intent(out) numCoor, numTriaConn, numQuadsConn, numBarsConn
  !f2py intent(out) numSurfTriaPtr, numSurfQuadsPtr, numCurveBarsPtr
  !f2py intent(out) numSurfNames, numCurveNames

  ! INPUTS
  character(128), intent(in) :: cgns_file
  integer(kind=intType), intent(in) :: comm

  ! OUTPUTS
  integer(kind=intType), intent(out) :: numCoor, numTriaConn, numQuadsConn, numBarsConn
  integer(kind=intType), intent(out) :: numSurfTriaPtr, numSurfQuadsPtr, numCurveBarsPtr
  integer(kind=intType), intent(out) :: numSurfNames, numCurveNames

  ! Working variables
  integer(kind=intType) :: index

  call readCGNSmain(cgns_file, comm, coor, triaConn, quadsConn, barsConn, &
                    surfTriaPtr, surfQuadsPtr, curveBarsPtr, &
                    surfNames, curveNames)

  ! Get the sizes of the allocated variables to return them to Python
  numCoor = size(coor, 2)
  numTriaConn = size(triaConn, 2)
  numQuadsConn = size(quadsConn, 2)
  numBarsConn = size(barsConn, 2)
  numSurfTriaPtr = size(surfTriaPtr)
  numSurfQuadsPtr = size(surfQuadsPtr)
  numCurveBarsPtr = size(curveBarsPtr)
  numSurfNames = size(surfNames)
  numCurveNames = size(curveNames)

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

subroutine retrieveData(numCoor, numTriaConn, numQuadsConn, numBarsConn, &
                        numSurfTriaPtr, numSurfQuadsPtr, numCurveBarsPtr, &
                        numSurfNames, numCurveNames, &
                        coorData, triaConnData, quadsConnData, barsConnData, &
                        surfTriaPtrData, surfQuadsPtrData, curveBarsPtrData, &
                        surfNamesData, curveNamesData)

  ! This is a function used to retrieve data from the allocatable arrays back to Python, so
  ! that Python does not need a direct access to the allocatable variables.

  ! INPUTS
  integer(kind=intType), intent(in) :: numCoor, numTriaConn, numQuadsConn, numBarsConn
  integer(kind=intType), intent(in) :: numSurfTriaPtr, numSurfQuadsPtr, numCurveBarsPtr
  integer(kind=intType), intent(in) :: numSurfNames, numCurveNames  

  ! OUTPUTS
  real(kind=realType), intent(out), dimension(3, numCoor) :: coorData
  integer(kind=intType), intent(out), dimension(3, numTriaConn) :: triaConnData
  integer(kind=intType), intent(out), dimension(4, numQuadsConn) :: quadsConnData
  integer(kind=intType), intent(out), dimension(2, numBarsConn) :: barsConnData
  integer(kind=intType), intent(out), dimension(numSurfTriaPtr) :: surfTriaPtrData
  integer(kind=intType), intent(out), dimension(numSurfQuadsPtr) :: surfQuadsPtrData
  integer(kind=intType), intent(out), dimension(numCurveBarsPtr) :: curveBarsPtrData
  character*32, intent(out), dimension(numSurfNames) :: surfNamesData
  character*32, intent(out), dimension(numCurveNames) :: curveNamesData

  ! EXECUTION
  
  ! Get data from the allocatable variables
  coorData = coor
  triaConnData = triaConn
  barsConnData = barsConn
  surfTriaPtrData = surfTriaPtr
  surfQuadsPtrData = surfQuadsPtr
  curveBarsPtrData = curveBarsPtr
  surfNamesData = surfNames
  curveNamesData = curveNames
  
end subroutine retrieveData

subroutine releaseMemory()

  ! This subroutine just deallocates memory used by the CGNSinterface.
  ! Remember to call this function after you copy the outputs in Python.

  use precision
  use CGNSGrid
  implicit none

  ! Working variables
  integer(kind=intType) :: iZone, sec, nZones

  ! Deallocate zones
  if (allocated(zones)) then
     do iZone=1, nZones
       do sec=1, zones(iZone)%nSections
          deallocate(zones(iZone)%sections(sec)%elemConn)
          deallocate(zones(iZone)%sections(sec)%elemPtr)
        end do
        deallocate(zones(iZone)%sections)
     end do
     deallocate(zones)
  end if

  ! Deallocate output variables
  if (allocated(coor)) then
     deallocate(coor, triaConn, quadsConn, barsConn, &
                surfTriaPtr, surfQuadsPtr, curveBarsPtr, &
                surfNames, curveNames)
  end if

end subroutine releaseMemory

end module CGNSapi
