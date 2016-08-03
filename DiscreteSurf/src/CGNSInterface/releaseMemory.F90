module releaseMemory

contains

subroutine release

  use precision
  use CGNSGrid
  use CGNSInterface
  implicit none

  ! Working variables
  integer(kind=intType) :: iZone, sec, nZones

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


end subroutine release

end module
