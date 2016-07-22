module readCGNS

contains

subroutine readCGNS_routine(cgns_file, coor, triaConn, quadsConn)
  use communication
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file

  ! Output
  integer(kind=intType), intent(out), dimension(:,:), allocatable :: triaConn, quadsConn
  real(kind=realType), intent(out), dimension(:,:), pointer :: coor

  ! Working
  integer(kind=intType) :: cg, ierr, i
  integer(kind=intType):: CellDim, PhysDim, nZones, base, nbases
  integer(kind=intType) :: nstructured, nunstructured, zoneType
  integer(kind=intType) :: iZone, sec, nElem, nConn, iEnd, iStart, nQuads, nTria
  integer(kind=intType) :: nTriaTotal, nQuadsTotal, iQuads, iTria, dims(2)

  character*32 :: baseName

  ! Set the default family names
  defaultFamName(BCAxisymmetricWedge) = 'axi'
  defaultFamName(BCDegenerateLine) = 'degenerate'
  defaultFamName(BCDegeneratePoint) ='degenerate'
  defaultFamName(BCDirichlet) = 'dirichlet'
  defaultFamName(BCExtrapolate) = 'extrap'
  defaultFamName(BCFarfield) = 'far'
  defaultFamName(BCGeneral) = 'general'
  defaultFamName(BCInflow) = 'inflow'
  defaultFamName(BCInflowSubsonic) = 'inflow'
  defaultFamName(BCInflowSupersonic) = 'inflow'
  defaultFamName(BCNeumann) = 'neumann'
  defaultFamName(BCOutflow) = 'outflow'
  defaultFamName(BCOutflowSubsonic) = 'outflow'
  defaultFamName(BCOutflowSupersonic)  ='outflow'
  defaultFamName(BCSymmetryPlane) = 'sym'
  defaultFamName(BCSymmetryPolar) = 'sympolar'
  defaultFamName(BCTunnelInflow) = 'inflow'
  defaultFamName(BCTunnelOutflow) = 'outflow'
  defaultFamName(BCWall) = 'wall'
  defaultFamName(BCWallInviscid) = 'wall'
  defaultFamName(BCWallViscous) = 'wall'
  defaultFamName(BCWallViscousHeatFlux) = 'wall'
  defaultFamName(BCWallViscousIsothermal) = 'wall'
  defaultFamName(FamilySpecified) = 'wall'

  ! Do the I/O that is common to both types of grids

  if (myid == 0) then
     print *, ' -> Reading CGNS File: ', cgns_file

     ! Open and get the number of zones:
     call cg_open_f(trim(cgns_file), CG_MODE_READ, cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     call cg_nbases_f(cg, nbases, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f

     if (nbases .gt. 1) then
        print *, ' ** Warning: pyWarpUstruct only reads the first base in a cgns file'
     end if

     base = 1_intType

     call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     if (cellDim .ne. 2 .or. PhysDim .ne. 3) then
        print *, 'The Cells must 3 dimensional'
        stop
     end if

     call cg_nzones_f(cg, base, nZones, ierr);
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     print *, '   -> Number of Zones:', nzones

     ! Determine if we have structured or unstructured zones. We can
     ! only deal with one or the other.
     nStructured = 0
     nUnstructured = 0
     do i=1, nZones
        call cg_zone_type_f(cg, base, i, zoneType, ierr)
        if (zoneType == Structured) then
           nStructured = nStructured + 1
        else if (zoneType == Unstructured) then
           nUnstructured = nUnstructured + 1
        end if
     end do


    call readUnstructuredCGNS(cg)

    nTriaTotal = 0
    nQuadsTotal = 0

    ! Loop over the zones to obtain the total number of trias and quads
    zoneLoop1: do iZone = 1,nZones
       do sec=1, zones(iZone)%nSections
         nTria = zones(iZone)%sections(sec)%nTria
         nQuads = zones(iZone)%sections(sec)%nQuads
         nTriaTotal = nTriaTotal + nTria
         nQuadsTotal = nQuadsTotal + nQuads
      end do
    end do zoneLoop1

    iTria = 1
    iQuads = 1
    allocate(triaConn(3, nTriaTotal))
    allocate(quadsConn(4, nQuadsTotal))

    ! Loop over the zones and read the nodes
    zoneLoop2: do iZone = 1,nZones
       do sec=1, zones(iZone)%nSections
         nElem = zones(iZone)%sections(sec)%nElem
         nTria = zones(iZone)%sections(sec)%nTria
         nQuads = zones(iZone)%sections(sec)%nQuads

         if (nTria + nQuads .ne. 0) then
           do i=2, nElem+1
             iStart = zones(iZone)%sections(sec)%elemPtr(i-1)
             iEnd = zones(iZone)%sections(sec)%elemPtr(i)
             nConn = iEnd - iStart
             if (nConn .eq. 3) then
               triaConn(:, iTria) = zones(iZone)%sections(sec)%elemConn(iStart:iEnd-1)
               iTria = iTria + 1
             else
               quadsConn(:, iQuads) = zones(iZone)%sections(sec)%elemConn(iStart:iEnd-1)
               iQuads = iQuads + 1
             end if
           end do
        end if
      end do
    end do zoneLoop2
  end if

  dims = shape(allNodes)

  allocate(coor(dims(1), dims(2)))
  coor => allNodes

end subroutine readCGNS_routine

end module readCGNS
