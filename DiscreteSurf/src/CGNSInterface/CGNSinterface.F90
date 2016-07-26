module CGNSinterface

contains

subroutine readCGNS(cgns_file, comm, coor, triaConn, quadsConn)

  ! This is the main subroutine that should be called to read a CGNS file
  ! INPUTS
  ! cgns_file: string -> name of the CGNS file
  ! comm: integer -> MPI communicator
  !
  ! OUTPUTS
  ! coor: real(3,numNodes) -> X,Y,Z coordinates of all nodes
  ! triaConn: real(3,numTria) -> Triangles connectivity
  ! quadsConn: real(4,numQuads) -> Quads connectivity

  use mpi
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  character*(*),intent(in) :: cgns_file
  integer, intent(in) :: comm

  ! Output
  integer(kind=intType), intent(out), dimension(:,:), allocatable :: triaConn, quadsConn
  real(kind=realType), intent(out), dimension(:,:), allocatable :: coor

  ! Working
  integer(kind=intType) :: cg, ierr, i, myid
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

  ! Determine current MPI rank
  call MPI_COMM_RANK(comm, myid, ierr)

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
        print *, 'The Cells must 2 dimensional'
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

     ! Call the read for unstructured meshes
     call readUnstructuredCGNS(cg, coor)

     ! From now on we need to create connectivity arrays in the format
     ! required by ADT. That is one array for triangle connectivities
     ! and another array for quad connectivities.

     ! Initialize counter for different types of elements
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

end subroutine readCGNS

subroutine readUnstructuredCGNS(cg,allNodes)

  use precision
  use communication
  use CGNSGrid

  implicit none
  include 'cgnslib_f.h'

  ! Input Arguments
  integer(kind=intType), intent(in) :: cg

  ! CGNS Variables
  integer(kind=intType) :: i, j, k, m, l, istart, iend, localsize, iProc
  integer(kind=intType):: ierr, base, dims(3), iZone
  integer(kind=intType):: nNodes, nCells
  integer(kind=intType) :: tmpSym, nSymm
  character*32 :: zoneName, bocoName, famName
  character*32 :: secName

  integer(kind=intType) :: nbocos, boco
  integer(kind=intType) :: nVertices,nElements, nzones
  integer(kind=intType) :: zoneType,dataType,sec,type
  integer(kind=intType) :: nSections,nElem,nConn
  integer(kind=intType) :: eBeg,eEnd,nBdry, parentFlag
  integer(kind=intType) :: bocoType,ptsettype,nbcelem
  integer(kind=intType) :: normalIndex,normalListFlag
  integer(kind=intType) :: nDataSet, tmpInt(2)
  integer(kind=intType) :: elementDataSize, curElement, eCount, nPnts
  real(kind=realType), dimension(:), allocatable :: coorX, coorY, coorZ
  integer(kind=intType), dimension(:), allocatable :: tmpConn
#ifdef USE_COMPLEX
  complex(kind=realType), intent(out), dimension(:, :), allocatable, target :: allNodes
  complex(kind=realType), intent(out), dimension(:, :), allocatable :: localNodes
  complex(kind=realType), dimension(:,:), pointer :: elemNodes
#else
  real(kind=realType), dimension(:, :), allocatable, target :: allNodes
  real(kind=realType), dimension(:, :), allocatable :: localNodes
  real(kind=realType), dimension(:,:), pointer :: elemNodes
#endif
  integer(kind=intType), dimension(:), allocatable :: surfaceNodes, localSurfaceNodes
  integer(kind=intType), dimension(:, :), allocatable :: sizes

  integer(kind=intType) :: status(MPI_STATUS_SIZE)
  type(sectionDataType), pointer :: secPtr
  integer(kind=intType), dimension(:), pointer :: elemConn, elemPtr

  integer(kind=intType) :: nElemNotFound, curElem, curElemSize, localIndex, zoneStart
  real(kind=realType) :: symmSum(3)
  logical :: isSurface

  ! ---------------------------------------
  !           Open CGNS File
  ! ---------------------------------------
  if (myID == 0) then

     base = 1_intType

     ! read the number of zones in the file
     call cg_nzones_f(cg, base, nzones, ierr)
     if (ierr .eq. ERROR) call cg_error_exit_f

     ! allocate the zonal storage
     allocate(zones(nZones))

     ! do a check that we have all unstructured zones
     do iZone = 1,nZones
        ! Check the zone type. This should be Unstructured.
        call cg_zone_type_f(cg,base, iZone, zoneType, ierr)
        if (ierr .eq. ERROR) call cg_error_exit_f

        if (zoneType .ne. Unstructured) then
           print *,"Unstructured zone expected...stopping"
           stop
        end if
     end do

     ! Now count up the total number of nodes
     nNodes = 0
     nCells = 0
     allocate(sizes(3, nZones))

     do iZone=1, nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        nVertices = dims(1)
        nElements = dims(2)
        zones(iZone)%nVertices = nVertices
        zones(iZone)%nElements = nElements
        zones(iZone)%name = zoneName
        sizes(:, iZone) = dims(1:3)
        nNodes = nNodes + dims(1)
        nCells = nCells + dims(2)
     end do

     ! Now we know the total number of nodes we can allocate the final
     ! required space and read them in.
     allocate(allNodes(3, nNodes))
     allocate(surfaceNodes(nNodes))
     surfaceNodes = 0
     zoneStart = 0

     ! Loop over the zones and read the nodes
     zoneLoop: do iZone = 1,nZones
        call cg_zone_read_f(cg, base, iZone, zoneName, dims, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        allocate(coorX(dims(1)), coorY(dims(1)), coorZ(dims(1)))

        ! Read the x,y,z-coordinates. Assume double precision
        call cg_coord_read_f(cg,base,iZone,"CoordinateX",&
             & realDouble, 1, dims(1), coorX, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_read_f(cg,base,iZone,"CoordinateY",&
             & realDouble, 1, dims(1), coorY, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_coord_read_f(cg,base,iZone,"CoordinateZ",&
             & realDouble, 1, dims(1), coorZ, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Now stack all of the zones in one array
        do i=1, dims(1)
#ifdef USE_COMPLEX
           allNodes(1, i+zoneStart) = cmplx(coorX(i), 0.0)
           allNodes(2, i+zoneStart) = cmplx(coorY(i), 0.0)
           allNodes(3, i+zoneStart) = cmplx(coorZ(i), 0.0)
#else
           allNodes(1, i+zoneStart) = coorX(i)
           allNodes(2, i+zoneStart) = coorY(i)
           allNodes(3, i+zoneStart) = coorZ(i)
#endif
        end do

        ! Now loop over the sections in the unstructured file and
        ! figure out the boundary nodes.

        ! Determine the total number of sections in this zone. This
        ! will involve surface sections and volume sections. We only
        ! care about the surface sections
        call cg_nsections_f(cg, base, iZone, nSections, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        zones(iZone)%nSections = nSections

        ! Allocate space for the sections
        allocate(zones(iZone)%sections(nSections))

        do sec=1, nSections
           ! Read the section and determine the element type and store
           ! some generic information.
           call cg_section_read_f(cg, base, iZone, sec, secName, &
                type, eBeg, eEnd, nBdry, parentFlag, ierr)
           if (ierr .eq. CG_ERROR) call cg_error_exit_f

           ! Nullify the elemPtr and elemConn, since it may not be allocated
           nullify(zones(iZone)%sections(sec)%elemConn, &
                zones(iZone)%sections(sec)%elemPtr)

           ! Number of elements on this section
           nElem = eEnd - eBeg + 1
           zones(iZone)%sections(sec)%nElem = nElem
           zones(iZone)%sections(sec)%name = secName
           zones(iZone)%sections(sec)%elemStart = eBeg
           zones(iZone)%sections(sec)%elemEnd = eEnd
           zones(iZone)%sections(sec)%nQuads = 0
           zones(iZone)%sections(sec)%nTria = 0

           ! Currently we can only deal with three and four sided
           ! surfaces. These can show up type == TRI_3 (all tris),
           ! type == QUAD_4 or type == MIXED. In the case of mixed we
           ! need to check each of the types individually.
           select case (type)

              ! all tris or quads
           case (TRI_3, QUAD_4)
              ! Firstly flag this section as a surface
              zones(iZone)%sections(sec)%isSurface = .True.

              ! Get the nodes per element "npe" (constant)
              call cg_npe_f(type, nConn, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Allocate space for the connectivities and the pointer
              allocate(zones(iZone)%sections(sec)%elemConn(nConn*nElem))
              allocate(zones(iZone)%sections(sec)%elemPtr(nElem+1))

              ! This is the actual connectivity real call.
              call cg_elements_read_f(cg, base, iZone, sec, &
                   zones(iZone)%sections(sec)%elemConn, NULL, ierr)
              zones(iZone)%sections(sec)%elemConn = &
                   zones(iZone)%sections(sec)%elemConn + zoneStart
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Set up the pointer which is simple in this case...it
              ! is just an arithematic list.
              zones(iZone)%sections(sec)%elemPtr(1) = 1
              do i=2,nElem+1
                 zones(iZone)%sections(sec)%elemPtr(i) = &
                      zones(iZone)%sections(sec)%elemPtr(i-1) + nConn
              end do

              if (nConn .eq. 3) then
                 zones(iZone)%sections(sec)%nTria = nElem
              else
                 zones(iZone)%sections(sec)%nQuads = nElem
              end if

           case (MIXED)

              ! This section *could* be a surface but we don't
              ! actually know yet.

              ! First get the number of values we need to read
              call cg_ElementDataSize_f(cg, base, iZone, sec, &
                   ElementDataSize , ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! Allocate a temporary for this
              allocate(tmpConn(ElementDataSize))

              ! Now read the 'connectivity'--- not really, it has the
              ! connectivity *AND* the element types.
              call cg_elements_read_f(cg, base, iZone, sec, &
                   tmpConn, NULL, ierr)
              if (ierr .eq. CG_ERROR) call cg_error_exit_f

              ! NOW we can figure out if this is actually a surface
              ! zone. Assume it is is, and flag and stop if not.
              zones(iZone)%sections(sec)%isSurface = .True.
              i = 1
              elementCheckLoop: do while (i <= ElementDataSize)
                 curElement = tmpConn(i)
                 select case (curElement)
                 case (TRI_3, QUAD_4)
                    ! We're OK...now determine how much we need to
                    ! increment the counter

                    call cg_npe_f(curElement, nConn, ierr)
                    if (ierr .eq. CG_ERROR) call cg_error_exit_f

                    i = i + nConn + 1 ! The plus 1 is for the element Type

                 case default
                    ! Not a surface or we can't read it.
                    zones(iZone)%sections(sec)%isSurface = .False.
                    exit elementCheckLoop
                 end select
              end do elementCheckLoop

              ! Only continue if we are now sure the section is only
              ! surfaces
              if (zones(iZone)%sections(sec)%isSurface) then

                 ! Allocate space for the connectivities and the
                 ! pointer. It should be safe to take this as
                 ! ElementDataSize - nElem.
                 allocate(zones(iZone)%sections(sec)%elemConn(ElementDataSize-nElem))
                 allocate(zones(iZone)%sections(sec)%elemPtr(nElem+1))

                 zones(iZone)%sections(sec)%elemConn = -999999999

                 ! Loop back over the section elements
                 ! again. Essentially we are taking the information in
                 ! tmpConn and splitting up into elemConn and elemPtr
                 i = 1
                 eCount = 1
                 zones(iZone)%sections(sec)%elemPtr(eCount) = 1
                 elementLoop: do while (i <= ElementDataSize)
                    curElement = tmpConn(i)
                    select case (curElement)
                    case (TRI_3, QUAD_4)
                       ! We're OK...now determine how much we need to
                       ! increment the counter

                       call cg_npe_f(curElement, nConn, ierr)
                       if (ierr .eq. CG_ERROR) call cg_error_exit_f

                       ! iStart is the start of the element in elemPtr
                       ! iEnd is one more than the end of the element
                       ! in elemPtr (or the start of the next one)

                       iStart = zones(iZone)%sections(sec)%elemPtr(eCount)
                       iEnd   = iStart + nConn
                       zones(iZone)%sections(sec)%elemPtr(eCount + 1) = iEnd

                       ! The connectivity is easy:
                       zones(iZone)%sections(sec)%elemConn(iStart:iEnd-1) = &
                            tmpConn(i+1: i+nConn) + zoneStart

                       ! Move to the next element
                       i = i + nConn + 1
                       eCount = eCount + 1

                       if (nConn .eq. 3) then
                          zones(iZone)%sections(sec)%nTria = &
                               zones(iZone)%sections(sec)%nTria + 1
                       else
                          zones(iZone)%sections(sec)%nQuads = &
                               zones(iZone)%sections(sec)%nQuads + 1
                       end if

                    case default
                       ! Not a surface or we can't read it.
                       zones(iZone)%sections(sec)%isSurface = .False.
                       exit elementLoop
                    end select
                 end do elementLoop

              end if

              ! Delete the temporary CGNS connectivity infor
              deallocate(tmpConn)

           case default
              ! Flag as not a surface
              zones(iZone)%sections(sec)%isSurface = .False.
           end select

        end do

        ! Now below we will be assuming that the elements are ordered
        ! sequentially between the muliple sections and that there are
        ! NO gaps. This is required for efficient searching. If this
        ! *isn't* the case, we will just stop. For future reference,
        ! we could make a pointer array to the sections that *is*
        ! sorted to get around this issue and make a compact pointer
        ! for the elements as well.
        do sec=1, nSections-1
           if ( zones(iZone)%sections(sec+1)%elemStart /= &
                zones(iZone)%sections(sec  )%elemEnd + 1) then
              print *,'The section element ranges are not in order. This cannot be handled.'
              stop
           end if
        end do

        ! Determine the number of boundary conditions:
        call cg_nbocos_f(cg, base, iZone, nBocos, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Allocate space for bocos
        allocate(zones(iZone)%bocos(nBocos))

        ! On the first pass we just want to figure out the number of
        ! *surface* nodes and the size we need for connectivity. This
        ! is somehwat "dumb" since we are just going to take the
        ! number of elements * number of nodes per element, which will
        ! obviously overpredict by a large fraction. That's not a huge
        ! deal.

        ! bocoLoop: do boco=1, nBocos
        !
        !    ! Read the info for this boundary condition.
        !    call cg_boco_info_f(cg, base, iZone, boco, boconame, bocotype, &
        !         ptsettype, nPnts, NormalIndex, NormalListFlag, datatype, &
        !         ndataset, ierr)
        !    if (ierr .eq. CG_ERROR) call cg_error_exit_f
        !
        !    ! We only allow ElementRange and ElementList...its just too
        !    ! hard to figure out what is going on with the point types.
        !    if (ptSetType == PointRange .or. ptSetType == pointList) then
        !       print *, 'pyWarpUStruct cannot handle boundary conditions &
        !            &defined by "PointRange" or "PointList". Please use &
        !            &boundary conditions defined by "ElementRange" or "ElementList" &
        !            &instead'
        !
        !       stop
        !    end if
        !
        !    if (ptSetType == ElementRange) then
        !       ! We have to read off the start and end elements, this
        !       ! will always require an array of lenght 2 (tmpInt)
        !
        !       call cg_boco_read_f(cg, base, iZone, boco, tmpInt, NULL, ierr)
        !       if (ierr .eq. CG_ERROR) call cg_error_exit_f
        !       nBCElem = tmpInt(2) - tmpInt(1) + 1
        !       allocate(zones(iZone)%bocos(boco)%BCElements(nBCElem))
        !
        !       ! We now need to fill it up trivially
        !       do i=1,nBCElem
        !          zones(iZone)%bocos(boco)%BCElements(i) = tmpInt(1) + i - 1
        !       end do
        !    else if (ptSetType == ElementList) then
        !
        !       ! nPnts is the actual number of elements we can allocate
        !       ! and read directly
        !
        !       allocate(zones(iZone)%bocos(boco)%BCElements(nPnts))
        !       call cg_boco_read_f(cg, base, iZone, boco, &
        !            zones(iZone)%bocos(boco)%BCElements, NULL, ierr)
        !       nBCElem = nPnts
        !    end if
        !
        !    ! Now we need to deduce the actual connectivity for this
        !    ! boundary condition. This will require querying the
        !    ! section info. Since the same BC could reference two
        !    ! sections (i think in theory)
        !
        !    ! On the first pass we find the sectionID and index of the
        !    ! element ON THAT SECTION. This require
        !
        !    ! We will assume at this point that we have at most quads,
        !    ! so we can determine the maximum size of elemPr and elemConn
        !    allocate(zones(iZone)%bocos(boco)%elemConn(4*nBCElem))
        !    allocate(zones(iZone)%bocos(boco)%elemPtr(nBCElem + 1))
        !    allocate(zones(iZOne)%bocos(boco)%elemNodes(3, 4*nBCElem))
        !
        !    ! Set pointers for easier reading
        !    elemConn => zones(iZone)%bocos(boco)%elemConn
        !    elemPtr  => zones(iZone)%bocos(boco)%elemPtr
        !    elemNodes=> zones(iZone)%bocos(boco)%elemNodes
        !    elemConn = -1
        !    elemPtr = 0
        !    elemPtr(1) = 1
        !    nElemNotFound = 0
        !    k = 1 ! Running index of elemConn, elemNodes
        !    l = 1 ! Running index of elemPtr
        !
        !    do i=1, size(zones(iZone)%bocos(boco)%bcElements)
        !
        !       curElem = zones(iZone)%bocos(boco)%bcElements(i)
        !       sec = -1
        !
        !       ! For now do a linear search
        !       setionLoop: do j=1,size(zones(iZone)%sections)
        !          if (curElem >= zones(iZone)%sections(j)%elemStart .and. &
        !               curElem <= zones(iZone)%sections(j)%elemEnd) then
        !             sec = j
        !             exit setionLOop
        !          end if
        !       end do setionLoop
        !
        !       if (sec == -1)  then
        !          nElemNotFound = nElemNotFound + 1
        !          ! No element...do not increment the k-counter
        !       else
        !          ! Set pointer for easier code reading
        !          secPtr => zones(iZone)%sections(sec)
        !
        !          ! localIndex is the index on the section
        !          localIndex = curElem - secPtr%elemStart + 1
        !
        !          ! iStart/iEnd are the indices on the section
        !          iStart = secPtr%elemPtr(localIndex)
        !          iEnd   = secPtr%elemPtr(localIndex+1) - 1
        !
        !          ! Number of nodes on this section.
        !          curElemSize =  iEnd - iStart + 1
        !
        !          ! Since we just blew up the global connectivity, the
        !          ! local connectivity is now easy...it is just
        !          ! 1,2,3,4...etc. elemPtr will determine the type of
        !          ! element (3 or 4 nodes). Now copy the actual nodes
        !          do m=iStart,iEnd
        !             elemNodes(:, k) = allNodes(:, secPtr%elemConn(m))
        !             elemConn(k) = k
        !             k = k + 1
        !          end do
        !
        !          ! Update the elemPtr too
        !          elemPtr(l+1) = elemPtr(l) + curElemSize
        !          l = l + 1
        !
        !          ! Flag surface nodes
        !          do m=iStart, iEnd
        !             surfaceNodes(secPtr%elemConn(m)) = 1
        !          end do
        !       end if
        !    end do
        !
        !    ! Now we set the *actual* number of nBCElem and nBCNodes
        !    zones(iZone)%bocos(boco)%nBCElem = nBCElem - nElemNotFound
        !    zones(iZOne)%bocos(boco)%nBCNodes = k-1
        !
        !    ! Save the boco name
        !    zones(iZone)%bocos(boco)%name =  bocoName
        !    zones(iZone)%bocos(boco)%type = bocoType
        !    ! By Default no family name
        !    zones(iZone)%bocos(boco)%family = ""
        !
        !    ! Read the family information
        !    call cg_goto_f(cg, base, ierr, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t",&
        !         boco, "end")
        !    if (ierr .eq. CG_ERROR) call cg_error_exit_f
        !
        !    if (ierr == 0) then ! Node exists
        !       call cg_famname_read_f(famName, ierr)
        !       if (ierr .eq. CG_ERROR) call cg_error_exit_f
        !
        !       zones(iZone)%bocos(boco)%family = famName
        !    end if
        ! end do bocoLoop

        ! And free the temporary arrays
        deallocate(coorX, coorY, coorZ)
        zoneStart = zoneStart + dims(1)

     end do zoneLoop

     call cg_close_f(cg, ierr)
     if (ierr .eq. CG_ERROR) call cg_error_exit_f
     deallocate(sizes)
  end if

  ! There is no MPI communication needed here because ADT only reads the
  ! data from the root processor.

end subroutine readUnstructuredCGNS

end module CGNSinterface
