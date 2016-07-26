program test

! This test uses CGNSinterface to get data and then run ADT on it

  ! EXTERNAL MODULES
  use precision
  use communication
  use CGNSinterface
  use adtAPI

  implicit none

  ! DECLARE VARIABLES

  integer(kind=intType) :: comm, ierr, ii
  integer(kind=intType), dimension(:,:), allocatable :: triaConn, quadsConn
  real(kind=realType), dimension(:,:), allocatable :: coor
  character(32) :: fileName

  integer(kind=intType) :: nNodes, nTria, nQuads, nPts, nInterpol
  integer(kind=intType), dimension(:), allocatable :: elementID, procID
  integer(kind=adtElementType), dimension(:), allocatable :: elementType
  real(kind=realType) :: BBox(3,2), coorPts(3,1), arrDonor(1,1)
  real(kind=realType), dimension(:), allocatable :: dist2
  real(kind=realType), dimension(:,:), allocatable :: uvw, arrInterpol
  real(kind=realType), dimension(:,:), allocatable :: allNorms, allxfs
  logical :: useBBox
  character(32) :: adtID

  !===============================
  ! USER INPUTS
  !===============================

  ! CGNS file input
  fileName = 'cube2.cgns'

  ! Bounding box parameters
  BBox = reshape([0., 0., 0., 0., 0., 0.],[3,2])
  useBBox = .false.

  ! Set ADT id
  adtID = 'testSurf'

  ! Define points to be projected
  coorPts = reshape([0.5, 0.5, 1.0], &
                     [3,1])

  !===============================
  ! EXECUTION
  !===============================

  ! Initialize MPI
  call mpi_init(ierr)

  ! Set communicator
  comm = MPI_COMM_WORLD

  ! Determine number of processors and current processor ID
  call MPI_COMM_SIZE(comm, nProc, ierr)
  call MPI_COMM_RANK(comm, myID, ierr)

  ! Read CGNS file
  call readCGNS(fileName, comm, coor, triaConn, quadsConn)

  ! Determine number of nodes, elements, and projected points
  nNodes = size(coor, 2)
  nTria = size(triaConn, 2)
  nQuads = size(quadsConn, 2)
  nPts = size(coorPts, 2)

  ! Setup ADT
  call adtBuildSurfaceADT(nTria, nQuads,   nNodes,    &
                          coor,  triaConn, quadsConn, &
                          BBox,  useBBox,  comm,      &
                          adtID)

  ! Allocate output arrays
  allocate(elementID(nPts), elementType(nPts), arrInterpol(1,nPts))
  allocate(procID(nPts), uvw(3, nPts), dist2(nPts))
  allocate(allxfs(3, nPts), allNorms(3, nPts))

  ! No additional variable should be interpolated
  nInterpol = 0
  arrDonor = 0

  ! Initialize minimum distance so far to high numbers
  dist2 = 1.0e5

  ! Find minimum distance to the given point
  call adtMinDistanceSearch(nPts,     nNodes,      coorPts,     adtID,     &
                            procID,   elementType, elementID, &
                            uvw,      dist2,       allxfs,    &
                            allNorms, nInterpol, &
                            arrDonor, arrInterpol)

  if (myID == 0) then
     print *,'dist2'
     print *,dist2
     print *,'uvw'
     do ii=1,3
         print *,uvw(ii,:)
     end do
     print *,'proc ID'
     print *,procID
     print *,'elem ID'
     print *,elementID
     print *,'global coord'
     print *,allxfs
     print *,'normal'
     print *,allNorms
  end if

  ! Finalize MPI
  call mpi_finalize(ierr)

end program test
