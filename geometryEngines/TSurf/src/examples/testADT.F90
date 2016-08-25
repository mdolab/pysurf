program test

  ! EXTERNAL MODULES
  use precision
  use adtAPI

  implicit none

  ! DECLARE VARIABLES
  integer(kind=intType) :: ierr, nNodes, nTria, nQuads, nPts, nInterpol, i
  integer(kind=intType) :: comm, numProcs, myID
  integer(kind=intType) :: triaConn(3,3), quadsConn(4,1)
  integer(kind=intType), dimension(:), allocatable :: elementID, procID
  integer(kind=adtElementType), dimension(:), allocatable :: elementType
  real(kind=realType) :: coor(3,6), BBox(3,2), coorPts(3,1), arrDonor(1,1)
  real(kind=realType), dimension(:), allocatable :: dist2
  real(kind=realType), dimension(:,:), allocatable :: allxfs
  real(kind=realType), dimension(:,:), allocatable :: uvw, arrInterpol
  logical :: useBBox
  character(32) :: adtID

  !===============================
  ! USER INPUTS
  !===============================

  ! Initialize MPI
  call mpi_init(ierr)

  ! Set communicator
  comm = MPI_COMM_WORLD

  ! Define nodes
  coor = reshape([0.0, 0.0, 0.0, &
                  2.0, 0.0, 0.0, &
                  4.0, 0.0, 0.0, &
                  1.0, 1.0, 0.0, &
                  3.0, 1.0, 0.0, &
                  5.0, 1.0, 0.0], &
                  [3,6])

  ! Define traingle elements
  triaConn(:,1) = [1,2,4]
  triaConn(:,2) = [2,3,5]
  triaConn(:,3) = [2,5,4]

  ! Define quad elements
  quadsConn(:,1) = [1,2,5,4]
  ! quadsConn(:,2) = [2,3,6,5]


  ! Bounding box parameters
  BBox = reshape([0., 0., 0., 0., 0., 0.],[3,2])
  useBBox = .false.

  ! Set ADT id
  adtID = 'testSurf'

  ! Define points to be projected
  coorPts = reshape([3., .5, .555], &
                     [3,1])

  !===============================
  ! EXECUTION
  !===============================

  ! Determine number of processors and current processor ID
  call MPI_COMM_SIZE(comm, numProcs, ierr)
  call MPI_COMM_RANK(comm, myID, ierr)
  if (myID == 0) then
     print *,'Number of processors: ',numProcs
     print *,'coor'
     do i=1,3
        print *,coor(i,:)
     end do
     print *,'coorPts'
     do i=1,3
        print *,coorPts(i,:)
     end do
  else
    coor = 0.
  end if

  ! Determine number of nodes, elements, and projected points
  nNodes = size(coor, 2)
  nTria = size(triaConn, 2)
  nQuads = 0
  nPts = size(coorPts, 2)

  ! Setup ADT
  call adtBuildSurfaceADT(nTria, nQuads,   nNodes,    &
                          coor,  triaConn, quadsConn, &
                          BBox,  useBBox,  comm,      &
                          adtID)

  ! Allocate output arrays
  allocate(elementID(nPts), elementType(nPts), arrInterpol(1,nPts))
  allocate(procID(nPts), uvw(3, nPts), dist2(nPts))
  allocate(allxfs(3, nPts))


  ! No additional variable should be interpolated
  nInterpol = 0
  arrDonor = 0

  ! Initialize minimum distance so far to high numbers
  dist2 = 1.0e5

  ! Find minimum distance to the given point
  call adtMinDistanceSearch(nPts,      nNodes,      coorPts,     adtID,     &
                            procID,    elementType, elementID, &
                            uvw,       dist2,       allxfs,    &
                            nInterpol, arrDonor,    arrInterpol)

  if (myID == 0) then
     print *,'dist2'
     print *,dist2
     print *,'uvw'
     do i=1,3
         print *,uvw(i,:)
     end do
     print *,'proc ID'
     print *,procID
     print *,'elem ID'
     print *,elementID
     print *,'global coord'
     print *,allxfs
  end if

  ! Finalize MPI
  call mpi_finalize(ierr)

end program test
