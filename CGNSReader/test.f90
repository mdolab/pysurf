program test

  ! EXTERNAL MODULES
  use precision
  use communication
  use readCGNS

  implicit none

  ! DECLARE VARIABLES

  !===============================
  ! USER INPUTS
  !===============================
  integer(kind=intType) :: comm, ierr

  integer(kind=intType), dimension(:,:), allocatable :: triaConn, quadsConn
  real(kind=realType), dimension(:,:), pointer :: coor

  ! Initialize MPI
  call mpi_init(ierr)

  ! Set communicator
  comm = MPI_COMM_WORLD

  !===============================
  ! EXECUTION
  !===============================

  call readCGNS_routine('cube2.cgns', coor, triaConn, quadsConn)

  ! Determine number of processors and current processor ID
  call MPI_COMM_SIZE(comm, nProc, ierr)
  call MPI_COMM_RANK(comm, myID, ierr)
  if (myID == 0) then
    print *,coor
    print *,''
    print *,triaConn
    print *,quadsConn
  end if

  ! Finalize MPI
  call mpi_finalize(ierr)

end program test
