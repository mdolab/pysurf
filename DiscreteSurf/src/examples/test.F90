program test

  ! EXTERNAL MODULES
  use precision
  use communication
  use CGNSinterface

  implicit none

  ! DECLARE VARIABLES

  !===============================
  ! USER INPUTS
  !===============================
  integer(kind=intType) :: comm, ierr, ii

  integer(kind=intType), dimension(:,:), allocatable :: triaConn, quadsConn, barsConn
  integer(kind=intType), dimension(:), allocatable :: surfTriaPtr, surfQuadsPtr, curveBarsPtr
  real(kind=realType), dimension(:,:), allocatable :: coor

  ! Initialize MPI
  call mpi_init(ierr)

  ! Set communicator
  comm = MPI_COMM_WORLD

  !===============================
  ! EXECUTION
  !===============================

  call readCGNSmain('cube2.cgns', comm, coor, triaConn, quadsConn, barsConn, &
                    surfTriaPtr, surfQuadsPtr, curveBarsPtr)

  print *,surfTriaPtr

  ! Determine number of processors and current processor ID
  call MPI_COMM_SIZE(comm, nProc, ierr)
  call MPI_COMM_RANK(comm, myID, ierr)
  if (myID == 0) then
     do ii=1,size(coor,2)
        print *,coor(:,ii)
     end do
    print *,''
    print *,triaConn
    print *,quadsConn
  end if

  ! Finalize MPI
  call mpi_finalize(ierr)

end program test
