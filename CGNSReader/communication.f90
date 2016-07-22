module communication

  use precision
  implicit none
  save

  ! myID:            My processor number in warp_comm_world.
  ! nProc:           The number of processors in warp_comm_world.

  integer(kind=intType) :: myID, nProc

end module communication
