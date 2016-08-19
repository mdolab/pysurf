module communication

  use precision
  implicit none
  save

  ! comm_world: The communicator of this processor group.
  ! comm_self : The single processor communicator 
  ! myID:       My processor number in hyp_comm_world.
  ! nProc:      The number of processors in hyp_comm_world.

  integer :: comm_world, comm_self, myID, nProc

end module communication
