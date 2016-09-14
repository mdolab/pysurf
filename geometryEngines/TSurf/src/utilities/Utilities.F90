module Utilities

use precision
implicit none

contains

!============================================================

subroutine condenseBarNodes_main(nNodes, nElem, distTol, &
                                 coor, barsConn, nUniqueNodes)

  ! This subroutine receives a list of bar FE which may have repeated points and
  ! Then condenses (merges) points that are closer than a given tolerance.
  ! This helps us getting a continuous line from the intersection algorithm results.
  ! Ney Secco 2016-08
  !
  ! INPUTS
  !
  ! distTol: real -> Distance tolerance used to determine if two nodes are effectively the same.
  !
  ! OUTPUTS
  !
  ! nUniqueNodes: integer -> Number of unique nodes.
  !
  ! INPUTS/OUTPUTS
  !
  ! coor: real(3,nNodes) -> Nodal coordinates. May have repeated coordinates at the beginning.
  !                         At the end, it will have zeros if the total number of nodes is reduced.
  !                         Then you should use nUniqueNodes to crop the matrix accordingly.
  !
  ! barsConn: int(2,nElements) -> Bar element connectivity array. Its values will be updated
  !           to use newCoor indices.


  implicit none

  ! INPUTS
  integer(kind=intType), intent(in) :: nNodes, nElem
  real(kind=realType), intent(in) :: distTol

  ! INPUTS/OUTPUTS
  real(kind=realType), dimension(3,nNodes), intent(inout) :: coor
  integer(kind=intType), dimension(2,nElem), intent(inout) :: barsConn

  ! OUTPUTS
  integer(kind=intType), intent(out) :: nUniqueNodes

  ! WORKING
  integer(kind=intType) :: nCopies
  integer(kind=intType) :: currNodeID, prevNodeID, link, elemID
  real(kind=realType), dimension(3) :: currCoor, prevCoor
  real(kind=realType) :: dist
  integer(kind=intType), dimension(size(coor,2)) :: linkOld2New

  ! EXECUTION

  ! Initialize number of unique nodes found so far.
  ! The first node is a unique one =P!
  nUniqueNodes = 1

  ! As the first node is unique, its old-to-new link should point to itself
  linkOld2New(1) = 1

  ! Loop over the nodes to find the unique ones
  do currNodeID = 2,nNodes

     ! Get coordinates of current node
     currCoor = coor(:,currNodeID)

     ! Now loop over the previous nodes to find if it is repeated
     do prevNodeID = 1,currNodeID-1

        ! Get coordinates of the previous node
        prevCoor = coor(:,prevNodeID)

        ! Compute distance between nodes
        call norm(currCoor - prevCoor, dist)

        ! Check if the distance is below the merging tolerance
        if (dist .le. distTol) then

           ! Update link array.
           ! The array linkOld2New will contain newCoor indices that correspond to each
           ! coor index.
           ! So the current node should use the same link as the previous node, as they will
           ! point to the same index of the new coordinate array.
           linkOld2New(currNodeID) = linkOld2New(prevNodeID)

           ! We can jump out of the prevNode do loop and check the next node.
           exit

        end if

     end do

     ! Check if we did not find any copy. In this case, we need to initialize a new
     ! unique node
     if (prevNodeID .eq. currNodeID) then

        ! Increase the number of unique nodes found so far
        nUniqueNodes = nUniqueNodes + 1

        ! Create new link
        linkOld2New(currNodeID) = nUniqueNodes

     end if

  end do

  ! Initialize number of nodes copied so far
  nCopies = 0

  ! We loop once again over the nodes so we can copy the unique values
  do currNodeID = 1,nNodes

     ! Get index of the current node in the new coordinate array
     link = linkOld2New(currNodeID)

     ! Check if the new link is already used
     if (link .gt. nCopies) then

        ! Get coordinates of current node
        currCoor = coor(:,currNodeID)

        ! Increment number of copies done so far
        nCopies = nCopies + 1

        ! Copy coordinates
        coor(:,nCopies) = currCoor

     end if

  end do

  ! Now zero out the unused point in coor
  coor(:,nCopies+1:nNodes) = 0.0

  ! Now the last step is updating the bars connectivity.
  ! Loop over the elements
  do elemID = 1,nElem

     ! Update connectivities
     barsConn(1,elemID) = linkOld2New(barsConn(1,elemID))
     barsConn(2,elemID) = linkOld2New(barsConn(2,elemID))

  end do

end subroutine condenseBarNodes_main


subroutine remesh_main(nNodes, nNewNodes, coor, barsConn, method, spacing, newCoor, newBarsConn)

  ! This function will redistribute the nodes along the curve to get
  ! better spacing among the nodes.
  ! This assumes that the FE data is ordered.
  !
  ! The first and last nodes will remain at the same place.
  !
  ! Consider using self.shift_end_nodes first so you can preserve the start and end points
  ! of a periodic curve.
  !
  ! INPUTS:
  !
  ! nNodes: integer -> Number of input nodes
  !
  ! nNewNodes: integer -> Number of new nodes desired in the new curve definition.
  !
  ! coor: real, shape(3, nNodes) -> Coordinates of the input curve
  !
  ! barsConn: real, shape(2, nNodes-1) -> Connectivity of the input bar elements
  !
  ! method: string -> Method used to interpolate new nodes with respect to
  !         existing ones. Check scipy.interpolate.interp1d for options.
  !
  ! spacing: string -> Desired spacing criteria for new nodes. Current options are:
  !          ['linear']
  !
  !
  ! OUTPUTS:
  !
  ! newCoor: real, shape(3, nNodes) -> Coordinates of the output curve
  !
  ! newBarsConn: real, shape(2, nNodes-1) -> Connectivity of the output bar elements
  !
  ! John Jasa, adapted from Ney Secco's Python 2016-09

  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: nNodes, nNewNodes
  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  integer(kind=intType), dimension(2,nNodes-1), intent(in) :: barsConn

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes), intent(out) :: newCoor
  integer(kind=intType), dimension(2,nNewNodes-1), intent(out) :: newBarsConn

  ! Working variables
  real(kind=realType), dimension(3,nNodes) :: nodeCoor
  real(kind=realType), dimension(nNodes) :: arcLength
  integer(kind=intType) :: elemID, prevNodeID, currNodeID, nElem
  real(kind=realType) :: dist

  real(kind=realType), dimension(nNewNodes) :: newArcLength
  real(kind=realType), dimension(3) :: node1, node2
  logical :: periodic

  nElem = nNodes - 1

  ! First we check if the FE data is ordered
  do elemID=3,nElem

    ! Get node indices
    prevNodeID = barsConn(2, elemID-1)
    currNodeID = barsConn(1, elemID)

    ! Check if the FE data is ordered
    if (prevNodeID .ne. currNodeID) then

      ! Print warning
      print *, 'WARNING: Could not remesh curve because it has unordered FE data.'
      print *, '         Call FEsort first.'
      return
    end if
  end do

  ! COMPUTE ARC-LENGTH

  ! We can proceed if FE data is ordered

  ! Store position of the first node (the other nodes will be covered in the loop)
  ! (the -1 is due Fortran indexing)
  nodeCoor(:,1) = coor(:,barsConn(1,1))

  ! Loop over each element to increment arcLength
  do elemID=1,nElem

    ! Get node positions (the -1 is due Fortran indexing)
    node1 = coor(:,barsConn(1,elemID))
    node2 = coor(:,barsConn(2,elemID))

    ! Compute distance between nodes
    call norm(node1 - node2, dist)

    ! Store nodal arc-length
    arcLength(elemID+1) = arcLength(elemID) + dist

    ! Store coordinates of the next node
    nodeCoor(:,elemID+1) = node2

  end do

  ! SAMPLING POSITION FOR NEW NODES

  ! Now that we know the initial and final arcLength, we can redistribute the
  ! parametric coordinates based on the used defined spacing criteria.
  ! These statements should initially create parametric coordinates in the interval
  ! [0.0, 1.0]. We will rescale it after the if statements.
  if (spacing .eq. 'linear') then
    call linspace(0.0, 1.0, nNewNodes, newArcLength)

  else if (spacing .eq. 'cosine') then
    call linspace(0.0, 3.141592653589793, nNewNodes, newArcLength)
    newArcLength = 0.5 * (1.0 - cos(newArcLength))
  end if

  ! Rescale newArcLength based on the final distance
  newArcLength = arcLength(nNodes) * newArcLength

  ! INTERPOLATE NEW NODES

  ! Now we sample the new coordinates based on the interpolation method given by the user

  ! Create interpolants for x, y, and z
  call interp1d(1, nNodes, arcLength, nodeCoor(1,:), nNewNodes, newArcLength, newCoor(1,:))
  call interp1d(1, nNodes, arcLength, nodeCoor(2,:), nNewNodes, newArcLength, newCoor(2,:))
  call interp1d(1, nNodes, arcLength, nodeCoor(3,:), nNewNodes, newArcLength, newCoor(3,:))

  ! ASSIGN NEW COORDINATES AND CONNECTIVITIES

  ! Check if the baseline curve is periodic
  if (barsConn(1, 1) .eq. barsConn(2, nNodes-1)) then
    periodic = .true.
  else
    periodic = .false.
  end if

  ! Generate new connectivity (the nodes are already in order so we just
  ! need to assign an ordered set to barsConn).
  do elemID=1,nNewNodes-1
    newBarsConn(1, elemID) = elemID
    newBarsConn(2, elemID) = elemID + 1
  end do

  ! We still need to keep periodicity if the original curve is periodic
  if (periodic) then
    newBarsConn(2, nNewNodes-1) = newBarsConn(1, 1)
  end if

end subroutine

!============================================================

subroutine norm(v, norm_)

  implicit none

  real(kind=realType), intent(in) :: v(3)
  real(kind=realType), intent(out) :: norm_
  real(kind=realType) :: dot_prod


  !norm = sqrt(dot_product(v, v))
  call dot(v, v, dot_prod)
  norm_ = dot_prod ** 0.5

end subroutine norm

subroutine dot(a, b, dot_prod)

  implicit none

  real(kind=realType), intent(in) :: a(3), b(3)
  real(kind=realType), intent(out) :: dot_prod

  dot_prod = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

end subroutine dot

subroutine linspace(l, k, n, z)

  implicit none

  !// Argument declarations
  integer(kind=intType), intent(in) :: n
  real(kind=realType), dimension(n), intent(out) :: z
  real(kind=realType), intent(in) :: l
  real(kind=realType), intent(in) :: k

  !// local variables
  integer(kind=intType) :: i
  real(kind=realType) :: d

  d = (k - l) / (n - 1)
  z(1) = l
  do i = 2, n-1
    z(i) = z(i-1) + d
  end do
  z(1) = l
  z(n) = k
  return

end subroutine linspace

subroutine interp1d ( m, data_num, t_data, p_data, interp_num, &
  t_interp, p_interp )

!*****************************************************************************80
!
!! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
!
!  Discussion:
!
!    From a space of M dimensions, we are given a sequence of
!    DATA_NUM points, which are presumed to be successive samples
!    from a curve of points P.
!
!    We are also given a parameterization of this data, that is,
!    an associated sequence of DATA_NUM values of a variable T.
!    The values of T are assumed to be strictly increasing.
!
!    Thus, we have a sequence of values P(T), where T is a scalar,
!    and each value of P is of dimension M.
!
!    We are then given INTERP_NUM values of T, for which values P
!    are to be produced, by linear interpolation of the data we are given.
!
!    Note that the user may request extrapolation.  This occurs whenever
!    a T_INTERP value is less than the minimum T_DATA or greater than the
!    maximum T_DATA.  In that case, linear extrapolation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
!    independent variable at the sample points.  The values of T_DATA
!    must be strictly increasing.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
!    dependent variables at the sample points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
!    at which interpolation is to be done.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
!    independent variable at the interpolation points.
!
!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
!    values of the dependent variables at the interpolation points.
!
  implicit none

  integer ( kind = intType ) data_num
  integer ( kind = intType ) m
  integer ( kind = intType ) interp_num

  integer ( kind = intType ) interp
  integer ( kind = intType ) left
  real ( kind = realType ) p_data(data_num)
  real ( kind = realType ) p_interp(interp_num)
  integer ( kind = intType ) right
  real ( kind = realType ) t
  real ( kind = realType ) t_data(data_num)
  real ( kind = realType ) t_interp(interp_num)

  do interp = 1, interp_num

    t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
    call r8vec_bracket ( data_num, t_data, t, left, right )

    p_interp(interp) = &
      ( ( t_data(right) - t                ) * p_data(left)   &
      + (                 t - t_data(left) ) * p_data(right) ) &
      / ( t_data(right)     - t_data(left) )

  end do

  return
end

subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = intType ) n

  integer ( kind = intType ) i
  integer ( kind = intType ) left
  integer ( kind = intType ) right
  real ( kind = realType ) x(n)
  real ( kind = realType ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end

end module Utilities
