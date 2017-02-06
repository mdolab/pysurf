module Utilities

use precision
implicit none

contains

!============================================================

subroutine condenseBarNodes_main(nNodes, nElem, distTol, &
                                 coor, barsConn, nUniqueNodes, linkOld2New)

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
  ! linkOld2New: integer(nNodes) -> This array maps the old coordinates (input coor) to the new
  !                                 set of coordinates (output coor).
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
  integer(kind=intType), dimension(nNodes), intent(out) :: linkOld2New

  ! WORKING
  integer(kind=intType) :: nCopies
  integer(kind=intType) :: currNodeID, prevNodeID, link, elemID
  real(kind=realType), dimension(3) :: currCoor, prevCoor
  real(kind=realType) :: dist
  integer(kind=intType), dimension(:), allocatable :: numAddedNodes
  real(kind=realType), dimension(:,:), allocatable :: newCoor

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

  ! Allocate an array that stores how many nodes will be merged to create a new node
  allocate(numAddedNodes(nUniqueNodes))
  numAddedNodes = 0

  ! Allocate array to store merged nodes
  allocate(newCoor(3,nUniqueNodes))
  newCoor = 0.0

  ! Initialize number of nodes copied so far
  nCopies = 0

  ! We loop once again over the nodes so we can add the coordinates of nodes to be merged.
  ! We will take the average later on.
  do currNodeID = 1,nNodes

     ! Get index of the current node in the new coordinate array
     link = linkOld2New(currNodeID)

     ! Add the current node to the corresponding location on the new curve
     newCoor(:,link) = newCoor(:,link) + coor(:,currNodeID)

     ! Increment number of nodes merged to this new location, so we can
     ! take the average later on
     numAddedNodes(link) = numAddedNodes(link) + 1

  end do

  ! Now reset the given set of coordinates, so we can add just the new nodes
  coor(:,nCopies+1:nNodes) = 0.0

  ! Take the average and copy the new nodes to the original coordinate array
  do currNodeID = 1,nUniqueNodes

     coor(:,currNodeID) = newCoor(:,currNodeID)/numAddedNodes(currNodeID)

  end do

  ! Now the last step is updating the bars connectivity.
  ! Loop over the elements
  do elemID = 1,nElem

     ! Update connectivities
     barsConn(1,elemID) = linkOld2New(barsConn(1,elemID))
     barsConn(2,elemID) = linkOld2New(barsConn(2,elemID))

  end do

end subroutine condenseBarNodes_main


subroutine remesh_main(nNodes, nElem, nNewNodes, coor, barsConn, method,&
  spacing, Sp1, Sp2, newCoor, newBarsConn)

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
  ! Sp1: real -> Desired initial spacing for tangent/hypTan remeshes.
  !
  ! Sp2: real -> Desired final spacing for tangent/hypTan remeshes.
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
  integer(kind=intType), intent(in) :: nNewNodes, nElem
  integer(kind=intType), intent(inout) :: nNodes
  character(32), intent(in) :: method, spacing
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  integer(kind=intType), dimension(2,nElem), intent(in) :: barsConn
  real(kind=realType), intent(in) :: Sp1, Sp2

  ! Output variables
  real(kind=realType), dimension(3,nNewNodes) :: newCoor
  integer(kind=intType), dimension(2,nNewNodes-1) :: newBarsConn

  ! Working variables
  real(kind=realType), dimension(3,nNodes) :: nodeCoor
  real(kind=realType), dimension(nNodes) :: arcLength
  integer(kind=intType) :: elemID, prevNodeID, currNodeID
  real(kind=realType) :: dist, zero, one, pi

  real(kind=realType), dimension(nNewNodes) :: newArcLength
  real(kind=realType), dimension(3) :: node1, node2

  nNodes = nElem + 1

  ! Initialize outputs
  newCoor = 0.0
  newBarsConn = 0

  ! First we check if the FE data is ordered
  do elemID=2,nElem

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
  arcLength(1) = 0.0

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

  zero = 0.
  one = 1.
  pi = 3.1415926535897932384626

  ! Now that we know the initial and final arcLength, we can redistribute the
  ! parametric coordinates based on the used defined spacing criteria.
  ! These statements should initially create parametric coordinates in the interval
  ! [0.0, 1.0]. We will rescale it after the if statements.
  if (spacing .eq. 'linear') then
    call linspace(zero, one, nNewNodes, newArcLength)
  else if (spacing .eq. 'cosine') then
    call linspace(zero, pi, nNewNodes, newArcLength)
    newArcLength = 0.5 * (1.0 - cos(newArcLength))
  else if (spacing .eq. 'hyptan') then
    call getHypTanDist(Sp1/arcLength(nElem+1), Sp2/arcLength(nElem+1), nNewNodes, newArcLength)
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

  ! Generate new connectivity (the nodes are already in order so we just
  ! need to assign an ordered set to barsConn).
  do elemID=1,nNewNodes-1
    newBarsConn(1, elemID) = elemID
    newBarsConn(2, elemID) = elemID + 1
  end do

end subroutine

!============================================================

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

!============================================================

subroutine dot(A, B, dot_)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(in) :: A(3), B(3)
  real(kind=realType), intent(out) :: dot_

  dot_ = sum(A*B)

end subroutine dot

!============================================================

subroutine norm(A, norm_)

  ! John Jasa - 2016-08

  implicit none

  real(kind=realType), intent(in) :: A(3)
  real(kind=realType), intent(out) :: norm_

  norm_ = sqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))

end subroutine norm

subroutine getHypTanDist(Sp1, Sp2, N, spacings)

  implicit none

  real(kind=realType), intent(in) :: Sp1, Sp2
  integer(kind=intType), intent(in) :: N

  real(kind=realType), intent(out), dimension(N) :: spacings

  integer(kind=intType) :: i
  real(kind=realType) :: b, step, b_out, b_out_step, f_prime, b_old
  real(kind=realType) :: A, U, R

  ! Manually created secant method for the above solve
  b = 4.
  step = 1.e-6
  do i=1,1000
    call findRootb(b, Sp1, Sp2, N, b_out)
    call findRootb(b-step, Sp1, Sp2, N, b_out_step)

    b_old = b
    f_prime = (b_out - b_out_step) / step
    b = b - b_out / f_prime

    if (abs(b_old - b) .lt. 1.e-10) then
      exit
    end if
  end do

  ! Compute parameter A
  A = sqrt(Sp1/Sp2)

  do i = 1,N
    R = float(i-1) / float(N-1) - .5
    U = 1 + tanh(b*R) / tanh(b/2)
    spacings(i) = U/(2*A + (1-A)*U)
  end do

end subroutine

subroutine findRootb(b, Sp1, Sp2, N, b_out)

  implicit none

  real(kind=realType), intent(in) :: b, Sp1, Sp2
  integer(kind=intType), intent(in) :: N

  real(kind=realType), intent(out) :: b_out

  b_out = sinh(b) - b/(N-1)/sqrt(Sp1*Sp2)

end subroutine


!============================================================
! BOUNDING BOX ROUTINES
!============================================================

subroutine computeBBox(coor, BBox)

  ! This subroutine computes the bounding box of all points given in coor.
  ! Ney Secco - 2016-08
  !
  ! INPUTS:
  !
  ! coor: real(3, nNodes) -> coordinates (x,y,z) of all nodes used in the
  !       triangulated surface definition.
  !
  ! OUTPUTS:
  !
  ! BBox: real(3, 2) -> coordinates defining the bounding box
  !                     [xmin, xmax]
  !                     [ymin, ymax]
  !                     [zmin, zmax]
  ! Ney Secco 2016-08

  implicit none

  ! INPUTS
  real(kind=realType), dimension(:,:), intent(in) :: coor

  ! OUTPUTS
  real(kind=realType), dimension(3,2), intent(out) :: BBox

  ! EXECUTION

  ! Get bounding values
  BBox(:, 1) = minval(coor, 2)
  BBox(:, 2) = maxval(coor, 2)

end subroutine computeBBox

!============================================================

subroutine computeBBoxPerElements(nNodes, nTria, nQuads, &
                                  coor, triaConn, quadsConn, &
                                  triaBBox, quadsBBox)

  ! This subroutine computes a bounding box for each element in triaConn and quadsConn.
  !
  ! INPUTS
  !
  ! coor: real(3,nNodes) -> Nodal coordinates
  !
  ! triaConn: integer(3,nTria) -> Connectivities of triangle elements
  !
  ! quadsConn: integer(4,nQuads) -> Connectivities of quad elements
  !
  ! OUTPUTS:
  !
  ! triaBBox: real(6,nTria) -> Bounding box coordinates for each triangle element.
  !                            (xmin, ymin, zmin, xmax, ymax, zmax).
  !
  ! quadsBBox: real(6,nQuads) -> Bounding box coordinates for each quad element.
  !                              (xmin, ymin, zmin, xmax, ymax, zmax).
  ! Ney Secco 2016-10

  implicit none

  ! INPUTS
  integer(kind=intType), intent(in) :: nNodes, nTria, nQuads
  real(kind=realType), dimension(3,nNodes), intent(in) :: coor
  integer(kind=intType), dimension(3,nTria), intent(in) :: triaConn
  integer(kind=intType), dimension(4,nQuads), intent(in) :: quadsConn

  ! OUTPUTS
  real(kind=realType), dimension(6,nTria), intent(out) :: triaBBox
  real(kind=realType), dimension(6,nQuads), intent(out) :: quadsBBox

  ! WORKING
  integer(kind=intType) :: elemID, nodeID
  real(kind=realType) :: triaCoor(3,3), quadsCoor(3,4), BBox(6)

  ! EXECUTION

  ! Loop over all triangles
  do elemID = 1,nTria

     ! Get coordinates of each node and store them in a single matrix
     triaCoor(:,1) = coor(:,triaConn(1,elemID))
     triaCoor(:,2) = coor(:,triaConn(2,elemID))
     triaCoor(:,3) = coor(:,triaConn(3,elemID))

     ! Assign min values (BBox(1:3)) and max values (BBox(4:6))
     ! based on the nodal coordinates.
     triaBBox(1:3, elemID) = minval(real(triaCoor), 2)
     triaBBox(4:6, elemID) = maxval(real(triaCoor), 2)

  end do

  ! Loop over all quads
  do elemID = 1,nQuads

     ! Get coordinates of each node and store them in a single matrix
     quadsCoor(:,1) = coor(:,quadsConn(1,elemID))
     quadsCoor(:,2) = coor(:,quadsConn(2,elemID))
     quadsCoor(:,3) = coor(:,quadsConn(3,elemID))
     quadsCoor(:,4) = coor(:,quadsConn(4,elemID))

     ! Assign min values (BBox(1:3)) and max values (BBox(4:6))
     ! based on the nodal coordinates.
     quadsBBox(1:3, elemID) = minval(real(quadsCoor), 2)
     quadsBBox(4:6, elemID) = maxval(real(quadsCoor), 2)

  end do

end subroutine computeBBoxPerElements

!============================================================

subroutine computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

  ! This subroutine gives a rectangle (BBoxAB) that correponds to the intersection
  ! of two bounding boxes (BBoxA, and BBoxB).
  ! If there is no intersection, we will return overlap = .false., but BBoxAB will
  ! still have meaningless value.
  ! Ney Secco - 2016-08
  !
  ! INPUTS
  !
  ! BBoxA: real(6) -> bounding box A coordinates (xmin, ymin, zmin, xmax, ymax, zmax).
  !
  ! BBoxB: real(6) -> bounding box B coordinates (xmin, ymin, zmin, xmax, ymax, zmax).
  !
  ! OUTPUTS
  !
  ! BBoxAB: real(6) -> coordinates of the A-B intersection bounding box.
  !         ATTENTION: This variable may still have meaningless values if
  !         overlap = .false.
  !
  ! overlap: logical -> logical indicating if there is an intersection
  !          between boxes A and B. If overlap = .false., disregard all
  !          values given in BBoxAB.

  implicit none

  ! INPUTS
  real(kind=realType), dimension(6), intent(in) :: BBoxA, BBoxB

  ! OUTPUTS
  real(kind=realType), dimension(6), intent(out) :: BBoxAB
  logical, intent(out) :: overlap

  ! WORKING
  real(kind=realType), dimension(6) :: bounds

  ! EXECUTION

  ! Check overlaps along each dimension
  ! We have an actual BBox intersection if there are overlaps in all dimensions

  ! X overlap
  call lineIntersectionInterval(BBoxA(1), BBoxA(4), BBoxB(1), BBoxB(4), &
       BBoxAB(1), BBoxAB(4), overlap)

  if (overlap) then
     ! Y overlap
     call lineIntersectionInterval(BBoxA(2), BBoxA(5), BBoxB(2), BBoxB(5), &
          BBoxAB(2), BBoxAB(5), overlap)

     if (overlap) then
        ! Z overlap
        call lineIntersectionInterval(BBoxA(3), BBoxA(6), BBoxB(3), BBoxB(6), &
             BBoxAB(3), BBoxAB(6), overlap)

     end if

  end if

  ! Determine the size of the edges of the bounding box
  bounds(1) = BBoxAB(4) - BBoxAB(1)
  bounds(2) = BBoxAB(5) - BBoxAB(2)
  bounds(3) = BBoxAB(6) - BBoxAB(3)
  bounds(4:6) = -bounds(1:3)

  ! Add buffer space to the bounding box.
  ! Although this slightly slows down the intersection algorithm, it is on the
  ! order of hundredths of seconds for the CRM case and solves a small bug.
  ! Will need to fix the filterElements code later so we can use a tighter
  ! bounding box.
  BBoxAB = BBoxAB - 0.01 * bounds

  ! Remember that overlap may be set to false in the Z overlap test

end subroutine computeBBoxIntersection

!============================================================

subroutine lineIntersectionInterval(xminA, xmaxA, xminB, xmaxB, xminAB, xmaxAB, overlap)

  ! This subroutine finds a 1D overlap interval between two lines.
  ! This step is used three times (one for each dimension) when
  ! computing bounding boxes intersection.
  ! If the bounding boxes do not overlap we return overlap = .false.
  ! Ney Secco - 2016-08

  implicit none

  ! INPUTS
  real(kind=realType), intent(in) :: xminA, xmaxA, xminB, xmaxB

  ! OUTPUTS
  real(kind=realType), intent(out) :: xminAB, xmaxAB
  logical, intent(out) :: overlap

  ! EXECUTION

  ! Get bounding interval
  xminAB = max(xminA, xminB)
  xmaxAB = min(xmaxA, xmaxB)

  ! Check if we actually have an overlap
  if (xminAB .gt. xmaxAB) then
     overlap = .false.
  else
     overlap = .true.
  end if

end subroutine lineIntersectionInterval

end module Utilities
