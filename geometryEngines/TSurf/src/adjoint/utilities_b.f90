!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
!
MODULE UTILITIES_B
  USE PRECISION
  IMPLICIT NONE

CONTAINS
!============================================================
  SUBROUTINE CONDENSEBARNODES_MAIN(nnodes, nelem, disttol, coor, &
&   barsconn, nuniquenodes, linkold2new)
    IMPLICIT NONE
! INPUTS
    INTEGER(kind=inttype), INTENT(IN) :: nnodes, nelem
    REAL(kind=realtype), INTENT(IN) :: disttol
! INPUTS/OUTPUTS
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(INOUT) :: coor
    INTEGER(kind=inttype), DIMENSION(2, nelem), INTENT(INOUT) :: &
&   barsconn
! OUTPUTS
    INTEGER(kind=inttype), INTENT(OUT) :: nuniquenodes
    INTEGER(kind=inttype), DIMENSION(nnodes), INTENT(OUT) :: linkold2new
! WORKING
    INTEGER(kind=inttype) :: ncopies
    INTEGER(kind=inttype) :: currnodeid, prevnodeid, link, elemid
    REAL(kind=realtype), DIMENSION(3) :: currcoor, prevcoor
    REAL(kind=realtype) :: dist
    INTEGER(kind=inttype), DIMENSION(:), ALLOCATABLE :: numaddednodes
    REAL(kind=realtype), DIMENSION(:, :), ALLOCATABLE :: newcoor
    REAL(kind=realtype), DIMENSION(3) :: arg1
! EXECUTION
! Initialize number of unique nodes found so far.
! The first node is a unique one =P!
    nuniquenodes = 1
! As the first node is unique, its old-to-new link should point to itself
    linkold2new(1) = 1
! Loop over the nodes to find the unique ones
    DO currnodeid=2,nnodes
! Get coordinates of current node
      currcoor = coor(:, currnodeid)
! Now loop over the previous nodes to find if it is repeated
      DO prevnodeid=1,currnodeid-1
! Get coordinates of the previous node
        prevcoor = coor(:, prevnodeid)
! Compute distance between nodes
        arg1(:) = currcoor - prevcoor
        CALL NORM(arg1(:), dist)
! Check if the distance is below the merging tolerance
        IF (dist .LE. disttol) THEN
! Update link array.
! The array linkOld2New will contain newCoor indices that correspond to each
! coor index.
! So the current node should use the same link as the previous node, as they will
! point to the same index of the new coordinate array.
          linkold2new(currnodeid) = linkold2new(prevnodeid)
          GOTO 100
        END IF
      END DO
! We can jump out of the prevNode do loop and check the next node.
! Check if we did not find any copy. In this case, we need to initialize a new
! unique node
 100  IF (prevnodeid .EQ. currnodeid) THEN
! Increase the number of unique nodes found so far
        nuniquenodes = nuniquenodes + 1
! Create new link
        linkold2new(currnodeid) = nuniquenodes
      END IF
    END DO
! Allocate an array that stores how many nodes will be merged to create a new node
    ALLOCATE(numaddednodes(nuniquenodes))
    numaddednodes = 0
! Allocate array to store merged nodes
    ALLOCATE(newcoor(3, nuniquenodes))
    newcoor = 0.0
! Initialize number of nodes copied so far
    ncopies = 0
! We loop once again over the nodes so we can add the coordinates of nodes to be merged.
! We will take the average later on.
    DO currnodeid=1,nnodes
! Get index of the current node in the new coordinate array
      link = linkold2new(currnodeid)
! Add the current node to the corresponding location on the new curve
      newcoor(:, link) = newcoor(:, link) + coor(:, currnodeid)
! Increment number of nodes merged to this new location, so we can
! take the average later on
      numaddednodes(link) = numaddednodes(link) + 1
    END DO
! Now reset the given set of coordinates, so we can add just the new nodes
    coor(:, ncopies+1:nnodes) = 0.0
! Take the average and copy the new nodes to the original coordinate array
    DO currnodeid=1,nuniquenodes
      coor(:, currnodeid) = newcoor(:, currnodeid)/numaddednodes(&
&       currnodeid)
    END DO
! Now the last step is updating the bars connectivity.
! Loop over the elements
    DO elemid=1,nelem
! Update connectivities
      barsconn(1, elemid) = linkold2new(barsconn(1, elemid))
      barsconn(2, elemid) = linkold2new(barsconn(2, elemid))
    END DO
  END SUBROUTINE CONDENSEBARNODES_MAIN
!  Differentiation of remesh_main in reverse (adjoint) mode:
!   gradient     of useful results: newcoor
!   with respect to varying inputs: coor newcoor
!   RW status of diff variables: coor:out newcoor:in-zero
  SUBROUTINE REMESH_MAIN_B(nnodes, nelem, nnewnodes, coor, coorb, &
&   barsconn, method, spacing, sp1, sp2, newcoor, newcoorb, newbarsconn)
    IMPLICIT NONE
! Input variables
    INTEGER(kind=inttype), INTENT(IN) :: nnewnodes, nelem
    INTEGER(kind=inttype), INTENT(INOUT) :: nnodes
    CHARACTER(len=32), INTENT(IN) :: method, spacing
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(IN) :: coor
    REAL(kind=realtype), DIMENSION(3, nnodes) :: coorb
    INTEGER(kind=inttype), DIMENSION(2, nelem), INTENT(IN) :: barsconn
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
! Output variables
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoor
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoorb
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes - 1) :: newbarsconn
! Working variables
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoor
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoorb
    REAL(kind=realtype), DIMENSION(nnodes) :: arclength
    REAL(kind=realtype), DIMENSION(nnodes) :: arclengthb
    INTEGER(kind=inttype) :: elemid, prevnodeid, currnodeid, ii, jj
    REAL(kind=realtype) :: dist, zero, one, pi, disttol
    REAL(kind=realtype) :: distb
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclength
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclengthb
    REAL(kind=realtype), DIMENSION(3) :: node1, node2, newnode, oldnode
    REAL(kind=realtype), DIMENSION(3) :: node1b, node2b
    REAL(kind=realtype), DIMENSION(3) :: distvec
    INTRINSIC COS
    INTRINSIC SQRT
    REAL(kind=realtype), DIMENSION(3) :: arg1
    REAL(kind=realtype), DIMENSION(3) :: arg1b
    REAL(kind=realtype) :: arg10
    REAL(kind=realtype) :: arg10b
    REAL(kind=realtype) :: arg2
    REAL(kind=realtype) :: arg2b
    INTEGER :: ad_count
    INTEGER :: i
    INTEGER :: branch
! Tolerance to avoid any interpolation if new node is too close to an original node
    disttol = 1e-7
    nnodes = nelem + 1
! Initialize outputs
    newcoor = 0.0
    ad_count = 1
! First we check if the FE data is ordered
    DO elemid=2,nelem
! Get node indices
      prevnodeid = barsconn(2, elemid-1)
      currnodeid = barsconn(1, elemid)
! Check if the FE data is ordered
      IF (prevnodeid .NE. currnodeid) THEN
        GOTO 100
      ELSE
        ad_count = ad_count + 1
      END IF
    END DO
    CALL PUSHCONTROL1B(0)
    CALL PUSHINTEGER4(ad_count)
! COMPUTE ARC-LENGTH
! We can proceed if FE data is ordered
! Store position of the first node (the other nodes will be covered in the loop)
! (the -1 is due Fortran indexing)
    nodecoor(:, 1) = coor(:, barsconn(1, 1))
    arclength(1) = 0.0
! Loop over each element to increment arcLength
    DO elemid=1,nelem
! Get node positions (the -1 is due Fortran indexing)
      node1 = coor(:, barsconn(1, elemid))
      node2 = coor(:, barsconn(2, elemid))
! Compute distance between nodes
      arg1(:) = node1 - node2
      CALL NORM(arg1(:), dist)
! Store nodal arc-length
      arclength(elemid+1) = arclength(elemid) + dist
! Store coordinates of the next node
      nodecoor(:, elemid+1) = node2
    END DO
! SAMPLING POSITION FOR NEW NODES
    zero = 0.
    one = 1.
    pi = 3.1415926535897932384626
! Now that we know the initial and final arcLength, we can redistribute the
! parametric coordinates based on the used defined spacing criteria.
! These statements should initially create parametric coordinates in the interval
! [0.0, 1.0]. We will rescale it after the if statements.
    IF (spacing .EQ. 'linear') THEN
      CALL PUSHCONTROL2B(0)
      CALL LINSPACE(zero, one, nnewnodes, newarclength)
    ELSE IF (spacing .EQ. 'cosine') THEN
      CALL PUSHCONTROL2B(1)
      CALL LINSPACE(zero, pi, nnewnodes, newarclength)
      newarclength = 0.5*(1.0-COS(newarclength))
    ELSE IF (spacing .EQ. 'hyptan') THEN
      arg10 = sp1/arclength(nelem+1)
      arg2 = sp2/arclength(nelem+1)
      CALL GETHYPTANDIST(arg10, arg2, nnewnodes, newarclength)
      CALL PUSHCONTROL2B(2)
    ELSE
      CALL PUSHCONTROL2B(3)
    END IF
! Rescale newArcLength based on the final distance
    CALL PUSHREAL4ARRAY(newarclength, realtype*nnewnodes/4)
    newarclength = arclength(nnodes)*newarclength
! INTERPOLATE NEW NODES
! Now we sample the new coordinates based on the interpolation method given by the user
! Create interpolants for x, y, and z
    CALL INTERP1D(1, nnodes, arclength, nodecoor(1, :), nnewnodes, &
&           newarclength, newcoor(1, :))
    CALL INTERP1D(1, nnodes, arclength, nodecoor(2, :), nnewnodes, &
&           newarclength, newcoor(2, :))
    CALL INTERP1D(1, nnodes, arclength, nodecoor(3, :), nnewnodes, &
&           newarclength, newcoor(3, :))
! NODE MERGING
! If a new node is very close to an old node, the derivative of this new node
! is undefined, because it is at the discontinuity between two elements.
! In this case, we will add some extra code to assign the position of these
! old nodes directly to the new nodes.
! We do this so that the AD codes uses the derivative seeds coming from
! the old nodes as well, avoiding the undefined derivative issue.
! Loop over the new nodes to see if they are too close to an old node
    DO ii=1,nnewnodes
! Get coordinates of the new node
      newnode = newcoor(:, ii)
! Loop over the old nodes
      DO jj=1,nnodes
! Get coordinates of the old node
        oldnode = nodecoor(:, jj)
! Compute the distance between nodes
        distvec = newnode - oldnode
        dist = SQRT(distvec(1)**2 + distvec(2)**2 + distvec(3)**2)
! Check if distance is below a threshold
        IF (dist .LT. disttol) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
    END DO
    nodecoorb = 0.0
    DO ii=nnewnodes,1,-1
      DO jj=nnodes,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          nodecoorb(:, jj) = nodecoorb(:, jj) + newcoorb(:, ii)
          newcoorb(:, ii) = 0.0
        END IF
      END DO
    END DO
    arclengthb = 0.0
    newarclengthb = 0.0
    CALL INTERP1D_B(1, nnodes, arclength, arclengthb, nodecoor(3, :), &
&             nodecoorb(3, :), nnewnodes, newarclength, newarclengthb, &
&             newcoor(3, :), newcoorb(3, :))
    CALL INTERP1D_B(1, nnodes, arclength, arclengthb, nodecoor(2, :), &
&             nodecoorb(2, :), nnewnodes, newarclength, newarclengthb, &
&             newcoor(2, :), newcoorb(2, :))
    CALL INTERP1D_B(1, nnodes, arclength, arclengthb, nodecoor(1, :), &
&             nodecoorb(1, :), nnewnodes, newarclength, newarclengthb, &
&             newcoor(1, :), newcoorb(1, :))
    CALL POPREAL4ARRAY(newarclength, realtype*nnewnodes/4)
    arclengthb(nnodes) = arclengthb(nnodes) + SUM(newarclength*&
&     newarclengthb)
    newarclengthb = arclength(nnodes)*newarclengthb
    CALL POPCONTROL2B(branch)
    IF (branch .GE. 2) THEN
      IF (branch .EQ. 2) THEN
        CALL GETHYPTANDIST_B(arg10, arg10b, arg2, arg2b, nnewnodes, &
&                      newarclength, newarclengthb)
        arclengthb(nelem+1) = arclengthb(nelem+1) - sp1*arg10b/arclength&
&         (nelem+1)**2 - sp2*arg2b/arclength(nelem+1)**2
      END IF
    END IF
    coorb = 0.0
    DO elemid=nelem,1,-1
      node2 = coor(:, barsconn(2, elemid))
      node2b = 0.0
      arclengthb(elemid) = arclengthb(elemid) + arclengthb(elemid+1)
      distb = arclengthb(elemid+1)
      arclengthb(elemid+1) = 0.0
      node1 = coor(:, barsconn(1, elemid))
      arg1(:) = node1 - node2
      CALL NORM_B0(arg1(:), arg1b(:), dist, distb)
      node2b = nodecoorb(:, elemid+1) - arg1b(:)
      nodecoorb(:, elemid+1) = 0.0
      node1b = 0.0
      node1b = arg1b(:)
      coorb(:, barsconn(2, elemid)) = coorb(:, barsconn(2, elemid)) + &
&       node2b
      coorb(:, barsconn(1, elemid)) = coorb(:, barsconn(1, elemid)) + &
&       node1b
    END DO
    coorb(:, barsconn(1, 1)) = coorb(:, barsconn(1, 1)) + nodecoorb(:, 1&
&     )
    GOTO 110
 100 CALL PUSHCONTROL1B(1)
    CALL PUSHINTEGER4(ad_count)
 110 CALL POPINTEGER4(ad_count)
    DO i=1,ad_count
      IF (i .EQ. 1) THEN
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) coorb = 0.0
      END IF
    END DO
    newcoorb = 0.0
  END SUBROUTINE REMESH_MAIN_B
  SUBROUTINE REMESH_MAIN(nnodes, nelem, nnewnodes, coor, barsconn, &
&   method, spacing, sp1, sp2, newcoor, newbarsconn)
    IMPLICIT NONE
! Input variables
    INTEGER(kind=inttype), INTENT(IN) :: nnewnodes, nelem
    INTEGER(kind=inttype), INTENT(INOUT) :: nnodes
    CHARACTER(len=32), INTENT(IN) :: method, spacing
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(IN) :: coor
    INTEGER(kind=inttype), DIMENSION(2, nelem), INTENT(IN) :: barsconn
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
! Output variables
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoor
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes - 1) :: newbarsconn
! Working variables
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoor
    REAL(kind=realtype), DIMENSION(nnodes) :: arclength
    INTEGER(kind=inttype) :: elemid, prevnodeid, currnodeid, ii, jj
    REAL(kind=realtype) :: dist, zero, one, pi, disttol
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclength
    REAL(kind=realtype), DIMENSION(3) :: node1, node2, newnode, oldnode
    REAL(kind=realtype), DIMENSION(3) :: distvec
    INTRINSIC COS
    INTRINSIC SQRT
    REAL(kind=realtype), DIMENSION(3) :: arg1
    REAL(kind=realtype) :: arg10
    REAL(kind=realtype) :: arg2
! Tolerance to avoid any interpolation if new node is too close to an original node
    disttol = 1e-7
    nnodes = nelem + 1
! Initialize outputs
    newcoor = 0.0
    newbarsconn = 0
! First we check if the FE data is ordered
    DO elemid=2,nelem
! Get node indices
      prevnodeid = barsconn(2, elemid-1)
      currnodeid = barsconn(1, elemid)
! Check if the FE data is ordered
      IF (prevnodeid .NE. currnodeid) THEN
! Print warning
        PRINT*, &
&    'WARNING: Could not remesh curve because it has unordered FE data.'
        PRINT*, '         Call FEsort first.'
        RETURN
      END IF
    END DO
! COMPUTE ARC-LENGTH
! We can proceed if FE data is ordered
! Store position of the first node (the other nodes will be covered in the loop)
! (the -1 is due Fortran indexing)
    nodecoor(:, 1) = coor(:, barsconn(1, 1))
    arclength(1) = 0.0
! Loop over each element to increment arcLength
    DO elemid=1,nelem
! Get node positions (the -1 is due Fortran indexing)
      node1 = coor(:, barsconn(1, elemid))
      node2 = coor(:, barsconn(2, elemid))
! Compute distance between nodes
      arg1(:) = node1 - node2
      CALL NORM(arg1(:), dist)
! Store nodal arc-length
      arclength(elemid+1) = arclength(elemid) + dist
! Store coordinates of the next node
      nodecoor(:, elemid+1) = node2
    END DO
! SAMPLING POSITION FOR NEW NODES
    zero = 0.
    one = 1.
    pi = 3.1415926535897932384626
! Now that we know the initial and final arcLength, we can redistribute the
! parametric coordinates based on the used defined spacing criteria.
! These statements should initially create parametric coordinates in the interval
! [0.0, 1.0]. We will rescale it after the if statements.
    IF (spacing .EQ. 'linear') THEN
      CALL LINSPACE(zero, one, nnewnodes, newarclength)
    ELSE IF (spacing .EQ. 'cosine') THEN
      CALL LINSPACE(zero, pi, nnewnodes, newarclength)
      newarclength = 0.5*(1.0-COS(newarclength))
    ELSE IF (spacing .EQ. 'hyptan') THEN
      arg10 = sp1/arclength(nelem+1)
      arg2 = sp2/arclength(nelem+1)
      CALL GETHYPTANDIST(arg10, arg2, nnewnodes, newarclength)
    END IF
! Rescale newArcLength based on the final distance
    newarclength = arclength(nnodes)*newarclength
! INTERPOLATE NEW NODES
! Now we sample the new coordinates based on the interpolation method given by the user
! Create interpolants for x, y, and z
    CALL INTERP1D(1, nnodes, arclength, nodecoor(1, :), nnewnodes, &
&           newarclength, newcoor(1, :))
    CALL INTERP1D(1, nnodes, arclength, nodecoor(2, :), nnewnodes, &
&           newarclength, newcoor(2, :))
    CALL INTERP1D(1, nnodes, arclength, nodecoor(3, :), nnewnodes, &
&           newarclength, newcoor(3, :))
! ASSIGN NEW COORDINATES AND CONNECTIVITIES
! Generate new connectivity (the nodes are already in order so we just
! need to assign an ordered set to barsConn).
    DO elemid=1,nnewnodes-1
      newbarsconn(1, elemid) = elemid
      newbarsconn(2, elemid) = elemid + 1
    END DO
! NODE MERGING
! If a new node is very close to an old node, the derivative of this new node
! is undefined, because it is at the discontinuity between two elements.
! In this case, we will add some extra code to assign the position of these
! old nodes directly to the new nodes.
! We do this so that the AD codes uses the derivative seeds coming from
! the old nodes as well, avoiding the undefined derivative issue.
! Loop over the new nodes to see if they are too close to an old node
    DO ii=1,nnewnodes
! Get coordinates of the new node
      newnode = newcoor(:, ii)
! Loop over the old nodes
      DO jj=1,nnodes
! Get coordinates of the old node
        oldnode = nodecoor(:, jj)
! Compute the distance between nodes
        distvec = newnode - oldnode
        dist = SQRT(distvec(1)**2 + distvec(2)**2 + distvec(3)**2)
! Check if distance is below a threshold
        IF (dist .LT. disttol) newcoor(:, ii) = nodecoor(:, jj)
! Repeat the old node to avoid indetermination in derivatives
      END DO
    END DO
  END SUBROUTINE REMESH_MAIN
!============================================================
  SUBROUTINE LINSPACE(l, k, n, z)
    IMPLICIT NONE
!// Argument declarations
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), DIMENSION(n), INTENT(OUT) :: z
    REAL(kind=realtype), INTENT(IN) :: l
    REAL(kind=realtype), INTENT(IN) :: k
!// local variables
    INTEGER(kind=inttype) :: i
    REAL(kind=realtype) :: d
    d = (k-l)/(n-1)
    z(1) = l
    DO i=2,n-1
      z(i) = z(i-1) + d
    END DO
    z(1) = l
    z(n) = k
    RETURN
  END SUBROUTINE LINSPACE
!  Differentiation of interp1d in reverse (adjoint) mode:
!   gradient     of useful results: p_interp p_data t_data t_interp
!   with respect to varying inputs: p_interp p_data t_data t_interp
  SUBROUTINE INTERP1D_B(m, data_num, t_data, t_datab, p_data, p_datab, &
&   interp_num, t_interp, t_interpb, p_interp, p_interpb)
    IMPLICIT NONE
    INTEGER(kind=inttype) :: data_num
    INTEGER(kind=inttype) :: m
    INTEGER(kind=inttype) :: interp_num
    INTEGER(kind=inttype) :: interp
    INTEGER(kind=inttype) :: left
    REAL(kind=realtype) :: p_data(data_num)
    REAL(kind=realtype) :: p_datab(data_num)
    REAL(kind=realtype) :: p_interp(interp_num)
    REAL(kind=realtype) :: p_interpb(interp_num)
    INTEGER(kind=inttype) :: right
    REAL(kind=realtype) :: t
    REAL(kind=realtype) :: tb
    REAL(kind=realtype) :: t_data(data_num)
    REAL(kind=realtype) :: t_datab(data_num)
    REAL(kind=realtype) :: t_interp(interp_num)
    REAL(kind=realtype) :: t_interpb(interp_num)
    REAL(kind=realtype) :: tempb0
    REAL(kind=realtype) :: tempb
    DO interp=1,interp_num
      t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
      CALL PUSHINTEGER4ARRAY(right, inttype/4)
      CALL PUSHINTEGER4ARRAY(left, inttype/4)
      CALL R8VEC_BRACKET(data_num, t_data, t, left, right)
    END DO
    DO interp=interp_num,1,-1
      t = t_interp(interp)
      tempb = p_interpb(interp)/(t_data(right)-t_data(left))
      tempb0 = -(((t_data(right)-t)*p_data(left)+(t-t_data(left))*p_data&
&       (right))*tempb/(t_data(right)-t_data(left)))
      t_datab(right) = t_datab(right) + tempb0 + p_data(left)*tempb
      tb = (p_data(right)-p_data(left))*tempb
      p_datab(left) = p_datab(left) + (t_data(right)-t)*tempb
      t_datab(left) = t_datab(left) - tempb0 - p_data(right)*tempb
      p_datab(right) = p_datab(right) + (t-t_data(left))*tempb
      p_interpb(interp) = 0.0
      CALL POPINTEGER4ARRAY(left, inttype/4)
      CALL POPINTEGER4ARRAY(right, inttype/4)
      t_interpb(interp) = t_interpb(interp) + tb
    END DO
  END SUBROUTINE INTERP1D_B
  SUBROUTINE INTERP1D(m, data_num, t_data, p_data, interp_num, t_interp&
&   , p_interp)
    IMPLICIT NONE
    INTEGER(kind=inttype) :: data_num
    INTEGER(kind=inttype) :: m
    INTEGER(kind=inttype) :: interp_num
    INTEGER(kind=inttype) :: interp
    INTEGER(kind=inttype) :: left
    REAL(kind=realtype) :: p_data(data_num)
    REAL(kind=realtype) :: p_interp(interp_num)
    INTEGER(kind=inttype) :: right
    REAL(kind=realtype) :: t
    REAL(kind=realtype) :: t_data(data_num)
    REAL(kind=realtype) :: t_interp(interp_num)
    DO interp=1,interp_num
      t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
      CALL R8VEC_BRACKET(data_num, t_data, t, left, right)
      p_interp(interp) = ((t_data(right)-t)*p_data(left)+(t-t_data(left)&
&       )*p_data(right))/(t_data(right)-t_data(left))
    END DO
    RETURN
  END SUBROUTINE INTERP1D
  SUBROUTINE R8VEC_BRACKET(n, x, xval, left, right)
    IMPLICIT NONE
    INTEGER(kind=inttype) :: n
    INTEGER(kind=inttype) :: i
    INTEGER(kind=inttype) :: left
    INTEGER(kind=inttype) :: right
    REAL(kind=realtype) :: x(n)
    REAL(kind=realtype) :: xval
    DO i=2,n-1
      IF (xval .LT. x(i)) THEN
        left = i - 1
        right = i
        RETURN
      END IF
    END DO
    left = n - 1
    right = n
    RETURN
  END SUBROUTINE R8VEC_BRACKET
!  Differentiation of dot in reverse (adjoint) mode:
!   gradient     of useful results: dot_ a b
!   with respect to varying inputs: a b
!============================================================
  SUBROUTINE DOT_B0(a, ab, b, bb, dot_, dot_b)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3), b(3)
    REAL(kind=realtype) :: ab(3), bb(3)
    REAL(kind=realtype) :: dot_
    REAL(kind=realtype) :: dot_b
    INTRINSIC SUM
    ab = ab + b*dot_b
    bb = bb + a*dot_b
  END SUBROUTINE DOT_B0
!============================================================
  SUBROUTINE DOT(a, b, dot_)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3), b(3)
    REAL(kind=realtype), INTENT(OUT) :: dot_
    INTRINSIC SUM
    dot_ = SUM(a*b)
  END SUBROUTINE DOT
!  Differentiation of norm in reverse (adjoint) mode:
!   gradient     of useful results: norm_
!   with respect to varying inputs: a
!============================================================
  SUBROUTINE NORM_B0(a, ab, norm_, norm_b)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3)
    REAL(kind=realtype) :: ab(3)
    REAL(kind=realtype) :: norm_
    REAL(kind=realtype) :: norm_b
    INTRINSIC SQRT
    REAL(kind=realtype) :: tempb
    ab = 0.0
    IF (a(1)**2 + a(2)**2 + a(3)**2 .EQ. 0.0) THEN
      tempb = 0.0
    ELSE
      tempb = norm_b/(2.0*SQRT(a(1)**2+a(2)**2+a(3)**2))
    END IF
    ab(1) = ab(1) + 2*a(1)*tempb
    ab(2) = ab(2) + 2*a(2)*tempb
    ab(3) = ab(3) + 2*a(3)*tempb
  END SUBROUTINE NORM_B0
!============================================================
  SUBROUTINE NORM(a, norm_)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3)
    REAL(kind=realtype), INTENT(OUT) :: norm_
    INTRINSIC SQRT
    norm_ = SQRT(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))
  END SUBROUTINE NORM
!  Differentiation of gethyptandist in reverse (adjoint) mode:
!   gradient     of useful results: spacings
!   with respect to varying inputs: sp1 sp2
  SUBROUTINE GETHYPTANDIST_B(sp1, sp1b, sp2, sp2b, n, spacings, &
&   spacingsb)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
    REAL(kind=realtype) :: sp1b, sp2b
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), DIMENSION(n) :: spacings
    REAL(kind=realtype), DIMENSION(n) :: spacingsb
    INTEGER(kind=inttype) :: i
    REAL(kind=realtype) :: b, step, b_out, b_out_step, f_prime, b_old
    REAL(kind=realtype) :: bb, b_outb, b_out_stepb, f_primeb
    REAL(kind=realtype) :: a, u, r
    REAL(kind=realtype) :: ab, ub
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC FLOAT
    INTRINSIC TANH
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: arg1b
    INTEGER :: ad_count
    INTEGER :: i0
    INTEGER :: branch
    REAL(kind=realtype) :: temp0
    REAL(kind=realtype) :: tempb0
    REAL(kind=realtype) :: tempb
    REAL(kind=realtype) :: abs0
    REAL(kind=realtype) :: temp
! Manually created secant method for the above solve
    b = 4.
    step = 1.e-6
    ad_count = 1
    DO i=1,1000
      CALL PUSHREAL4ARRAY(b_out, realtype/4)
      CALL FINDROOTB(b, sp1, sp2, n, b_out)
      arg1 = b - step
      CALL FINDROOTB(arg1, sp1, sp2, n, b_out_step)
      b_old = b
      CALL PUSHREAL4ARRAY(f_prime, realtype/4)
      f_prime = (b_out-b_out_step)/step
      CALL PUSHREAL4ARRAY(b, realtype/4)
      b = b - b_out/f_prime
      IF (b_old - b .GE. 0.) THEN
        abs0 = b_old - b
      ELSE
        abs0 = -(b_old-b)
      END IF
      IF (abs0 .LT. 1.e-10) THEN
        GOTO 100
      ELSE
        ad_count = ad_count + 1
      END IF
    END DO
    CALL PUSHCONTROL1B(0)
    CALL PUSHINTEGER4(ad_count)
    GOTO 110
 100 CALL PUSHCONTROL1B(1)
    CALL PUSHINTEGER4(ad_count)
! Compute parameter A
 110 a = SQRT(sp1/sp2)
    DO i=1,n
      CALL PUSHREAL4ARRAY(r, realtype/4)
      r = FLOAT(i-1)/FLOAT(n-1) - .5
      CALL PUSHREAL4ARRAY(u, realtype/4)
      u = 1 + TANH(b*r)/TANH(b/2)
    END DO
    ab = 0.0
    bb = 0.0
    DO i=n,1,-1
      temp0 = 2*a + (-a+1)*u
      tempb0 = -(u*spacingsb(i)/temp0**2)
      ub = (1-a)*tempb0 + spacingsb(i)/temp0
      ab = ab + (2-u)*tempb0
      spacingsb(i) = 0.0
      CALL POPREAL4ARRAY(u, realtype/4)
      temp = TANH(b/2)
      bb = bb + ((1.0-TANH(r*b)**2)*r/temp-(1.0-TANH(b/2)**2)*TANH(r*b)/&
&       (temp**2*2))*ub
      CALL POPREAL4ARRAY(r, realtype/4)
    END DO
    IF (sp1/sp2 .EQ. 0.0) THEN
      tempb = 0.0
    ELSE
      tempb = ab/(2.0*SQRT(sp1/sp2)*sp2)
    END IF
    sp1b = tempb
    sp2b = -(sp1*tempb/sp2)
    CALL POPINTEGER4(ad_count)
    DO 120 i0=1,ad_count
      IF (i0 .EQ. 1) THEN
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) GOTO 120
      END IF
      f_primeb = b_out*bb/f_prime**2
      CALL POPREAL4ARRAY(b, realtype/4)
      b_outb = f_primeb/step - bb/f_prime
      CALL POPREAL4ARRAY(f_prime, realtype/4)
      b_out_stepb = -(f_primeb/step)
      arg1 = b - step
      arg1b = 0.0
      CALL FINDROOTB_B(arg1, arg1b, sp1, sp1b, sp2, sp2b, n, b_out_step&
&                , b_out_stepb)
      bb = bb + arg1b
      CALL POPREAL4ARRAY(b_out, realtype/4)
      CALL FINDROOTB_B(b, bb, sp1, sp1b, sp2, sp2b, n, b_out, b_outb)
 120 CONTINUE
  END SUBROUTINE GETHYPTANDIST_B
  SUBROUTINE GETHYPTANDIST(sp1, sp2, n, spacings)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), DIMENSION(n), INTENT(OUT) :: spacings
    INTEGER(kind=inttype) :: i
    REAL(kind=realtype) :: b, step, b_out, b_out_step, f_prime, b_old
    REAL(kind=realtype) :: a, u, r
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC FLOAT
    INTRINSIC TANH
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: abs0
! Manually created secant method for the above solve
    b = 4.
    step = 1.e-6
    DO i=1,1000
      CALL FINDROOTB(b, sp1, sp2, n, b_out)
      arg1 = b - step
      CALL FINDROOTB(arg1, sp1, sp2, n, b_out_step)
      b_old = b
      f_prime = (b_out-b_out_step)/step
      b = b - b_out/f_prime
      IF (b_old - b .GE. 0.) THEN
        abs0 = b_old - b
      ELSE
        abs0 = -(b_old-b)
      END IF
      IF (abs0 .LT. 1.e-10) GOTO 100
    END DO
! Compute parameter A
 100 a = SQRT(sp1/sp2)
    DO i=1,n
      r = FLOAT(i-1)/FLOAT(n-1) - .5
      u = 1 + TANH(b*r)/TANH(b/2)
      spacings(i) = u/(2*a+(1-a)*u)
    END DO
  END SUBROUTINE GETHYPTANDIST
!  Differentiation of findrootb in reverse (adjoint) mode:
!   gradient     of useful results: b_out sp1 sp2 b
!   with respect to varying inputs: sp1 sp2 b
  SUBROUTINE FINDROOTB_B(b, bb, sp1, sp1b, sp2, sp2b, n, b_out, b_outb)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: b, sp1, sp2
    REAL(kind=realtype) :: bb, sp1b, sp2b
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype) :: b_out
    REAL(kind=realtype) :: b_outb
    INTRINSIC SINH
    INTRINSIC SQRT
    REAL(kind=realtype) :: tempb0
    REAL(kind=realtype) :: tempb
    REAL(kind=realtype) :: temp
    temp = SQRT(sp1*sp2)
    tempb = -(b_outb/((n-1)*temp))
    IF (sp1*sp2 .EQ. 0.0) THEN
      tempb0 = 0.0
    ELSE
      tempb0 = -(b*tempb/(temp**2*2.0))
    END IF
    bb = bb + tempb + COSH(b)*b_outb
    sp1b = sp1b + sp2*tempb0
    sp2b = sp2b + sp1*tempb0
  END SUBROUTINE FINDROOTB_B
  SUBROUTINE FINDROOTB(b, sp1, sp2, n, b_out)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: b, sp1, sp2
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), INTENT(OUT) :: b_out
    INTRINSIC SINH
    INTRINSIC SQRT
    b_out = SINH(b) - b/(n-1)/SQRT(sp1*sp2)
  END SUBROUTINE FINDROOTB
!============================================================
! BOUNDING BOX ROUTINES
!============================================================
  SUBROUTINE COMPUTEBBOX(coor, bbox)
    IMPLICIT NONE
! INPUTS
    REAL(kind=realtype), DIMENSION(:, :), INTENT(IN) :: coor
! OUTPUTS
    REAL(kind=realtype), DIMENSION(3, 2), INTENT(OUT) :: bbox
    INTRINSIC MINVAL
    INTRINSIC MAXVAL
! EXECUTION
! Get bounding values
    bbox(:, 1) = MINVAL(coor, 2)
    bbox(:, 2) = MAXVAL(coor, 2)
  END SUBROUTINE COMPUTEBBOX
!============================================================
  SUBROUTINE COMPUTEBBOXPERELEMENTS(nnodes, ntria, nquads, coor, &
&   triaconn, quadsconn, triabbox, quadsbbox)
    IMPLICIT NONE
! INPUTS
    INTEGER(kind=inttype), INTENT(IN) :: nnodes, ntria, nquads
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(IN) :: coor
    INTEGER(kind=inttype), DIMENSION(3, ntria), INTENT(IN) :: triaconn
    INTEGER(kind=inttype), DIMENSION(4, nquads), INTENT(IN) :: quadsconn
! OUTPUTS
    REAL(kind=realtype), DIMENSION(6, ntria), INTENT(OUT) :: triabbox
    REAL(kind=realtype), DIMENSION(6, nquads), INTENT(OUT) :: quadsbbox
! WORKING
    INTEGER(kind=inttype) :: elemid, nodeid
    REAL(kind=realtype) :: triacoor(3, 3), quadscoor(3, 4), bbox(6)
    INTRINSIC REAL
    INTRINSIC MINVAL
    INTRINSIC MAXVAL
    REAL, DIMENSION(3, 3) :: arg1
    REAL, DIMENSION(3, 4) :: arg10
! EXECUTION
! Loop over all triangles
    DO elemid=1,ntria
! Get coordinates of each node and store them in a single matrix
      triacoor(:, 1) = coor(:, triaconn(1, elemid))
      triacoor(:, 2) = coor(:, triaconn(2, elemid))
      triacoor(:, 3) = coor(:, triaconn(3, elemid))
! Assign min values (BBox(1:3)) and max values (BBox(4:6))
! based on the nodal coordinates.
      arg1(:, :) = REAL(triacoor)
      triabbox(1:3, elemid) = MINVAL(arg1(:, :), 2)
      arg1(:, :) = REAL(triacoor)
      triabbox(4:6, elemid) = MAXVAL(arg1(:, :), 2)
    END DO
! Loop over all quads
    DO elemid=1,nquads
! Get coordinates of each node and store them in a single matrix
      quadscoor(:, 1) = coor(:, quadsconn(1, elemid))
      quadscoor(:, 2) = coor(:, quadsconn(2, elemid))
      quadscoor(:, 3) = coor(:, quadsconn(3, elemid))
      quadscoor(:, 4) = coor(:, quadsconn(4, elemid))
! Assign min values (BBox(1:3)) and max values (BBox(4:6))
! based on the nodal coordinates.
      arg10(:, :) = REAL(quadscoor)
      quadsbbox(1:3, elemid) = MINVAL(arg10(:, :), 2)
      arg10(:, :) = REAL(quadscoor)
      quadsbbox(4:6, elemid) = MAXVAL(arg10(:, :), 2)
    END DO
  END SUBROUTINE COMPUTEBBOXPERELEMENTS
!============================================================
  SUBROUTINE COMPUTEBBOXINTERSECTION(bboxa, bboxb, bboxab, overlap)
    IMPLICIT NONE
! Remember that overlap may be set to false in the Z overlap test
! INPUTS
    REAL(kind=realtype), DIMENSION(6), INTENT(IN) :: bboxa, bboxb
! OUTPUTS
    REAL(kind=realtype), DIMENSION(6), INTENT(OUT) :: bboxab
    LOGICAL, INTENT(OUT) :: overlap
! WORKING
    REAL(kind=realtype), DIMENSION(6) :: bounds
! EXECUTION
! Check overlaps along each dimension
! We have an actual BBox intersection if there are overlaps in all dimensions
! X overlap
    CALL LINEINTERSECTIONINTERVAL(bboxa(1), bboxa(4), bboxb(1), bboxb(4)&
&                           , bboxab(1), bboxab(4), overlap)
    IF (overlap) THEN
! Y overlap
      CALL LINEINTERSECTIONINTERVAL(bboxa(2), bboxa(5), bboxb(2), bboxb(&
&                             5), bboxab(2), bboxab(5), overlap)
      IF (overlap) CALL LINEINTERSECTIONINTERVAL(bboxa(3), bboxa(6), &
&                                          bboxb(3), bboxb(6), bboxab(3)&
&                                          , bboxab(6), overlap)
! Z overlap
    END IF
! Determine the size of the edges of the bounding box
    bounds(1) = bboxab(4) - bboxab(1)
    bounds(2) = bboxab(5) - bboxab(2)
    bounds(3) = bboxab(6) - bboxab(3)
    bounds(4:6) = -bounds(1:3)
! Add buffer space to the bounding box.
! Although this slightly slows down the intersection algorithm, it is on the
! order of hundredths of seconds for the CRM case and solves a small bug.
! Will need to fix the filterElements code later so we can use a tighter
! bounding box.
    bboxab = bboxab - 0.01*bounds
  END SUBROUTINE COMPUTEBBOXINTERSECTION
!============================================================
  SUBROUTINE LINEINTERSECTIONINTERVAL(xmina, xmaxa, xminb, xmaxb, xminab&
&   , xmaxab, overlap)
    IMPLICIT NONE
! INPUTS
    REAL(kind=realtype), INTENT(IN) :: xmina, xmaxa, xminb, xmaxb
! OUTPUTS
    REAL(kind=realtype), INTENT(OUT) :: xminab, xmaxab
    LOGICAL, INTENT(OUT) :: overlap
    INTRINSIC MAX
    INTRINSIC MIN
    IF (xmina .LT. xminb) THEN
      xminab = xminb
    ELSE
      xminab = xmina
    END IF
    IF (xmaxa .GT. xmaxb) THEN
      xmaxab = xmaxb
    ELSE
      xmaxab = xmaxa
    END IF
! Check if we actually have an overlap
    IF (xminab .GT. xmaxab) THEN
      overlap = .false.
    ELSE
      overlap = .true.
    END IF
  END SUBROUTINE LINEINTERSECTIONINTERVAL
END MODULE UTILITIES_B
