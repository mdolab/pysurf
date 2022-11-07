!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
MODULE UTILITIES_D
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
        CALL NORM(currcoor - prevcoor, dist)
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

!  Differentiation of remesh_main in forward (tangent) mode:
!   variations   of useful results: newcoor
!   with respect to varying inputs: coor
!   RW status of diff variables: coor:in newcoor:out
  SUBROUTINE REMESH_MAIN_D(nnodes, nelem, nnewnodes, coor, coord, &
&   barsconn, method, spacing, sp1, sp2, newcoor, newcoord, newbarsconn)
    IMPLICIT NONE
! Input variables
    INTEGER(kind=inttype), INTENT(IN) :: nnewnodes, nelem
    INTEGER(kind=inttype), INTENT(INOUT) :: nnodes
    CHARACTER(len=32), INTENT(IN) :: method, spacing
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(IN) :: coor
    REAL(kind=realtype), DIMENSION(3, nnodes), INTENT(IN) :: coord
    INTEGER(kind=inttype), DIMENSION(2, nelem), INTENT(IN) :: barsconn
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
! Output variables
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoor
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoord
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes-1) :: newbarsconn
! Working variables
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoor
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoord
    REAL(kind=realtype), DIMENSION(nnodes) :: arclength
    REAL(kind=realtype), DIMENSION(nnodes) :: arclengthd
    INTEGER(kind=inttype) :: elemid, prevnodeid, currnodeid, ii, jj
    REAL(kind=realtype) :: dist, zero, one, pi, disttol
    REAL(kind=realtype) :: distd
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclength
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclengthd
    REAL(kind=realtype), DIMENSION(3) :: node1, node2, newnode, oldnode
    REAL(kind=realtype), DIMENSION(3) :: node1d, node2d
    REAL(kind=realtype), DIMENSION(3) :: distvec
    INTRINSIC COS
    INTRINSIC SQRT
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: temp
    REAL(kind=realtype) :: temp0
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
      IF (prevnodeid .NE. currnodeid) GOTO 100
    END DO
! COMPUTE ARC-LENGTH
! We can proceed if FE data is ordered
! Store position of the first node (the other nodes will be covered in the loop)
! (the -1 is due Fortran indexing)
    nodecoord = 0.0_8
    nodecoord(:, 1) = coord(:, barsconn(1, 1))
    nodecoor(:, 1) = coor(:, barsconn(1, 1))
    arclength(1) = 0.0
    arclengthd = 0.0_8
! Loop over each element to increment arcLength
    DO elemid=1,nelem
! Get node positions (the -1 is due Fortran indexing)
      node1d = coord(:, barsconn(1, elemid))
      node1 = coor(:, barsconn(1, elemid))
      node2d = coord(:, barsconn(2, elemid))
      node2 = coor(:, barsconn(2, elemid))
! Compute distance between nodes
      CALL NORM_D0(node1 - node2, node1d - node2d, dist, distd)
! Store nodal arc-length
      arclengthd(elemid+1) = arclengthd(elemid) + distd
      arclength(elemid+1) = arclength(elemid) + dist
! Store coordinates of the next node
      nodecoord(:, elemid+1) = node2d
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
      newarclengthd = 0.0_8
    ELSE IF (spacing .EQ. 'cosine') THEN
      CALL LINSPACE(zero, pi, nnewnodes, newarclength)
      newarclength = 0.5*(1.0-COS(newarclength))
      newarclengthd = 0.0_8
    ELSE IF (spacing .EQ. 'hyptan') THEN
      temp = sp2/arclength(nelem+1)
      temp0 = sp1/arclength(nelem+1)
      CALL GETHYPTANDIST_D(sp1/arclength(nelem+1), -(temp0*arclengthd(&
&                    nelem+1)/arclength(nelem+1)), sp2/arclength(nelem+1&
&                    ), -(temp*arclengthd(nelem+1)/arclength(nelem+1)), &
&                    nnewnodes, newarclength, newarclengthd)
    ELSE
      newarclengthd = 0.0_8
    END IF
! Rescale newArcLength based on the final distance
    newarclengthd = newarclength*arclengthd(nnodes) + arclength(nnodes)*&
&     newarclengthd
    newarclength = arclength(nnodes)*newarclength
! INTERPOLATE NEW NODES
! Now we sample the new coordinates based on the interpolation method given by the user
! Create interpolants for x, y, and z
    newcoord = 0.0_8
    CALL INTERP1D_D(1, nnodes, arclength, arclengthd, nodecoor(1, :), &
&             nodecoord(1, :), nnewnodes, newarclength, newarclengthd, &
&             newcoor(1, :), newcoord(1, :))
    CALL INTERP1D_D(1, nnodes, arclength, arclengthd, nodecoor(2, :), &
&             nodecoord(2, :), nnewnodes, newarclength, newarclengthd, &
&             newcoor(2, :), newcoord(2, :))
    CALL INTERP1D_D(1, nnodes, arclength, arclengthd, nodecoor(3, :), &
&             nodecoord(3, :), nnewnodes, newarclength, newarclengthd, &
&             newcoor(3, :), newcoord(3, :))
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
        arg1 = distvec(1)**2 + distvec(2)**2 + distvec(3)**2
        dist = SQRT(arg1)
! Check if distance is below a threshold
        IF (dist .LT. disttol) THEN
! Repeat the old node to avoid indetermination in derivatives
          newcoord(:, ii) = nodecoord(:, jj)
          newcoor(:, ii) = nodecoor(:, jj)
        END IF
      END DO
    END DO
    GOTO 110
! Print warning
 100 PRINT*, &
&    'WARNING: Could not remesh curve because it has unordered FE data.'
    PRINT*, '         Call FEsort first.'
    newcoord = 0.0_8
    RETURN
 110 CONTINUE
  END SUBROUTINE REMESH_MAIN_D

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
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes-1) :: newbarsconn
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
    REAL(kind=realtype) :: arg1
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
      CALL NORM(node1 - node2, dist)
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
      CALL GETHYPTANDIST(sp1/arclength(nelem+1), sp2/arclength(nelem+1)&
&                  , nnewnodes, newarclength)
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
        arg1 = distvec(1)**2 + distvec(2)**2 + distvec(3)**2
        dist = SQRT(arg1)
! Check if distance is below a threshold
        IF (dist .LT. disttol) THEN
! Repeat the old node to avoid indetermination in derivatives
          newcoor(:, ii) = nodecoor(:, jj)
        END IF
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

!  Differentiation of interp1d in forward (tangent) mode:
!   variations   of useful results: p_interp
!   with respect to varying inputs: p_interp p_data t_data t_interp
  SUBROUTINE INTERP1D_D(m, data_num, t_data, t_datad, p_data, p_datad, &
&   interp_num, t_interp, t_interpd, p_interp, p_interpd)
    IMPLICIT NONE
    INTEGER(kind=inttype) :: data_num
    INTEGER(kind=inttype) :: m
    INTEGER(kind=inttype) :: interp_num
    INTEGER(kind=inttype) :: interp
    INTEGER(kind=inttype) :: left
    REAL(kind=realtype) :: p_data(data_num)
    REAL(kind=realtype) :: p_datad(data_num)
    REAL(kind=realtype) :: p_interp(interp_num)
    REAL(kind=realtype) :: p_interpd(interp_num)
    INTEGER(kind=inttype) :: right
    REAL(kind=realtype) :: t
    REAL(kind=realtype) :: td
    REAL(kind=realtype) :: t_data(data_num)
    REAL(kind=realtype) :: t_datad(data_num)
    REAL(kind=realtype) :: t_interp(interp_num)
    REAL(kind=realtype) :: t_interpd(interp_num)
    REAL(kind=realtype) :: temp
    DO interp=1,interp_num
      td = t_interpd(interp)
      t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
      CALL R8VEC_BRACKET(data_num, t_data, t, left, right)
      temp = ((t_data(right)-t)*p_data(left)+(t-t_data(left))*p_data(&
&       right))/(t_data(right)-t_data(left))
      p_interpd(interp) = (p_data(left)*(t_datad(right)-td)+(t_data(&
&       right)-t)*p_datad(left)+p_data(right)*(td-t_datad(left))+(t-&
&       t_data(left))*p_datad(right)-temp*(t_datad(right)-t_datad(left))&
&       )/(t_data(right)-t_data(left))
      p_interp(interp) = temp
    END DO
    RETURN
  END SUBROUTINE INTERP1D_D

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

!  Differentiation of dot in forward (tangent) mode:
!   variations   of useful results: dot_
!   with respect to varying inputs: a b
!============================================================
  SUBROUTINE DOT_D0(a, ad, b, bd, dot_, dot_d)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3), b(3)
    REAL(kind=realtype), INTENT(IN) :: ad(3), bd(3)
    REAL(kind=realtype), INTENT(OUT) :: dot_
    REAL(kind=realtype), INTENT(OUT) :: dot_d
    INTRINSIC SUM
    dot_d = SUM(b*ad + a*bd)
    dot_ = SUM(a*b)
  END SUBROUTINE DOT_D0

!============================================================
  SUBROUTINE DOT(a, b, dot_)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3), b(3)
    REAL(kind=realtype), INTENT(OUT) :: dot_
    INTRINSIC SUM
    dot_ = SUM(a*b)
  END SUBROUTINE DOT

!  Differentiation of norm in forward (tangent) mode:
!   variations   of useful results: norm_
!   with respect to varying inputs: a
!============================================================
  SUBROUTINE NORM_D0(a, ad, norm_, norm_d)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3)
    REAL(kind=realtype), INTENT(IN) :: ad(3)
    REAL(kind=realtype), INTENT(OUT) :: norm_
    REAL(kind=realtype), INTENT(OUT) :: norm_d
    INTRINSIC SQRT
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: arg1d
    REAL(kind=realtype) :: temp
    arg1d = 2*a(1)*ad(1) + 2*a(2)*ad(2) + 2*a(3)*ad(3)
    arg1 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    temp = SQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      norm_d = 0.0_8
    ELSE
      norm_d = arg1d/(2.0*temp)
    END IF
    norm_ = temp
  END SUBROUTINE NORM_D0

!============================================================
  SUBROUTINE NORM(a, norm_)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3)
    REAL(kind=realtype), INTENT(OUT) :: norm_
    INTRINSIC SQRT
    REAL(kind=realtype) :: arg1
    arg1 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    norm_ = SQRT(arg1)
  END SUBROUTINE NORM

!  Differentiation of gethyptandist in forward (tangent) mode:
!   variations   of useful results: spacings
!   with respect to varying inputs: sp1 sp2
  SUBROUTINE GETHYPTANDIST_D(sp1, sp1d, sp2, sp2d, n, spacings, &
&   spacingsd)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: sp1, sp2
    REAL(kind=realtype), INTENT(IN) :: sp1d, sp2d
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), DIMENSION(n), INTENT(OUT) :: spacings
    REAL(kind=realtype), DIMENSION(n), INTENT(OUT) :: spacingsd
    INTEGER(kind=inttype) :: i
    REAL(kind=realtype) :: b, step, b_out, b_out_step, f_prime, b_old
    REAL(kind=realtype) :: bd, b_outd, b_out_stepd, f_primed
    REAL(kind=realtype) :: a, u, r
    REAL(kind=realtype) :: ad, ud
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC FLOAT
    INTRINSIC TANH
    REAL(kind=realtype) :: abs0
    REAL(kind=realtype) :: temp
    REAL(kind=realtype) :: temp0
! Manually created secant method for the above solve
    b = 4.
    step = 1.e-6
    bd = 0.0_8
    DO i=1,1000
      CALL FINDROOTB_D(b, bd, sp1, sp1d, sp2, sp2d, n, b_out, b_outd)
      CALL FINDROOTB_D(b - step, bd, sp1, sp1d, sp2, sp2d, n, b_out_step&
&                , b_out_stepd)
      b_old = b
      f_primed = (b_outd-b_out_stepd)/step
      f_prime = (b_out-b_out_step)/step
      bd = bd - (b_outd-b_out*f_primed/f_prime)/f_prime
      b = b - b_out/f_prime
      IF (b_old - b .GE. 0.) THEN
        abs0 = b_old - b
      ELSE
        abs0 = -(b_old-b)
      END IF
      IF (abs0 .LT. 1.e-10) EXIT
    END DO
! Compute parameter A
    temp = sp1/sp2
    temp0 = SQRT(temp)
    IF (temp .EQ. 0.0) THEN
      ad = 0.0_8
    ELSE
      ad = (sp1d-temp*sp2d)/(2.0*temp0*sp2)
    END IF
    a = temp0
    spacingsd = 0.0_8
    DO i=1,n
      r = FLOAT(i-1)/FLOAT(n-1) - .5
      temp0 = TANH(b/2)
      temp = TANH(r*b)/temp0
      ud = ((1.0-TANH(r*b)**2)*r-temp*(1.0-TANH(b/2)**2)/2)*bd/temp0
      u = temp + 1
      temp0 = 2*a + (-a+1)*u
      spacingsd(i) = (ud-u*((2-u)*ad+(1-a)*ud)/temp0)/temp0
      spacings(i) = u/temp0
    END DO
  END SUBROUTINE GETHYPTANDIST_D

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
    REAL(kind=realtype) :: abs0
! Manually created secant method for the above solve
    b = 4.
    step = 1.e-6
    DO i=1,1000
      CALL FINDROOTB(b, sp1, sp2, n, b_out)
      CALL FINDROOTB(b - step, sp1, sp2, n, b_out_step)
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

!  Differentiation of findrootb in forward (tangent) mode:
!   variations   of useful results: b_out
!   with respect to varying inputs: sp1 sp2 b
  SUBROUTINE FINDROOTB_D(b, bd, sp1, sp1d, sp2, sp2d, n, b_out, b_outd)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: b, sp1, sp2
    REAL(kind=realtype), INTENT(IN) :: bd, sp1d, sp2d
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), INTENT(OUT) :: b_out
    REAL(kind=realtype), INTENT(OUT) :: b_outd
    INTRINSIC SINH
    INTRINSIC SQRT
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: arg1d
    REAL(kind=realtype) :: result1
    REAL(kind=realtype) :: result1d
    REAL(kind=realtype) :: temp
    arg1d = sp2*sp1d + sp1*sp2d
    arg1 = sp1*sp2
    temp = SQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.0_8
    ELSE
      result1d = arg1d/(2.0*temp)
    END IF
    result1 = temp
    temp = b/((n-1)*result1)
    b_outd = COSH(b)*bd - (bd-temp*(n-1)*result1d)/((n-1)*result1)
    b_out = SINH(b) - temp
  END SUBROUTINE FINDROOTB_D

  SUBROUTINE FINDROOTB(b, sp1, sp2, n, b_out)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: b, sp1, sp2
    INTEGER(kind=inttype), INTENT(IN) :: n
    REAL(kind=realtype), INTENT(OUT) :: b_out
    INTRINSIC SINH
    INTRINSIC SQRT
    REAL(kind=realtype) :: arg1
    REAL(kind=realtype) :: result1
    arg1 = sp1*sp2
    result1 = SQRT(arg1)
    b_out = SINH(b) - b/(n-1)/result1
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
      IF (overlap) THEN
! Z overlap
        CALL LINEINTERSECTIONINTERVAL(bboxa(3), bboxa(6), bboxb(3), &
&                               bboxb(6), bboxab(3), bboxab(6), overlap)
      END IF
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

END MODULE UTILITIES_D

