!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
!
MODULE UTILITIES_D
  USE PRECISION
  IMPLICIT NONE

CONTAINS
!============================================================
  SUBROUTINE CONDENSEBARNODES_MAIN(nnodes, nelem, disttol, coor, &
&   barsconn, nuniquenodes)
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
! WORKING
    INTEGER(kind=inttype) :: ncopies
    INTEGER(kind=inttype) :: currnodeid, prevnodeid, link, elemid
    REAL(kind=realtype), DIMENSION(3) :: currcoor, prevcoor
    REAL(kind=realtype) :: dist
    INTEGER(kind=inttype), DIMENSION(SIZE(coor, 2)) :: linkold2new
    INTRINSIC SIZE
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
! Initialize number of nodes copied so far
    ncopies = 0
! We loop once again over the nodes so we can copy the unique values
    DO currnodeid=1,nnodes
! Get index of the current node in the new coordinate array
      link = linkold2new(currnodeid)
! Check if the new link is already used
      IF (link .GT. ncopies) THEN
! Get coordinates of current node
        currcoor = coor(:, currnodeid)
! Increment number of copies done so far
        ncopies = ncopies + 1
! Copy coordinates
        coor(:, ncopies) = currcoor
      END IF
    END DO
! Now zero out the unused point in coor
    coor(:, ncopies+1:nnodes) = 0.0
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
!   with respect to varying inputs: coor newcoor
!   RW status of diff variables: coor:in newcoor:in-out
  SUBROUTINE REMESH_MAIN_D(nnodes, nnewnodes, coor, coord, barsconn, &
&   method, spacing, newcoor, newcoord, newbarsconn)
    IMPLICIT NONE
! Input variables
    INTEGER(kind=inttype), INTENT(IN) :: nnodes, nnewnodes
    CHARACTER(len=32), INTENT(IN) :: method, spacing
    REAL(kind=realtype), DIMENSION(3, nnodes) :: coor
    REAL(kind=realtype), DIMENSION(3, nnodes) :: coord
    INTEGER(kind=inttype), DIMENSION(2, nnodes - 1) :: barsconn
! Output variables
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoor
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoord
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes - 1) :: newbarsconn
! Working variables
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoor
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoord
    REAL(kind=realtype), DIMENSION(nnodes) :: arclength
    REAL(kind=realtype), DIMENSION(nnodes) :: arclengthd
    INTEGER(kind=inttype) :: elemid, prevnodeid, currnodeid, nelem
    REAL(kind=realtype) :: dist
    REAL(kind=realtype) :: distd
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclength
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclengthd
    REAL(kind=realtype), DIMENSION(3) :: node1, node2
    REAL(kind=realtype), DIMENSION(3) :: node1d, node2d
    LOGICAL :: periodic
    INTRINSIC COS
    nelem = nnodes - 1
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
    nodecoord = 0.0
    nodecoord(:, 1) = coord(:, barsconn(1, 1))
    nodecoor(:, 1) = coor(:, barsconn(1, 1))
    arclengthd = 0.0
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
! Now that we know the initial and final arcLength, we can redistribute the
! parametric coordinates based on the used defined spacing criteria.
! These statements should initially create parametric coordinates in the interval
! [0.0, 1.0]. We will rescale it after the if statements.
    IF (spacing .EQ. 'linear') THEN
      CALL LINSPACE(0.0_8, 1.0_8, nnewnodes, newarclength)
    ELSE IF (spacing .EQ. 'cosine') THEN
      CALL LINSPACE(0.0_8, 3.141592653589793_8, nnewnodes, newarclength)
      newarclength = 0.5*(1.0-COS(newarclength))
    END IF
! Rescale newArcLength based on the final distance
    newarclengthd = newarclength*arclengthd(nnodes)
    newarclength = arclength(nnodes)*newarclength
! INTERPOLATE NEW NODES
! Now we sample the new coordinates based on the interpolation method given by the user
! Create interpolants for x, y, and z
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
! Check if the baseline curve is periodic
    IF (barsconn(1, 1) .EQ. barsconn(2, nnodes-1)) THEN
      periodic = .true.
    ELSE
      periodic = .false.
    END IF
! Generate new connectivity (the nodes are already in order so we just
! need to assign an ordered set to barsConn).
    DO elemid=1,nnewnodes-1
      newbarsconn(1, elemid) = elemid
      newbarsconn(2, elemid) = elemid + 1
    END DO
! We still need to keep periodicity if the original curve is periodic
    IF (periodic) newbarsconn(2, nnewnodes-1) = newbarsconn(1, 1)
    GOTO 110
! Print warning
 100 PRINT*, &
&    'WARNING: Could not remesh curve because it has unordered FE data.'
    PRINT*, '         Call FEsort first.'
    RETURN
 110 CONTINUE
  END SUBROUTINE REMESH_MAIN_D
  SUBROUTINE REMESH_MAIN(nnodes, nnewnodes, coor, barsconn, method, &
&   spacing, newcoor, newbarsconn)
    IMPLICIT NONE
! Input variables
    INTEGER(kind=inttype), INTENT(IN) :: nnodes, nnewnodes
    CHARACTER(len=32), INTENT(IN) :: method, spacing
    REAL(kind=realtype), DIMENSION(3, nnodes) :: coor
    INTEGER(kind=inttype), DIMENSION(2, nnodes - 1) :: barsconn
! Output variables
    REAL(kind=realtype), DIMENSION(3, nnewnodes) :: newcoor
    INTEGER(kind=inttype), DIMENSION(2, nnewnodes - 1) :: newbarsconn
! Working variables
    REAL(kind=realtype), DIMENSION(3, nnodes) :: nodecoor
    REAL(kind=realtype), DIMENSION(nnodes) :: arclength
    INTEGER(kind=inttype) :: elemid, prevnodeid, currnodeid, nelem
    REAL(kind=realtype) :: dist
    REAL(kind=realtype), DIMENSION(nnewnodes) :: newarclength
    REAL(kind=realtype), DIMENSION(3) :: node1, node2
    LOGICAL :: periodic
    INTRINSIC COS
    nelem = nnodes - 1
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
! Now that we know the initial and final arcLength, we can redistribute the
! parametric coordinates based on the used defined spacing criteria.
! These statements should initially create parametric coordinates in the interval
! [0.0, 1.0]. We will rescale it after the if statements.
    IF (spacing .EQ. 'linear') THEN
      CALL LINSPACE(0.0_8, 1.0_8, nnewnodes, newarclength)
    ELSE IF (spacing .EQ. 'cosine') THEN
      CALL LINSPACE(0.0_8, 3.141592653589793_8, nnewnodes, newarclength)
      newarclength = 0.5*(1.0-COS(newarclength))
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
! Check if the baseline curve is periodic
    IF (barsconn(1, 1) .EQ. barsconn(2, nnodes-1)) THEN
      periodic = .true.
    ELSE
      periodic = .false.
    END IF
! Generate new connectivity (the nodes are already in order so we just
! need to assign an ordered set to barsConn).
    DO elemid=1,nnewnodes-1
      newbarsconn(1, elemid) = elemid
      newbarsconn(2, elemid) = elemid + 1
    END DO
! We still need to keep periodicity if the original curve is periodic
    IF (periodic) newbarsconn(2, nnewnodes-1) = newbarsconn(1, 1)
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
    DO interp=1,interp_num
      td = t_interpd(interp)
      t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
      CALL R8VEC_BRACKET(data_num, t_data, t, left, right)
      p_interpd(interp) = (((t_datad(right)-td)*p_data(left)+(t_data(&
&       right)-t)*p_datad(left)+(td-t_datad(left))*p_data(right)+(t-&
&       t_data(left))*p_datad(right))*(t_data(right)-t_data(left))-((&
&       t_data(right)-t)*p_data(left)+(t-t_data(left))*p_data(right))*(&
&       t_datad(right)-t_datad(left)))/(t_data(right)-t_data(left))**2
      p_interp(interp) = ((t_data(right)-t)*p_data(left)+(t-t_data(left)&
&       )*p_data(right))/(t_data(right)-t_data(left))
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
!============================================================
  SUBROUTINE DOT(a, b, dot_)
    IMPLICIT NONE
    REAL(kind=realtype), INTENT(IN) :: a(3), b(3)
    REAL(kind=realtype), INTENT(OUT) :: dot_
    dot_ = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
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
    arg1d = ad(1)*a(1) + a(1)*ad(1) + ad(2)*a(2) + a(2)*ad(2) + ad(3)*a(&
&     3) + a(3)*ad(3)
    arg1 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    IF (arg1 .EQ. 0.0) THEN
      norm_d = 0.0
    ELSE
      norm_d = arg1d/(2.0*SQRT(arg1))
    END IF
    norm_ = SQRT(arg1)
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
END MODULE UTILITIES_D
