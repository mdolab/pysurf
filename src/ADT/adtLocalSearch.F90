!
!     ******************************************************************
!     *                                                                *
!     * File:          adtLocalSearch.f90                              *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-10-2006                                      *
!     * Last modified: 03-17-2006                                      *
!     *                                                                *
!     ******************************************************************
!
module adtLocalSearch
!
!     ******************************************************************
!     *                                                                *
!     * Module which contains the subroutines to perform the local     *
!     * searches, i.e. the tree traversals.                            *
!     *                                                                *
!     ******************************************************************
!
    use adtUtils
    use adtProjections
    implicit none

    !=================================================================

contains

    !===============================================================

    subroutine containmentTreeSearch(jj, coor, &
                                     intInfo, uvw, &
                                     arrDonor, nCoor, &
                                     nInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the actual containment search in the   *
!       * local tree. It is a local routine in the sense that no       *
!       * communication is involved.                                   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * jj:        Entry in the array ADTs, whose ADT must be        *
!       *            searched.                                         *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates of these points.                  *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 2D integer array, in which the following output     *
!       *          be be stored:                                       *
!       *          intInfo(1,:): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2,:): The element type of the element.      *
!       *          intInfo(3,:): The element ID of the element in the  *
!       *                        connectivity.                         *
!       * uvw:     2D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          and the interpolated data:                          *
!       *          uvw(1, :): Parametric u-weight.                     *
!       *          uvw(2, :): Parametric v-weight.                     *
!       *          uvw(3, :): Parametric w-weight.                     *
!       *          uvw(4:,:): Interpolated solution, if desired. It is *
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, jj
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(:, :), intent(in) :: coor
        real(kind=realType), dimension(:, :), intent(in) :: arrDonor

        integer(kind=intType), dimension(:, :), intent(out) :: intInfo
        real(kind=realType), dimension(:, :), intent(out) :: uvw
!
!       Local parameters used in the Newton algorithm.
!
        integer(kind=intType), parameter :: iterMax = 15
        real(kind=realType), parameter :: adtEps = 1.e-25_realType
        real(kind=realType), parameter :: thresConv = 1.e-10_realType
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront
        integer(kind=intType) :: i, nNodeElement

        integer(kind=intType), dimension(8) :: n

        integer(kind=intType), dimension(:), pointer :: BB
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew

        real(kind=realType) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
        real(kind=realType) :: oneMinusU, oneMinusV, oneMinusW
        real(kind=realType) :: oneMinusUMinusV
        real(kind=realType) :: a11, a12, a13, a21, a22, a23
        real(kind=realType) :: a31, a32, a33, val

        real(kind=realType), dimension(3) :: x, f
        real(kind=realType), dimension(8) :: weight
        real(kind=realType), dimension(3, 2:8) :: xn

        real(kind=realType), dimension(:, :), pointer :: xBBox

        logical :: elementFound

        type(adtLeafType), dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox => ADTs(jj)%xBBox
        ADTree => ADTs(jj)%ADTree

        ! Initial allocation of the arrays for the tree traversal.

        nAllocBB = 10
        nAllocFront = 25

        allocate (BB(nAllocBB), frontLeaves(nAllocFront), &
                  frontLeavesNew(nAllocFront), stat=ierr)
        if (ierr /= 0) &
          call adtTerminate(jj, "containmentTreeSearch", &
                            "Memory allocation failure for BB, &
                            &frontLeaves and frontLeavesNew.")

        ! Loop over the number of coordinates to be treated.

        coorLoop: do nn = 1, nCoor

            ! Initialize the processor ID to -1 to indicate that no
            ! corresponding volume element is found.

            intInfo(1, nn) = -1
!
!         **************************************************************
!         *                                                            *
!         * Part 1. Traverse the tree and determine the target         *
!         *         bounding boxes, which may contain the element.     *
!         *                                                            *
!         **************************************************************
!
            ! Start at the root, i.e. set the front leaf to the root leaf.
            ! Also initialize the number of possible bounding boxes to 0.

            nBB = 0

            nFrontLeaves = 1
            frontLeaves(1) = 1

            treeTraversalLoop: do

                ! Initialize the number of leaves for the new front, i.e.
                ! the front of the next round, to 0.

                nFrontLeavesNew = 0

                ! Loop over the leaves of the current front.

                currentFrontLoop: do ii = 1, nFrontLeaves

                    ! Store the ID of the leaf a bit easier and loop over
                    ! its two children.

                    ll = frontLeaves(ii)

                    childrenLoop: do mm = 1, 2

                        ! Determine whether this child contains a bounding box
                        ! or a leaf of the next level.

                        kk = ADTree(ll)%children(mm)
                        terminalTest: if (kk < 0) then

                            ! Child contains a bounding box. Check if the
                            ! coordinate is inside the bounding box.

                            kk = -kk
                            if (coor(1, nn) >= xBBox(1, kk) .and. &
                                coor(1, nn) <= xBBox(4, kk) .and. &
                                coor(2, nn) >= xBBox(2, kk) .and. &
                                coor(2, nn) <= xBBox(5, kk) .and. &
                                coor(3, nn) >= xBBox(3, kk) .and. &
                                coor(3, nn) <= xBBox(6, kk)) then

                                ! Coordinate is inside the bounding box. Store the
                                ! bounding box in the list of possible candidates.

                                if (nBB == nAllocBB) &
                                    call reallocPlus(BB, nAllocBB, 10, jj)

                                nBB = nBB + 1
                                BB(nBB) = kk
                            end if

                            else terminalTest

                            ! Child contains a leaf. Check if the coordinate is
                            ! inside the bounding box of the leaf.

                            if (coor(1, nn) >= ADTree(kk)%xMin(1) .and. &
                                coor(1, nn) <= ADTree(kk)%xMax(4) .and. &
                                coor(2, nn) >= ADTree(kk)%xMin(2) .and. &
                                coor(2, nn) <= ADTree(kk)%xMax(5) .and. &
                                coor(3, nn) >= ADTree(kk)%xMin(3) .and. &
                                coor(3, nn) <= ADTree(kk)%xMax(6)) then

                                ! Coordinate is inside the leaf. Store the leaf in
                                ! the list for the new front.

                                if (nFrontLeavesNew == nAllocFront) then
                                    i = nAllocFront
                                    call reallocPlus(frontLeavesNew, i, 25, jj)
                                    call reallocPlus(frontLeaves, nAllocFront, 25, jj)
                                end if

                                nFrontLeavesNew = nFrontLeavesNew + 1
                                frontLeavesNew(nFrontLeavesNew) = kk

                            end if

                        end if terminalTest

                    end do childrenLoop

                end do currentFrontLoop

                ! End of the loop over the current front. If the new front
                ! is empty the entire tree has been traversed and an exit is
                ! made from the corresponding loop.

                if (nFrontLeavesNew == 0) exit treeTraversalLoop

                ! Copy the data of the new front leaves into the current
                ! front for the next round.

                nFrontLeaves = nFrontLeavesNew
                do ll = 1, nFrontLeaves
                    frontLeaves(ll) = frontLeavesNew(ll)
                end do

            end do treeTraversalLoop
!
!         **************************************************************
!         *                                                            *
!         * Part 2: Loop over the selected bounding boxes and check if *
!         *         the corresponding elements contain the point.      *
!         *                                                            *
!         **************************************************************
!
            elementFound = .false.

            BBoxLoop: do mm = 1, nBB

                ! Determine the element type stored in this bounding box.

                kk = BB(mm)
                select case (ADTs(jj)%elementType(kk))

                case (adtTetrahedron)

                    ! Element is a tetrahedron.
                    ! Compute the coordinates relative to node 1.

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%tetraConn(1, ll)

                    do i = 2, 4
                        n(i) = ADTs(jj)%tetraConn(i, ll)

                        xn(1, i) = ADTs(jj)%coor(1, n(i)) - ADTs(jj)%coor(1, n(1))
                        xn(2, i) = ADTs(jj)%coor(2, n(i)) - ADTs(jj)%coor(2, n(1))
                        xn(3, i) = ADTs(jj)%coor(3, n(i)) - ADTs(jj)%coor(3, n(1))
                    end do

                    x(1) = coor(1, nn) - ADTs(jj)%coor(1, n(1))
                    x(2) = coor(2, nn) - ADTs(jj)%coor(2, n(1))
                    x(3) = coor(3, nn) - ADTs(jj)%coor(3, n(1))

                    ! Determine the matrix for the linear transformation
                    ! from the standard element to the current element.

                    a11 = xn(1, 2); a12 = xn(1, 3); a13 = xn(1, 4)
                    a21 = xn(2, 2); a22 = xn(2, 3); a23 = xn(2, 4)
                    a31 = xn(3, 2); a32 = xn(3, 3); a33 = xn(3, 4)

                    ! Compute the determinant. Make sure that it is not zero
                    ! and invert the value.

                    val = a11 * (a22 * a33 - a32 * a23) + a21 * (a13 * a32 - a12 * a33) &
                          + a31 * (a12 * a23 - a13 * a22)
                    val = sign(adtOne, val) / max(abs(val), adtEps)

                    ! Compute the u, v, w weights for the given coordinate.

                    u = val * ((a22 * a33 - a23 * a32) * x(1) &
                               + (a13 * a32 - a12 * a33) * x(2) &
                               + (a12 * a23 - a13 * a22) * x(3))
                    v = val * ((a23 * a31 - a21 * a33) * x(1) &
                               + (a11 * a33 - a13 * a31) * x(2) &
                               + (a13 * a21 - a11 * a23) * x(3))
                    w = val * ((a21 * a32 - a22 * a31) * x(1) &
                               + (a12 * a31 - a11 * a32) * x(2) &
                               + (a11 * a22 - a12 * a21) * x(3))

                    ! Check if the coordinate is inside the tetrahedron.
                    ! If so, set elementFound to .true. and determine the
                    ! interpolation weights.

                    if (u >= adtZero .and. v >= adtZero .and. &
                        w >= adtZero .and. (u + v + w) <= adtOne) then
                        elementFound = .true.

                        ! Set the number of interpolation nodes to 4 and
                        ! determine the interpolation weights.

                        nNodeElement = 4

                        weight(1) = adtOne - u - v - w
                        weight(2) = u
                        weight(3) = v
                        weight(4) = w
                    end if

                    !=========================================================

                case (adtPyramid)

                    ! Element is a pyramid.
                    ! Compute the coordinates relative to node 1.

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%pyraConn(1, ll)

                    do i = 2, 5
                        n(i) = ADTs(jj)%pyraConn(i, ll)

                        xn(1, i) = ADTs(jj)%coor(1, n(i)) - ADTs(jj)%coor(1, n(1))
                        xn(2, i) = ADTs(jj)%coor(2, n(i)) - ADTs(jj)%coor(2, n(1))
                        xn(3, i) = ADTs(jj)%coor(3, n(i)) - ADTs(jj)%coor(3, n(1))
                    end do

                    x(1) = coor(1, nn) - ADTs(jj)%coor(1, n(1))
                    x(2) = coor(2, nn) - ADTs(jj)%coor(2, n(1))
                    x(3) = coor(3, nn) - ADTs(jj)%coor(3, n(1))

                    ! Modify the coordinates of node 3, such that it
                    ! corresponds to the weights of the u*v term in the
                    ! transformation.

                    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
                    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
                    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)

                    ! Set the starting values of u, v and w such that it is
                    ! somewhere in the middle of the element. In this way the
                    ! Jacobian matrix is always regular, even if the element
                    ! is degenerate.

                    u = adtHalf; v = adtHalf; w = adtHalf

                    ! The Newton algorithm to determine the parametric
                    ! weights u, v and w for the given coordinate.

                    NewtonPyra: do ll = 1, iterMax

                        ! Compute the RHS.

                        uv = u * v
                        oneMinusW = adtOne - w

                        f(1) = oneMinusW * (xn(1, 2) * u + xn(1, 4) * v + xn(1, 3) * uv) &
                               + xn(1, 5) * w - x(1)
                        f(2) = oneMinusW * (xn(2, 2) * u + xn(2, 4) * v + xn(2, 3) * uv) &
                               + xn(2, 5) * w - x(2)
                        f(3) = oneMinusW * (xn(3, 2) * u + xn(3, 4) * v + xn(3, 3) * uv) &
                               + xn(3, 5) * w - x(3)

                        ! Compute the Jacobian.

                        a11 = oneMinusW * (xn(1, 2) + xn(1, 3) * v)
                        a12 = oneMinusW * (xn(1, 4) + xn(1, 3) * u)
                        a13 = xn(1, 5) - xn(1, 2) * u - xn(1, 4) * v - xn(1, 3) * uv

                        a21 = oneMinusW * (xn(2, 2) + xn(2, 3) * v)
                        a22 = oneMinusW * (xn(2, 4) + xn(2, 3) * u)
                        a23 = xn(2, 5) - xn(2, 2) * u - xn(2, 4) * v - xn(2, 3) * uv

                        a31 = oneMinusW * (xn(3, 2) + xn(3, 3) * v)
                        a32 = oneMinusW * (xn(3, 4) + xn(3, 3) * u)
                        a33 = xn(3, 5) - xn(3, 2) * u - xn(3, 4) * v - xn(3, 3) * uv

                        ! Compute the determinant. Make sure that it is not zero
                        ! and invert the value. The cut off is needed to be able
                        ! to handle exceptional cases for degenerate elements.

                        val = a11 * (a22 * a33 - a32 * a23) + a21 * (a13 * a32 - a12 * a33) &
                              + a31 * (a12 * a23 - a13 * a22)
                        val = sign(adtOne, val) / max(abs(val), adtEps)

                        ! Compute the new values of u, v and w.

                        du = val * ((a22 * a33 - a23 * a32) * f(1) &
                                    + (a13 * a32 - a12 * a33) * f(2) &
                                    + (a12 * a23 - a13 * a22) * f(3))
                        dv = val * ((a23 * a31 - a21 * a33) * f(1) &
                                    + (a11 * a33 - a13 * a31) * f(2) &
                                    + (a13 * a21 - a11 * a23) * f(3))
                        dw = val * ((a21 * a32 - a22 * a31) * f(1) &
                                    + (a12 * a31 - a11 * a32) * f(2) &
                                    + (a11 * a22 - a12 * a21) * f(3))

                        u = u - du; v = v - dv; w = w - dw

                        ! Exit the loop if the update of the parametric
                        ! weights is below the threshold

                        val = sqrt(du * du + dv * dv + dw * dw)
                        if (val <= thresConv) exit NewtonPyra
                    end do NewtonPyra

                    ! Check if the coordinate is inside the pyramid.
                    ! If so, set elementFound to .true. and determine the
                    ! interpolation weights.

                    if (u >= adtZero .and. v >= adtZero .and. &
                        w >= adtZero .and. (u + w) <= adtOne .and. &
                        (v + w) <= adtOne) then
                        elementFound = .true.

                        ! Set the number of interpolation nodes to 5 and
                        ! determine the interpolation weights.

                        nNodeElement = 5

                        oneMinusU = adtOne - u
                        oneMinusV = adtOne - v
                        oneMinusW = adtOne - w

                        weight(1) = oneMinusU * oneMinusV * oneMinusW
                        weight(2) = u * oneMinusV * oneMinusW
                        weight(3) = u * v * oneMinusW
                        weight(4) = oneMinusU * v * oneMinusW
                        weight(5) = w
                    end if

                    !=========================================================

                case (adtPrism)

                    ! Element is a prism.
                    ! Compute the coordinates relative to node 1.

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%prismsConn(1, ll)

                    do i = 2, 6
                        n(i) = ADTs(jj)%prismsConn(i, ll)

                        xn(1, i) = ADTs(jj)%coor(1, n(i)) - ADTs(jj)%coor(1, n(1))
                        xn(2, i) = ADTs(jj)%coor(2, n(i)) - ADTs(jj)%coor(2, n(1))
                        xn(3, i) = ADTs(jj)%coor(3, n(i)) - ADTs(jj)%coor(3, n(1))
                    end do

                    x(1) = coor(1, nn) - ADTs(jj)%coor(1, n(1))
                    x(2) = coor(2, nn) - ADTs(jj)%coor(2, n(1))
                    x(3) = coor(3, nn) - ADTs(jj)%coor(3, n(1))

                    ! Modify the coordinates of node 5 and 6, such that they
                    ! correspond to the weights of the u*w and v*w term in the
                    ! transformation respectively.

                    xn(1, 5) = xn(1, 5) - xn(1, 2) - xn(1, 4)
                    xn(2, 5) = xn(2, 5) - xn(2, 2) - xn(2, 4)
                    xn(3, 5) = xn(3, 5) - xn(3, 2) - xn(3, 4)

                    xn(1, 6) = xn(1, 6) - xn(1, 3) - xn(1, 4)
                    xn(2, 6) = xn(2, 6) - xn(2, 3) - xn(2, 4)
                    xn(3, 6) = xn(3, 6) - xn(3, 3) - xn(3, 4)

                    ! Set the starting values of u, v and w such that it is
                    ! somewhere in the middle of the element. In this way the
                    ! Jacobian matrix is always regular, even if the element
                    ! is degenerate.

                    u = adtFourth; v = adtFourth; w = adtHalf

                    ! The Newton algorithm to determine the parametric
                    ! weights u, v and w for the given coordinate.

                    NewtonPrisms: do ll = 1, iterMax

                        ! Compute the RHS.

                        uw = u * w; vw = v * w

                        f(1) = xn(1, 2) * u + xn(1, 3) * v + xn(1, 4) * w &
                               + xn(1, 5) * uw + xn(1, 6) * vw - x(1)
                        f(2) = xn(2, 2) * u + xn(2, 3) * v + xn(2, 4) * w &
                               + xn(2, 5) * uw + xn(2, 6) * vw - x(2)
                        f(3) = xn(3, 2) * u + xn(3, 3) * v + xn(3, 4) * w &
                               + xn(3, 5) * uw + xn(3, 6) * vw - x(3)

                        ! Compute the Jacobian.

                        a11 = xn(1, 2) + xn(1, 5) * w
                        a12 = xn(1, 3) + xn(1, 6) * w
                        a13 = xn(1, 4) + xn(1, 5) * u + xn(1, 6) * v

                        a21 = xn(2, 2) + xn(2, 5) * w
                        a22 = xn(2, 3) + xn(2, 6) * w
                        a23 = xn(2, 4) + xn(2, 5) * u + xn(2, 6) * v

                        a31 = xn(3, 2) + xn(3, 5) * w
                        a32 = xn(3, 3) + xn(3, 6) * w
                        a33 = xn(3, 4) + xn(3, 5) * u + xn(3, 6) * v

                        ! Compute the determinant. Make sure that it is not zero
                        ! and invert the value. The cut off is needed to be able
                        ! to handle exceptional cases for degenerate elements.

                        val = a11 * (a22 * a33 - a32 * a23) + a21 * (a13 * a32 - a12 * a33) &
                              + a31 * (a12 * a23 - a13 * a22)
                        val = sign(adtOne, val) / max(abs(val), adtEps)

                        ! Compute the new values of u, v and w.

                        du = val * ((a22 * a33 - a23 * a32) * f(1) &
                                    + (a13 * a32 - a12 * a33) * f(2) &
                                    + (a12 * a23 - a13 * a22) * f(3))
                        dv = val * ((a23 * a31 - a21 * a33) * f(1) &
                                    + (a11 * a33 - a13 * a31) * f(2) &
                                    + (a13 * a21 - a11 * a23) * f(3))
                        dw = val * ((a21 * a32 - a22 * a31) * f(1) &
                                    + (a12 * a31 - a11 * a32) * f(2) &
                                    + (a11 * a22 - a12 * a21) * f(3))

                        u = u - du; v = v - dv; w = w - dw

                        ! Exit the loop if the update of the parametric
                        ! weights is below the threshold

                        val = sqrt(du * du + dv * dv + dw * dw)
                        if (val <= thresConv) exit NewtonPrisms
                    end do NewtonPrisms

                    ! Check if the coordinate is inside the prism.
                    ! If so, set elementFound to .true. and determine the
                    ! interpolation weights.

                    if (u >= adtZero .and. v >= adtZero .and. &
                        w >= adtZero .and. w <= adtOne .and. &
                        (u + v) <= adtOne) then
                        elementFound = .true.

                        ! Set the number of interpolation nodes to 6 and
                        ! determine the interpolation weights.

                        nNodeElement = 6

                        oneMinusUminusV = adtOne - u - v
                        oneMinusW = adtOne - w

                        weight(1) = oneMinusUminusV * oneMinusW
                        weight(2) = u * oneMinusW
                        weight(3) = v * oneMinusW
                        weight(4) = oneMinusUminusV * w
                        weight(5) = u * w
                        weight(6) = v * w
                    end if

                    !=========================================================

                case (adtHexahedron)

                    ! Element is a hexahedron.
                    ! Compute the coordinates relative to node 1.

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%hexaConn(1, ll)

                    do i = 2, 8
                        n(i) = ADTs(jj)%hexaConn(i, ll)

                        xn(1, i) = ADTs(jj)%coor(1, n(i)) - ADTs(jj)%coor(1, n(1))
                        xn(2, i) = ADTs(jj)%coor(2, n(i)) - ADTs(jj)%coor(2, n(1))
                        xn(3, i) = ADTs(jj)%coor(3, n(i)) - ADTs(jj)%coor(3, n(1))
                    end do

                    x(1) = coor(1, nn) - ADTs(jj)%coor(1, n(1))
                    x(2) = coor(2, nn) - ADTs(jj)%coor(2, n(1))
                    x(3) = coor(3, nn) - ADTs(jj)%coor(3, n(1))

                    ! Modify the coordinates of node 3, 6, 8 and 7 such that
                    ! they correspond to the weights of the u*v, u*w, v*w and
                    ! u*v*w term in the transformation respectively.

                    xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) &
                               - xn(1, 3) - xn(1, 6) - xn(1, 8)
                    xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) &
                               - xn(2, 3) - xn(2, 6) - xn(2, 8)
                    xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) &
                               - xn(3, 3) - xn(3, 6) - xn(3, 8)

                    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
                    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
                    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)

                    xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
                    xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
                    xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)

                    xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
                    xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
                    xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)

                    ! Set the starting values of u, v and w such that it is
                    ! somewhere in the middle of the element. In this way the
                    ! Jacobian matrix is always regular, even if the element
                    ! is degenerate.

                    u = adtHalf; v = adtHalf; w = adtHalf

                    ! The Newton algorithm to determine the parametric
                    ! weights u, v and w for the given coordinate.

                    NewtonHexa: do ll = 1, iterMax

                        ! Compute the RHS.

                        uv = u * v; uw = u * w; vw = v * w; wvu = u * v * w

                        f(1) = xn(1, 2) * u + xn(1, 4) * v + xn(1, 5) * w &
                               + xn(1, 3) * uv + xn(1, 6) * uw + xn(1, 8) * vw &
                               + xn(1, 7) * wvu - x(1)
                        f(2) = xn(2, 2) * u + xn(2, 4) * v + xn(2, 5) * w &
                               + xn(2, 3) * uv + xn(2, 6) * uw + xn(2, 8) * vw &
                               + xn(2, 7) * wvu - x(2)
                        f(3) = xn(3, 2) * u + xn(3, 4) * v + xn(3, 5) * w &
                               + xn(3, 3) * uv + xn(3, 6) * uw + xn(3, 8) * vw &
                               + xn(3, 7) * wvu - x(3)

                        ! Compute the Jacobian.

                        a11 = xn(1, 2) + xn(1, 3) * v + xn(1, 6) * w + xn(1, 7) * vw
                        a12 = xn(1, 4) + xn(1, 3) * u + xn(1, 8) * w + xn(1, 7) * uw
                        a13 = xn(1, 5) + xn(1, 6) * u + xn(1, 8) * v + xn(1, 7) * uv

                        a21 = xn(2, 2) + xn(2, 3) * v + xn(2, 6) * w + xn(2, 7) * vw
                        a22 = xn(2, 4) + xn(2, 3) * u + xn(2, 8) * w + xn(2, 7) * uw
                        a23 = xn(2, 5) + xn(2, 6) * u + xn(2, 8) * v + xn(2, 7) * uv

                        a31 = xn(3, 2) + xn(3, 3) * v + xn(3, 6) * w + xn(3, 7) * vw
                        a32 = xn(3, 4) + xn(3, 3) * u + xn(3, 8) * w + xn(3, 7) * uw
                        a33 = xn(3, 5) + xn(3, 6) * u + xn(3, 8) * v + xn(3, 7) * uv

                        ! Compute the determinant. Make sure that it is not zero
                        ! and invert the value. The cut off is needed to be able
                        ! to handle exceptional cases for degenerate elements.

                        val = a11 * (a22 * a33 - a32 * a23) + a21 * (a13 * a32 - a12 * a33) &
                              + a31 * (a12 * a23 - a13 * a22)
                        val = sign(adtOne, val) / max(abs(val), adtEps)

                        ! Compute the new values of u, v and w.

                        du = val * ((a22 * a33 - a23 * a32) * f(1) &
                                    + (a13 * a32 - a12 * a33) * f(2) &
                                    + (a12 * a23 - a13 * a22) * f(3))
                        dv = val * ((a23 * a31 - a21 * a33) * f(1) &
                                    + (a11 * a33 - a13 * a31) * f(2) &
                                    + (a13 * a21 - a11 * a23) * f(3))
                        dw = val * ((a21 * a32 - a22 * a31) * f(1) &
                                    + (a12 * a31 - a11 * a32) * f(2) &
                                    + (a11 * a22 - a12 * a21) * f(3))

                        u = u - du; v = v - dv; w = w - dw

                        ! Exit the loop if the update of the parametric
                        ! weights is below the threshold

                        val = sqrt(du * du + dv * dv + dw * dw)
                        if (val <= thresConv) exit NewtonHexa
                    end do NewtonHexa

                    ! Check if the coordinate is inside the hexahedron.
                    ! If so, set elementFound to .true. and determine the
                    ! interpolation weights.

                    if (u >= adtZero .and. u <= adtOne .and. &
                        v >= adtZero .and. v <= adtOne .and. &
                        w >= adtZero .and. w <= adtOne) then
                        elementFound = .true.

                        ! Set the number of interpolation nodes to 8 and
                        ! determine the interpolation weights.

                        nNodeElement = 8

                        oneMinusU = adtOne - u
                        oneMinusV = adtOne - v
                        oneMinusW = adtOne - w

                        weight(1) = oneMinusU * oneMinusV * oneMinusW
                        weight(2) = u * oneMinusV * oneMinusW
                        weight(3) = u * v * oneMinusW
                        weight(4) = oneMinusU * v * oneMinusW
                        weight(5) = oneMinusU * oneMinusV * w
                        weight(6) = u * oneMinusV * w
                        weight(7) = u * v * w
                        weight(8) = oneMinusU * v * w
                    end if

                end select

                ! If the coordinate is inside the element store all the
                ! necessary information and exit the loop over the target
                ! bounding boxes.

                if (elementFound) then

                    ! The processor, element type and local element ID.

                    intInfo(1, nn) = ADTs(jj)%myID
                    intInfo(2, nn) = ADTs(jj)%elementType(kk)
                    intInfo(3, nn) = ADTs(jj)%elementID(kk)

                    ! The parametric weights.

                    uvw(1, nn) = u
                    uvw(2, nn) = v
                    uvw(3, nn) = w

                    ! The interpolated solution.

                    do ll = 1, nInterpol
                        ii = 3 + ll
                        uvw(ii, nn) = weight(1) * arrDonor(ll, n(1))
                        do i = 2, nNodeElement
                            uvw(ii, nn) = uvw(ii, nn) + weight(i) * arrDonor(ll, n(i))
                        end do
                    end do

                    ! And exit the loop over the bounding boxes.

                    exit BBoxLoop
                end if

            end do BBoxLoop

        end do coorLoop

        ! Release the memory allocated in this routine.

        deallocate (BB, frontLeaves, frontLeavesNew, stat=ierr)
        if (ierr /= 0) &
          call adtTerminate(jj, "containmentTreeSearch", &
                            "Deallocation failure for BB, &
                            &frontLeaves and frontLeavesNew.")

    end subroutine containmentTreeSearch

    !***************************************************************
    !***************************************************************

    subroutine minDistanceTreeSearch(jj, coor, &
                                     intInfo, uvw, &
                                     allxfs, arrDonor, &
                                     nCoor, nInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the actual minimum distance search in  *
!       * the local tree. It is a local routine in the sense that no   *
!       * communication is involved.                                   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * jj:        Entry in the array ADTs, whose ADT must be        *
!       *            searched.                                         *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates and the currently stored minimum  *
!       *            distance squared of these points:                 *
!       *            coor(1,;): Coordinate 1.                          *
!       *            coor(2,;): Coordinate 2.                          *
!       *            coor(3,;): Coordinate 3.                          *
!       *            coor(4,;): The currently stored minimum distance  *
!       *            squared.                                          *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 2D integer array, in which the following output     *
!       *          will be stored:                                     *
!       *          intInfo(1,:): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2,:): The element type of the element.      *
!       *          intInfo(3,:): The element ID of the element in the  *
!       *                        the connectivity.                     *
!       * uvw:     2D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          as well as the new distance squared and the         *
!       *          interpolated solution:                              *
!       *          uvw(1, :): Parametric u-weight.                     *
!       *          uvw(2, :): Parametric v-weight.                     *
!       *          uvw(3, :): Parametric w-weight.                     *
!       *          uvw(4, :): The new distance squared.                *
!       *          uvw(5:,:): Interpolated solution, if desired. It is *
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, jj
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(:, :), intent(in) :: coor
        real(kind=realType), dimension(:, :), intent(in) :: arrDonor

        integer(kind=intType), dimension(:, :), intent(out) :: intInfo
        real(kind=realType), dimension(:, :), intent(out) :: uvw
        real(kind=realType), dimension(:, :), intent(out) :: allxfs

!
!       Local parameters used in the Newton algorithm.
!
        integer(kind=intType), parameter :: iterMax = 15
        real(kind=realType), parameter :: adtEps = 1.e-25_realType
        real(kind=realType), parameter :: thresConv = 1.e-10_realType
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn, activeLeaf
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront, nNodeElement
        integer(kind=intType) :: i, kkk

        integer(kind=intType), dimension(8) :: n, m

        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew

        real(kind=realType) :: dx, dy, dz, d1, d2, val
        real(kind=realType) :: u, v, w
        real(kind=realType) :: uu, vv, ww

        real(kind=realType), dimension(2) :: dd
        real(kind=realType), dimension(3) :: x1, x2, x3, x4, x, xf
        real(kind=realType), dimension(3) :: vf, vt, a, b, an, bn
        real(kind=realType), dimension(3) :: chi
        real(kind=realType), dimension(3) :: xf_
        real(kind=realType), dimension(8) :: weight, donorData

        real(kind=realType), dimension(:, :), pointer :: xBBox

        logical :: elementFound

        type(adtBBoxTargetType), dimension(:), pointer :: BB
        type(adtLeafType), dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox => ADTs(jj)%xBBox
        ADTree => ADTs(jj)%ADTree

        ! Initial allocation of the arrays for the tree traversal as well
        ! as the stack array used in the qsort routine. The latter is
        ! done, because the qsort routine is called for every coordinate
        ! and therefore it is more efficient to allocate the stack once
        ! rather than over and over again. The disadvantage of course is
        ! that an essentially local variable, stack, is now stored in
        ! adtData.

        nAllocBB = 10
        nAllocFront = 25
        nStack = 100

        allocate (stack(nStack), BB(nAllocBB), frontLeaves(nAllocFront), &
                  frontLeavesNew(nAllocFront), stat=ierr)
        if (ierr /= 0) &
          call adtTerminate(jj, "minDistanceTreeSearch", &
                            "Memory allocation failure for stack, BB, &
                            &etc.")

        ! Loop over the number of coordinates to be treated.

        coorLoop: do nn = 1, nCoor

            ! Initialize the processor ID to -1 to indicate that no
            ! corresponding volume element is found and the new minimum
            ! distance squared to the old value.

            intInfo(1, nn) = -1
            uvw(4, nn) = coor(4, nn)
!
!         **************************************************************
!         *                                                            *
!         * Part 1. Determine the possible minimum distance squared to *
!         *         the root leaf. If larger than the current distance *
!         *         there is no need to search this tree.              *
!         *                                                            *
!         **************************************************************
!
            if (coor(1, nn) < ADTree(1)%xMin(1)) then
                dx = coor(1, nn) - ADTree(1)%xMin(1)
            else if (coor(1, nn) > ADTree(1)%xMax(4)) then
                dx = coor(1, nn) - ADTree(1)%xMax(4)
            else
                dx = adtZero
            end if

            if (coor(2, nn) < ADTree(1)%xMin(2)) then
                dy = coor(2, nn) - ADTree(1)%xMin(2)
            else if (coor(2, nn) > ADTree(1)%xMax(5)) then
                dy = coor(2, nn) - ADTree(1)%xMax(5)
            else
                dy = adtZero
            end if

            if (coor(3, nn) < ADTree(1)%xMin(3)) then
                dz = coor(3, nn) - ADTree(1)%xMin(3)
            else if (coor(3, nn) > ADTree(1)%xMax(6)) then
                dz = coor(3, nn) - ADTree(1)%xMax(6)
            else
                dz = adtZero
            end if

            ! Continue with the next coordinate if the possible distance
            ! squared to the root leaf is larger than the currently stored
            ! value.

            if ((dx * dx + dy * dy + dz * dz) >= uvw(4, nn)) cycle
!
!         **************************************************************
!         *                                                            *
!         * Part 2. Find a likely bounding box, which minimizes the    *
!         *         guaranteed distance.                               *
!         *                                                            *
!         **************************************************************
!
            activeLeaf = 1

            ! Traverse the tree until a terminal leaf is found.

            treeTraversal1: do

                ! Exit the loop when a terminal leaf has been found.
                ! This is indicated by a negative value of activeLeaf.

                if (activeLeaf < 0) exit treeTraversal1

                ! Determine the guaranteed distance squared for both children.

                do mm = 1, 2

                    ! Determine whether the child contains a bounding box or
                    ! another leaf of the tree.

                    ll = ADTree(activeLeaf)%children(mm)
                    if (ll > 0) then

                        ! Child contains a leaf. Determine the guaranteed distance
                        ! vector to the leaf.

                        d1 = abs(coor(1, nn) - ADTree(ll)%xMin(1))
                        d2 = abs(coor(1, nn) - ADTree(ll)%xMax(4))
                        dx = max(d1, d2)

                        d1 = abs(coor(2, nn) - ADTree(ll)%xMin(2))
                        d2 = abs(coor(2, nn) - ADTree(ll)%xMax(5))
                        dy = max(d1, d2)

                        d1 = abs(coor(3, nn) - ADTree(ll)%xMin(3))
                        d2 = abs(coor(3, nn) - ADTree(ll)%xMax(6))
                        dz = max(d1, d2)

                    else

                        ! Child contains a bounding box. Determine the guaranteed
                        ! distance vector to it.

                        ll = -ll

                        d1 = abs(coor(1, nn) - xBBox(1, ll))
                        d2 = abs(coor(1, nn) - xBBox(4, ll))
                        dx = max(d1, d2)

                        d1 = abs(coor(2, nn) - xBBox(2, ll))
                        d2 = abs(coor(2, nn) - xBBox(5, ll))
                        dy = max(d1, d2)

                        d1 = abs(coor(3, nn) - xBBox(3, ll))
                        d2 = abs(coor(3, nn) - xBBox(6, ll))
                        dz = max(d1, d2)

                    end if

                    ! Compute the guaranteed distance squared for this child.

                    dd(mm) = dx * dx + dy * dy + dz * dz

                end do

                ! Determine which will be the next leaf in the tree traversal.
                ! This will be the leaf which has the minimum guaranteed
                ! distance. In case of ties take the left leaf, because this
                ! leaf may have more children.

                if (dd(1) <= dd(2)) then
                    activeLeaf = ADTree(activeLeaf)%children(1)
                else
                    activeLeaf = ADTree(activeLeaf)%children(2)
                end if

            end do treeTraversal1

            ! Store the minimum of the just computed guaranteed distance
            ! squared and the currently stored value in uvw.

            uvw(4, nn) = min(uvw(4, nn), dd(1), dd(2))
!
!         **************************************************************
!         *                                                            *
!         * Part 3. Find the bounding boxes whose possible minimum     *
!         *         distance is less than the currently stored value.  *
!         *                                                            *
!         **************************************************************
!
            ! In part 1 it was already tested that the possible distance
            ! squared of the root leaf was less than the current value.
            ! Therefore initialize the current front to the root leaf and
            ! set the number of bounding boxes to 0.

            nBB = 0

            nFrontLeaves = 1
            frontLeaves(1) = 1

            ! Second tree traversal. Now to find all possible bounding
            ! box candidates.

            treeTraversal2: do

                ! Initialize the number of leaves for the new front, i.e.
                ! the front of the next round, to 0.

                nFrontLeavesNew = 0

                ! Loop over the leaves of the current front.

                currentFrontLoop: do ii = 1, nFrontLeaves

                    ! Store the ID of the leaf a bit easier and loop over
                    ! its two children.

                    ll = frontLeaves(ii)

                    childrenLoop: do mm = 1, 2

                        ! Determine whether this child contains a bounding box
                        ! or a leaf of the next level.

                        kk = ADTree(ll)%children(mm)
                        terminalTest: if (kk < 0) then

                            ! Child contains a bounding box. Determine the possible
                            ! minimum distance squared to this bounding box.

                            kk = -kk

                            if (coor(1, nn) < xBBox(1, kk)) then
                                dx = coor(1, nn) - xBBox(1, kk)
                            else if (coor(1, nn) > xBBox(4, kk)) then
                                dx = coor(1, nn) - xBBox(4, kk)
                            else
                                dx = adtZero
                            end if

                            if (coor(2, nn) < xBBox(2, kk)) then
                                dy = coor(2, nn) - xBBox(2, kk)
                            else if (coor(2, nn) > xBBox(5, kk)) then
                                dy = coor(2, nn) - xBBox(5, kk)
                            else
                                dy = adtZero
                            end if

                            if (coor(3, nn) < xBBox(3, kk)) then
                                dz = coor(3, nn) - xBBox(3, kk)
                            else if (coor(3, nn) > xBBox(6, kk)) then
                                dz = coor(3, nn) - xBBox(6, kk)
                            else
                                dz = adtZero
                            end if

                            d2 = dx * dx + dy * dy + dz * dz

                            ! If this distance squared is less than the current
                            ! value, store this bounding box as a target.

                            testStoreBBox: if (d2 < uvw(4, nn)) then

                                ! Check if the memory must be reallocated.

                                if (nBB == nAllocBB) &
                                    call reallocBBoxTargetTypePlus(BB, nAllocBB, &
                                                                   10, jj)

                                ! Update the counter and store the data.

                                nBB = nBB + 1
                                BB(nBB)%ID = kk
                                BB(nBB)%posDist2 = d2

                                ! Although in step 2, i.e. the first tree traversal,
                                ! the guaranteed distance squared to a bounding box
                                ! has already been computed, this has been done only
                                ! for a likely candidate and not for all the possible
                                ! candidates. As this test is relatively cheap, do it
                                ! now for this bounding box.

                                d1 = abs(coor(1, nn) - xBBox(1, kk))
                                d2 = abs(coor(1, nn) - xBBox(4, kk))
                                dx = max(d1, d2)

                                d1 = abs(coor(2, nn) - xBBox(2, kk))
                                d2 = abs(coor(2, nn) - xBBox(5, kk))
                                dy = max(d1, d2)

                                d1 = abs(coor(3, nn) - xBBox(3, kk))
                                d2 = abs(coor(3, nn) - xBBox(6, kk))
                                dz = max(d1, d2)

                                d2 = dx * dx + dy * dy + dz * dz
                                uvw(4, nn) = min(uvw(4, nn), d2)

                            end if testStoreBBox

                            else terminalTest

                            ! Child contains a leaf. Compute the possible minimum
                            ! distance squared to the current coordinate.

                            if (coor(1, nn) < ADTree(kk)%xMin(1)) then
                                dx = coor(1, nn) - ADTree(kk)%xMin(1)
                            else if (coor(1, nn) > ADTree(kk)%xMax(4)) then
                                dx = coor(1, nn) - ADTree(kk)%xMax(4)
                            else
                                dx = adtZero
                            end if

                            if (coor(2, nn) < ADTree(kk)%xMin(2)) then
                                dy = coor(2, nn) - ADTree(kk)%xMin(2)
                            else if (coor(2, nn) > ADTree(kk)%xMax(5)) then
                                dy = coor(2, nn) - ADTree(kk)%xMax(5)
                            else
                                dy = adtZero
                            end if

                            if (coor(3, nn) < ADTree(kk)%xMin(3)) then
                                dz = coor(3, nn) - ADTree(kk)%xMin(3)
                            else if (coor(3, nn) > ADTree(kk)%xMax(6)) then
                                dz = coor(3, nn) - ADTree(kk)%xMax(6)
                            else
                                dz = adtZero
                            end if

                            d2 = dx * dx + dy * dy + dz * dz

                            ! If this distance squared is less than the current
                            ! value, store this leaf in the new front.

                            testStoreLeave: if (d2 < uvw(4, nn)) then

                                ! Check if enough memory has been allocated and
                                ! store the leaf.

                                if (nFrontLeavesNew == nAllocFront) then
                                    i = nAllocFront
                                    call reallocPlus(frontLeavesNew, i, 25, jj)
                                    call reallocPlus(frontLeaves, nAllocFront, 25, jj)
                                end if

                                nFrontLeavesNew = nFrontLeavesNew + 1
                                frontLeavesNew(nFrontLeavesNew) = kk

                                ! Compute the guaranteed distance squared to this leaf.
                                ! It may be less than the currently stored value.

                                d1 = abs(coor(1, nn) - ADTree(kk)%xMin(1))
                                d2 = abs(coor(1, nn) - ADTree(kk)%xMax(4))
                                dx = max(d1, d2)

                                d1 = abs(coor(2, nn) - ADTree(kk)%xMin(2))
                                d2 = abs(coor(2, nn) - ADTree(kk)%xMax(5))
                                dy = max(d1, d2)

                                d1 = abs(coor(3, nn) - ADTree(kk)%xMin(3))
                                d2 = abs(coor(3, nn) - ADTree(kk)%xMax(6))
                                dz = max(d1, d2)

                                d2 = dx * dx + dy * dy + dz * dz
                                uvw(4, nn) = min(uvw(4, nn), d2)

                            end if testStoreLeave

                        end if terminalTest
                    end do childrenLoop
                end do currentFrontLoop

                ! End of the loop over the current front. If the new front
                ! is empty the entire tree has been traversed and an exit is
                ! made from the corresponding loop.

                if (nFrontLeavesNew == 0) exit treeTraversal2

                ! Copy the data of the new front leaves into the current
                ! front for the next round.

                nFrontLeaves = nFrontLeavesNew
                do ll = 1, nFrontLeaves
                    frontLeaves(ll) = frontLeavesNew(ll)
                end do

            end do treeTraversal2

            ! Sort the target bounding boxes in increasing order such that
            ! the one with the smallest possible distance is first.

            call qsortBBoxTargets(BB, nBB, jj)
!
!         **************************************************************
!         *                                                            *
!         * Part 4: Loop over the selected bounding boxes and check if *
!         *         the corresponding element minimizes the distance.  *
!         *                                                            *
!         **************************************************************
!
            elementFound = .false.

            BBoxLoop: do mm = 1, nBB

                ! Exit the loop if the possible minimum distance of this
                ! bounding box is not smaller than the current value.
                ! Remember that BB has been sorted in increasing order.

                if (uvw(4, nn) <= BB(mm)%posDist2) exit BBoxLoop

                ! Determine the element type stored in this bounding box.

                kk = BB(mm)%ID
                select case (ADTs(jj)%elementType(kk))

                case (adtTriangle)

                    ! =============================================
                    ! =============================================
                    ! =============================================
                    ! =============================================

                    ! This is a surface element, so set the parametric weight
                    ! w to zero.

                    w = adtZero

                    ! Determine the 3 nodes which completely describe
                    ! the traingle face

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%triaConn(1, ll)
                    n(2) = ADTs(jj)%triaConn(2, ll)
                    n(3) = ADTs(jj)%triaConn(3, ll)

                    x1(1) = ADTs(jj)%coor(1, n(1))
                    x1(2) = ADTs(jj)%coor(2, n(1))
                    x1(3) = ADTs(jj)%coor(3, n(1))

                    x2(1) = ADTs(jj)%coor(1, n(2))
                    x2(2) = ADTs(jj)%coor(2, n(2))
                    x2(3) = ADTs(jj)%coor(3, n(2))

                    x3(1) = ADTs(jj)%coor(1, n(3))
                    x3(2) = ADTs(jj)%coor(2, n(3))
                    x3(3) = ADTs(jj)%coor(3, n(3))

                    ! Get the point that should be projected

                    x = coor(1:3, nn)

                    ! Call projection algorithm
                    ! triaProjection defined in adtProjections.F90

                    call triaProjection(x1, x2, x3, x, xf, u, v, val)

                    ! If the distance squared is less than the current value
                    ! store the wall distance and interpolation info and
                    ! indicate that an element was found.

                    if (val < uvw(4, nn)) then
                        uvw(4, nn) = val
                        nNodeElement = 3
                        elementFound = .true.

                        kkk = kk; uu = u; vv = v; ww = w
                        m(1) = n(1); m(2) = n(2); m(3) = n(3)

                        ! triaWeights defined in adtProjections.F90
                        call triaWeights(u, v, weight)

                        xf_(1) = xf(1)
                        xf_(2) = xf(2)
                        xf_(3) = xf(3)
                    end if

                    ! =============================================
                    ! =============================================
                    ! =============================================
                    ! =============================================

                case (adtQuadrilateral)

                    ! Temporary implementation. I'm waiting for Juan to come
                    ! up with his more sophisticated algorithm.

                    ! This is a surface element, so set the parametric weight
                    ! w to zero.

                    w = adtZero

                    ! Determine the 4 vectors which completely describe
                    ! the quadrilateral face

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%quadsConn(1, ll)
                    n(2) = ADTs(jj)%quadsConn(2, ll)
                    n(3) = ADTs(jj)%quadsConn(3, ll)
                    n(4) = ADTs(jj)%quadsConn(4, ll)

                    x1(1) = ADTs(jj)%coor(1, n(1))
                    x1(2) = ADTs(jj)%coor(2, n(1))
                    x1(3) = ADTs(jj)%coor(3, n(1))

                    x2(1) = ADTs(jj)%coor(1, n(2))
                    x2(2) = ADTs(jj)%coor(2, n(2))
                    x2(3) = ADTs(jj)%coor(3, n(2))

                    x3(1) = ADTs(jj)%coor(1, n(3))
                    x3(2) = ADTs(jj)%coor(2, n(3))
                    x3(3) = ADTs(jj)%coor(3, n(3))

                    x4(1) = ADTs(jj)%coor(1, n(4))
                    x4(2) = ADTs(jj)%coor(2, n(4))
                    x4(3) = ADTs(jj)%coor(3, n(4))

                    ! Get the point that should be projected

                    x = coor(1:3, nn)

                    ! Call projection algorithm
                    ! quadProjection defined in adtProjections.F90
                    call quadProjection(x1, x2, x3, x4, &
                                        x, &
                                        xf, u, v, val)

                    ! If the distance squared is less than the current value
                    ! store the wall distance and interpolation info and
                    ! indicate that an element was found.

                    if (val < uvw(4, nn)) then
                        uvw(4, nn) = val
                        nNodeElement = 4
                        elementFound = .true.

                        kkk = kk; uu = u; vv = v; ww = w
                        m(1) = n(1); m(2) = n(2); m(3) = n(3); m(4) = n(4)

                        ! quadWeights defined in adtProjections.F90
                        call quadWeights(u, v, weight)

                        xf_(1) = xf(1)
                        xf_(2) = xf(2)
                        xf_(3) = xf(3)

                    end if

                    !=========================================================

                case (adtTetrahedron)
                    call adtTerminate(jj, "minDistanceTreeSearch", &
                                      "Minimum distance search for &
                                      &tetrahedra not implemented yet")

                    !===========================================================

                case (adtPyramid)
                    call adtTerminate(jj, "minDistanceTreeSearch", &
                                      "Minimum distance search for &
                                      &pyramids not implemented yet")

                    !===========================================================

                case (adtPrism)
                    call adtTerminate(jj, "minDistanceTreeSearch", &
                                      "Minimum distance search for &
                                      &prisms not implemented yet")

                    !===========================================================

                case (adtHexahedron)

                    ! Determine the element ID and the corresponding
                    ! 8 node ID's.

                    ll = ADTs(jj)%elementID(kk)
                    n(1) = ADTs(jj)%hexaConn(1, ll)
                    n(2) = ADTs(jj)%hexaConn(2, ll)
                    n(3) = ADTs(jj)%hexaConn(3, ll)
                    n(4) = ADTs(jj)%hexaConn(4, ll)
                    n(5) = ADTs(jj)%hexaConn(5, ll)
                    n(6) = ADTs(jj)%hexaConn(6, ll)
                    n(7) = ADTs(jj)%hexaConn(7, ll)
                    n(8) = ADTs(jj)%hexaConn(8, ll)

                    ! Call the subroutine minD2Hexa to do the work.

                    call minD2Hexa(coor(1:3, nn), &
                                   ADTs(jj)%coor(1:3, n(1)), &
                                   ADTs(jj)%coor(1:3, n(2)), &
                                   ADTs(jj)%coor(1:3, n(3)), &
                                   ADTs(jj)%coor(1:3, n(4)), &
                                   ADTs(jj)%coor(1:3, n(5)), &
                                   ADTs(jj)%coor(1:3, n(6)), &
                                   ADTs(jj)%coor(1:3, n(7)), &
                                   ADTs(jj)%coor(1:3, n(8)), &
                                   val, chi, ierr)

                    ! If the distance squared is less than the current value
                    ! store the wall distance and interpolation info and
                    ! indicate that an element was found.

                    if (val < uvw(4, nn)) then
                        uvw(4, nn) = val
                        nNodeElement = 8
                        elementFound = .true.

                        kkk = kk; 
                        uu = chi(1); vv = chi(2); ww = chi(3)

                        m(1) = n(1); m(2) = n(2); m(3) = n(3); m(4) = n(4)
                        m(5) = n(5); m(6) = n(6); m(7) = n(7); m(8) = n(8)

                        weight(1) = (adtOne - uu) * (adtOne - vv) * (adtOne - ww)
                        weight(2) = uu * (adtOne - vv) * (adtOne - ww)
                        weight(3) = uu * vv * (adtOne - ww)
                        weight(4) = (adtOne - uu) * vv * (adtOne - ww)
                        weight(5) = (adtOne - uu) * (adtOne - vv) * ww
                        weight(6) = uu * (adtOne - vv) * ww
                        weight(7) = uu * vv * ww
                        weight(8) = (adtOne - uu) * vv * ww
                    end if

                end select

            end do BBoxLoop

            ! Check if an element was found. As all the minimum distance
            ! searches are initialized by the calling routine (to support
            ! periodic searches) this is not always the case.

            if (elementFound) then

                ! Store the interpolation information for this point.
                ! First the integer info, i.e. the processor ID, element type
                ! and local element ID.

                intInfo(1, nn) = ADTs(jj)%myID
                intInfo(2, nn) = ADTs(jj)%elementType(kkk)
                intInfo(3, nn) = ADTs(jj)%elementID(kkk)

                ! The parametric weights. Note that the wall distance
                ! squared, stored in the 4th position of uvw, already
                ! contains the correct value.

                uvw(1, nn) = uu
                uvw(2, nn) = vv
                uvw(3, nn) = ww

                ! The interpolated solution, if needed.

                do ll = 1, nInterpol
                    ii = 4 + ll
                    ! Gather data from nodes
                    donorData = adtZero
                    do i = 1, nNodeElement
                        donorData(i) = arrDonor(ll, m(i))
                    end do
                    ! Dot function defined in utilities.F90
                    call dotProd(weight, donorData, uvw(ii, nn))

                    ! Old code
                    !uvw(ii,nn) = weight(1)*arrDonor(ll,m(1))
                    !do i=2,nNodeElement
                    !  uvw(ii,nn) = uvw(ii,nn) + weight(i)*arrDonor(ll,m(i))
                    !enddo

                end do

                ! Store global coordinate of projected point
                allxfs(1, nn) = xf_(1)
                allxfs(2, nn) = xf_(2)
                allxfs(3, nn) = xf_(3)
            end if

        end do coorLoop

        ! Release the memory allocated in this routine.

        deallocate (stack, BB, frontLeaves, frontLeavesNew, stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(jj, "minDistanceTreeSearch", &
                              "Deallocation failure for stack, BB, etc.")

    end subroutine minDistanceTreeSearch

!===============================================================

    subroutine intersectionTreeSearch(jj, inpBBox, &
                                      intInfo, BBoxPtr, nBBox)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the detects elements whose bounding    *
!       * boxes intersect the user provided bounding boxes (inpBBox).  *
!       * This narrows down, for instance, candidates in intersection  *
!       * algorithms.                                                  *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * jj:        Entry in the array ADTs, whose ADT must be        *
!       *            searched.                                         *
!       * nBBox:     Number of bounding boxes for which the            *
!       *            intersecting elements must be determined.         *
!       * inpBBox:   real(6,nBBox) containing the corners of each      *
!       *            bounding box (xmin, ymin, zmin, xmax, ymax, zmax).*
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 2D integer array, in which the following output     *
!       *          will be stored:                                     *
!       *          intInfo(1,:): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2,:): The element type of the element.      *
!       *          intInfo(3,:): The element ID of the element in the  *
!       *                        connectivity.                         *
!       *          This array will hold all elements that intersects   *
!       *          any of the given bounding boxes. Then you need to   *
!       *          use BBoxPtr to find which slice of the intInfo      *
!       *          array corresponds to a specific bounding box given  *
!       *          in inpBBox.                                         *
!       *                                                              *
!       * BBoxPtr: integer(nBBox+1). Pointers used to slice intInfo.     *
!       *          The indices from BBoxPtr(i) to BBoxPtr(i+1)-1       *
!       *          belong to inpBBox(i).                               *
!       *                                                              *
!       ****************************************************************
!
        use utilities ! ../utilities/utilities.F90
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nBBox, jj

        real(kind=realType), dimension(6, nBBox), intent(in) :: inpBBox

        integer(kind=intType), dimension(:, :), allocatable, intent(out) :: intInfo

        integer(kind=intType), dimension(nBBox + 1), intent(out) :: BBoxPtr
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn, i
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront, nAllocIntInfo
        integer(kind=intType) :: nBBoxTotal, startID, endID

        integer(kind=intType), dimension(:), pointer :: BB
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew

        real(kind=realType), dimension(6) :: BBoxA, BBoxB, BBoxAB

        real(kind=realType), dimension(:, :), pointer :: xBBox

        logical :: overlap

        type(adtLeafType), dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox => ADTs(jj)%xBBox
        ADTree => ADTs(jj)%ADTree

        ! Initial allocation of the arrays for the tree traversal.

        nAllocBB = 10
        nAllocFront = 25
        nAllocIntInfo = 30

        allocate (BB(nAllocBB), frontLeaves(nAllocFront), &
                  frontLeavesNew(nAllocFront), intInfo(3, nAllocIntInfo), &
                  stat=ierr)
        if (ierr /= 0) &
             call adtTerminate(jj, "containmentTreeSearch", &
             "Memory allocation failure for BB, &
             &frontLeaves and frontLeavesNew.")

        ! Initialize the processor ID to -1 to indicate that no
        ! corresponding bounding box is found.
        intInfo = -1

        ! Initialize counter of total number of BBoxes flagged so far
        nBBoxTotal = 0

        ! Loop over the number of bounding boxes to be treated.

        coorLoop: do nn = 1, nBBox

            ! Get corners of the current bounding box
            BBoxA = inpBBox(:, nn)

!
!         **************************************************************
!         *                                                            *
!         * Part 1. Traverse the tree and determine the target         *
!         *         elements, which may intersect the given BBox.      *
!         *                                                            *
!         **************************************************************
!
            ! Start at the root, i.e. set the front leaf to the root leaf.
            ! Also initialize the number of possible bounding boxes to 0.

            nBB = 0

            nFrontLeaves = 1
            frontLeaves(1) = 1

            treeTraversalLoop: do

                ! Initialize the number of leaves for the new front, i.e.
                ! the front of the next round, to 0.

                nFrontLeavesNew = 0

                ! Loop over the leaves of the current front.

                currentFrontLoop: do ii = 1, nFrontLeaves

                    ! Store the ID of the leaf a bit easier and loop over
                    ! its two children.

                    ll = frontLeaves(ii)

                    childrenLoop: do mm = 1, 2

                        ! Determine whether this child contains a bounding box
                        ! or a leaf of the next level.

                        kk = ADTree(ll)%children(mm)

                        terminalTest: if (kk < 0) then ! We have a terminal leaf with a BBox!

                            ! Child contains a bounding box. Check if the
                            ! coordinate is inside the bounding box.

                            kk = -kk ! This is the BBox ID

                            ! Get coordinates of the candidate bounding box
                            BBoxB = xBBox(:, kk)

                            ! Check bounding box intersection
                            ! computeBBoxIntersection defined in utilities.F90 of the utilities foldes
                            call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

                            if (overlap) then

                                ! The bounding boxes intersect. Store the
                                ! bounding box in the list of possible candidates.
                                ! But first check if we need to increase the BBox array
                                if (nBB == nAllocBB) &
                                    call reallocPlus(BB, nAllocBB, 10, jj)

                                ! Store the new candidate BBox
                                nBB = nBB + 1
                                BB(nBB) = kk
                            end if

                            else terminalTest ! We do not have a terminal leaf because kk is positive

                            ! Child contains a leaf. Check if the current BBox intersects the
                            ! leaf bounding box.

                            ! Get coordinates of the leaf BBox
                            BBoxB(1:3) = ADTree(kk)%xMin(1:3)
                            BBoxB(4:6) = ADTree(kk)%xMax(4:6)

                            ! Check bounding box intersection
                            ! computeBBoxIntersection defined in utilities.F90 of the utilities foldes
                            call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

                            if (overlap) then

                                ! Our bounding box intersects the leaf BBox. Store the leaf in
                                ! the list for the new front.

                                if (nFrontLeavesNew == nAllocFront) then
                                    i = nAllocFront
                                    call reallocPlus(frontLeavesNew, i, 25, jj)
                                    call reallocPlus(frontLeaves, nAllocFront, 25, jj)
                                end if

                                nFrontLeavesNew = nFrontLeavesNew + 1
                                frontLeavesNew(nFrontLeavesNew) = kk

                            end if

                        end if terminalTest

                    end do childrenLoop

                end do currentFrontLoop

                ! End of the loop over the current front. If the new front
                ! is empty the entire tree has been traversed and an exit is
                ! made from the corresponding loop.

                if (nFrontLeavesNew == 0) exit treeTraversalLoop

                ! Copy the data of the new front leaves into the current
                ! front for the next round.

                nFrontLeaves = nFrontLeavesNew
                do ll = 1, nFrontLeaves
                    frontLeaves(ll) = frontLeavesNew(ll)
                end do

            end do treeTraversalLoop
!
!         **************************************************************
!         *                                                            *
!         * Part 2: Store the new bounding boxes in the output array   *
!         *                                                            *
!         **************************************************************
!

            do ii = 1, nBB

                ! Get BBox ID that we stored in part 1
                kk = BB(ii)

                ! Check if we need to reallocate intInfo.
                ! reallocPlus2D defined in adtUtils.F90
                if (nBBoxTotal + ii .gt. nAllocIntInfo) then
                    call reallocPlus2D(intInfo, nAllocIntInfo, 30, jj)
                end if

                ! Store proc ID
                intInfo(1, nBBoxTotal + ii) = ADTs(jj)%myID

                ! Store element type
                intInfo(2, nBBoxTotal + ii) = ADTs(jj)%elementType(kk)

                ! Store element ID
                intInfo(3, nBBoxTotal + ii) = ADTs(jj)%elementID(kk)

            end do

            ! Store pointer info for the BBox we are currently checking.
            BBoxPtr(nn) = nBBoxTotal + 1

            ! Update total number of BBox found so far
            nBBoxTotal = nBBoxTotal + nBB

        end do coorLoop

        ! Add the final pointer to BBoxPtr (which will point to the end of ther array)
        BBoxPtr(nBBox + 1) = nBBoxTotal + 1

        ! Now we crop intInfo to have the correct number of bounding boxes
        call reallocPlus2D(intInfo, nAllocIntInfo, nBBoxTotal - nAllocIntInfo, jj)

    end subroutine intersectionTreeSearch

end module adtLocalSearch
