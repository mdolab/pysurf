!
!     ******************************************************************
!     *                                                                *
!     * File:          curveSearch.F90                                 *
!     * Author:        John Jasa                                       *
!     * Starting date: 07-27-2016                                      *
!     * Last modified: 11-08-2016                                      *
!     *                                                                *
!     ******************************************************************
!
module curveSearchAPI

  use precision
  implicit none

  !=================================================================

contains

  subroutine minDistanceCurve(nxyz, nCoor, nBars, xyz, coor, barsConn, &
                              allProjPoints, allTangents, allDist2, allElemIDs, &
                              curveMask)
    !
    !       ****************************************************************
    !       * John Jasa                                                    *
    !       *                                                              *
    !       ****************************************************************
    !

    use curveUtils
    implicit none

    !f2py intent(in) nxyz, nCoor, nBars, xyz, coor, barsConn
    !f2py intent(inout) allProjPoints, allTangents, allDist2, allElemIDs
    !f2py intent(out) curveMask

    !
    !       Subroutine arguments.
    !
    ! Input
    integer(kind=intType), intent(in) :: nxyz, nCoor, nBars
    real(kind=realType), dimension(3,nxyz), intent(in) :: xyz
    real(kind=realType), dimension(3,nCoor), intent(in) :: coor
    integer(kind=intType), dimension(2,nBars), intent(in) :: barsConn

    ! Input/Output
    real(kind=realType), dimension(3,nxyz), intent(inout) :: allProjPoints
    real(kind=realType), dimension(3,nxyz), intent(inout) :: allTangents
    real(kind=realType), dimension(nxyz), intent(inout) :: allDist2
    integer(kind=intType), dimension(nxyz), intent(inout) :: allElemIDs

    ! Outputs
    integer(kind=intType), dimension(nxyz), intent(out) :: curveMask

    ! Working
    real(kind=realType) :: x(3), x21(3), vec(3)
    real(kind=realType) :: x1(3), x2(3), p1(3), u, mag2
    real(kind=realType) :: xf(3), dist2, newDist2, projPoint(3)
    real(kind=realType) :: tangent(3)
    integer(kind=intType) :: i, j, jj

    !===============================================================

    ! Initialize curveMask. This array initially have zeroes. If we find a better projection
    ! for the ii-th node on this curve, we will set curveMask[ii] to 1.
    curveMask = 0

    ! Loop over all points to be projected
    do i = 1,nxyz

       ! Obtain the raw point to be projected and its current minimum
       ! distance from a curve
       x = xyz(:, i)
       dist2 = allDist2(i)

       ! Initialize jj, which should be the index of the element that is
       ! best suited for the projection
       jj = 0

       ! Loop over all individual bar elements that define the curve
       do j = 1,nBars

          ! Get the absolute nodal points from the coordinate array
          x1 = coor(:, barsConn(1, j))
          x2 = coor(:, barsConn(2, j))

          ! Call projection function
          call barProjection(x1,x2,x,xf,u)

          ! Compute the new projected point coordinates in the global frame
          vec = xf - x
          newDist2 = dot_product(vec, vec)

          ! If the distance for this new projected point is less than the
          ! provided distance, save the new point
          if (newDist2 .lt. dist2) then
             dist2 = newDist2
             projPoint = xf
             jj = j
          end if

       end do

       ! If we could not find a better candidate in this curve, then
       ! jj still has its initial value. Otherwise, we found a better candidate
       ! and we need to replace values.
       if (jj .gt. 0) then

          ! Compute tangent direction of the best bar element
          x1 = coor(:, barsConn(1, jj))
          x2 = coor(:, barsConn(2, jj))
          call computeTangent(x1, x2, tangent)

          ! Store the saved projected point in the full array
          allDist2(i) = dist2
          allProjPoints(:, i) = projPoint
          allTangents(:, i) = tangent
          allElemIDs(i) = jj
          curveMask(i) = 1

       end if

    end do

  end subroutine minDistanceCurve

  !===========================================================
  !===========================================================
  !===========================================================

  subroutine minDistanceCurve_d(nxyz, nCoor, nBars, &
                                xyz, xyzd, coor, coord, barsConn, &
                                allProjPoints, allProjPointsd, &
                                allTangents, allTangentsd, &
                                allElemIDs, curveMask)
    
    ! This is the interface to the forward mode AD.
    ! The user should run the primal function first in order to obtain the correct values for
    ! allElemIDs. This way we can avoid another search for the closest projection.
    !
    ! Ney Secco 2016-11

    use curveUtils_d
    implicit none

    !f2py intent(in) nxyz, nCoor, nBars, xyz, xyzd, coor, coord, barsConn
    !f2py intent(in) allProjPoints, allTangents, allElemIDs, curveMask
    !f2py intent(inout) allProjPointsd, allTangentsd

    !
    !       Subroutine arguments.
    !
    ! Input
    integer(kind=intType), intent(in) :: nxyz, nCoor, nBars
    real(kind=realType), dimension(3,nxyz), intent(in) :: xyz
    real(kind=realType), dimension(3,nxyz), intent(in) :: xyzd
    real(kind=realType), dimension(3,nCoor), intent(in) :: coor
    real(kind=realType), dimension(3,nCoor), intent(in) :: coord
    integer(kind=intType), dimension(2,nBars), intent(in) :: barsConn
    real(kind=realType), dimension(3,nxyz), intent(in) :: allProjPoints
    real(kind=realType), dimension(3,nxyz), intent(in) :: allTangents
    integer(kind=intType), dimension(nxyz), intent(in) :: allElemIDs
    integer(kind=intType), dimension(nxyz), intent(in) :: curveMask

    ! Input/Output
    real(kind=realType), dimension(3,nxyz), intent(inout) :: allProjPointsd
    real(kind=realType), dimension(3,nxyz), intent(inout) :: allTangentsd

    ! Working
    real(kind=realType) :: x(3), xf(3), x1(3), x2(3), tangent(3)
    real(kind=realType) :: xd(3), xfd(3), x1d(3), x2d(3), tangentd(3)
    real(kind=realType) :: u

    integer(kind=intType) :: i, barID, node1, node2

    ! EXECUTION

    ! Loop over all points to be projected
    do i = 1,nxyz

       ! Check if the point was actually projected to this curve
       if (curveMask(i) == 1) then

          ! Recover data from the primal run
          x = xyz(:, i)
          xf = allProjPoints(:, i)
          tangent = allTangents(:, i)
          barID = allElemIDs(i)
          
          ! Get information of the bar element that has the projection
          node1 = barsConn(1, barID)
          node2 = barsConn(2, barID)
          x1 = coor(:, node1)
          x2 = coor(:, node2)
          
          ! Get derivative seeds
          xd = xyzd(:, i)
          x1d = coord(:, node1)
          x2d = coord(:, node2)
          
          ! Compute projection derivatives
          call barProjection_d(x1,x1d,x2,x2d,x,xd,xf,xfd,u)
          
          ! Compute tangent derivatives
          call computeTangent_d(x1, x1d, x2, x2d, tangent, tangentd)

          ! Store derivatives in the full array
          allProjPointsd(:, i) = xfd
          allTangentsd(:, i) = tangentd

       end if

    end do

  end subroutine minDistanceCurve_d

  !===========================================================
  !===========================================================
  !===========================================================

  subroutine minDistanceCurve_b(nxyz, nCoor, nBars, &
                                xyz, xyzb, coor, coorb, barsConn, &
                                allProjPoints, allProjPointsb, &
                                allTangents, allTangentsb, &
                                allElemIDs, curveMask)
    
    ! This is the interface to the reverse mode AD.
    ! The user should run the primal function first in order to obtain the correct values for
    ! allElemIDs. This way we can avoid another search for the closest projection.
    !
    ! ATTENTION:
    ! This code DOES NOT ACCUMULATE DERIVATIVE SEEDS. It just gives the updates. So you
    ! have to add them later on to continue the AD propagation.
    !
    ! Ney Secco 2016-11

    use curveUtils_b
    implicit none

    !f2py intent(in) nxyz, nCoor, nBars, xyz, coor, barsConn
    !f2py intent(in) allProjPoints, allProjPointsb, allTangents, allTangentsb, allElemIDs, curveMask
    !f2py intent(out) xyzb, coorb

    !
    !       Subroutine arguments.
    !
    ! Input
    integer(kind=intType), intent(in) :: nxyz, nCoor, nBars
    real(kind=realType), dimension(3,nxyz), intent(in) :: xyz
    real(kind=realType), dimension(3,nCoor), intent(in) :: coor
    integer(kind=intType), dimension(2,nBars), intent(in) :: barsConn
    real(kind=realType), dimension(3,nxyz), intent(in) :: allProjPoints
    real(kind=realType), dimension(3,nxyz), intent(in) :: allProjPointsb
    real(kind=realType), dimension(3,nxyz), intent(in) :: allTangents
    real(kind=realType), dimension(3,nxyz), intent(in) :: allTangentsb
    integer(kind=intType), dimension(nxyz), intent(in) :: allElemIDs
    integer(kind=intType), dimension(nxyz), intent(in) :: curveMask

    ! Input/Output
    real(kind=realType), dimension(3,nxyz), intent(out) :: xyzb
    real(kind=realType), dimension(3,nCoor), intent(out) :: coorb

    ! Working
    real(kind=realType) :: x(3), xf(3), x1(3), x2(3), tangent(3)
    real(kind=realType) :: xb(3), xfb(3), x1b(3), x2b(3), tangentb(3)
    real(kind=realType) :: x1b_tang(3), x2b_tang(3)
    real(kind=realType) :: u

    integer(kind=intType) :: i, barID, node1, node2

    ! EXECUTION

    ! Initialize derivative seeds
    xyzb = 0.0
    coorb = 0.0

    ! Loop over all points to be projected
    do i = 1,nxyz

       ! Check if the point was actually projected to this curve
       if (curveMask(i) == 1) then

          ! Recover data from the primal run
          x = xyz(:, i)
          xf = allProjPoints(:, i)
          tangent = allTangents(:, i)
          barID = allElemIDs(i)
          
          ! Get information of the bar element that has the projection
          node1 = barsConn(1, barID)
          node2 = barsConn(2, barID)
          x1 = coor(:, node1)
          x2 = coor(:, node2)
          
          ! Get derivative seeds
          xfb = allProjPointsb(:, i)
          tangentb = allTangentsb(:, i)
          
          ! Backpropagate derivatives due to tangent calculation
          x1b_tang = 0.0
          x2b_tang = 0.0
          call computeTangent_b(x1, x1b_tang, x2, x2b_tang, tangent, tangentb)


          ! Call projection function
          x1b = 0.0
          x2b = 0.0
          xb = 0.0
          call barProjection_b(x1,x1b,x2,x2b,x,xb,xf,xfb,u)
          
          ! Store derivatives in the full array
          xyzb(:, i) = xyzb(:, i) + xb
          coorb(:, node1) = coorb(:, node1) + x1b + x1b_tang
          coorb(:, node2) = coorb(:, node2) + x2b + x2b_tang

       end if

    end do

  end subroutine minDistanceCurve_b

    !===============================================================

end module curveSearchAPI
