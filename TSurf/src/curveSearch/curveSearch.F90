!
!     ******************************************************************
!     *                                                                *
!     * File:          curveSearch.F90                                 *
!     * Author:        John Jasa                                       *
!     * Starting date: 07-27-2016                                      *
!     * Last modified: 07-27-2016                                      *
!     *                                                                *
!     ******************************************************************
!
      module curveProj

      use precision
      implicit none

      !=================================================================

      contains

        subroutine minDistanceCurve(nxyz, nCoor, nBars, xyz, coor, barsConn, allProjPoints, tangents, allDist)
!
!       ****************************************************************
!       *                                                              *
!       *                                                              *
!       ****************************************************************
!
        implicit none

        !f2py intent(in) nxyz, nCoor, nBars, xyz, coor, barsConn
        !f2py intent(out) projPoints, tangents

!
!       Subroutine arguments.
!
        ! Input
        integer(kind=intType), intent(in) :: nxyz, nCoor, nBars
        real(kind=realType), dimension(3,nxyz), intent(in) :: xyz
        real(kind=realType), dimension(3,nCoor), intent(in) :: coor
        integer(kind=intType), dimension(2,nBars), intent(in) :: barsConn

        ! Output
        real(kind=realType), dimension(3,nxyz), intent(inout) :: allProjPoints
        real(kind=realType), dimension(3,nxyz), intent(inout) :: tangents
        real(kind=realType), dimension(nxyz), intent(inout) :: allDist

        ! Working
        real(kind=realType) :: pt(3), dist, x21(3), vec(3)
        real(kind=realType) :: x1(3), x2(3), p1(3), u, mag2
        real(kind=realType) :: newProjPoint(3), newDist, projPoint(3)
        real(kind=realType) :: nodalTangents(3, nBars)
        integer(kind=intType) :: i, j, jj

        !===============================================================

        do j = 1,nBars

          ! Get the absolute nodal points from the coordinate array
          x1 = coor(:, barsConn(1, j))
          x2 = coor(:, barsConn(2, j))

          ! Get the relative vectors for the bar element normalize
          x21 = x2 - x1
          nodalTangents(:, j) = x21 / sqrt(dot_product(x21, x21))
        end do

        ! Loop over all points to be projected
        do i = 1,nxyz

          ! Obtain the raw point to be projected and its current minimum
          ! distance from a curve
          pt = xyz(:, i)
          dist = allDist(i)

          ! Loop over all individual bar elements that define the curve
          do j = 1,nBars

            ! Get the absolute nodal points from the coordinate array
            x1 = coor(:, barsConn(1, j))
            x2 = coor(:, barsConn(2, j))

            ! Get the relative vectors for the bar element and the point
            x21 = x2 - x1
            mag2 = dot_product(x21, x21)
            p1 =  pt - x1

            ! Compute the amount of the point that projects onto the element
            u = dot_product(x21, p1) / mag2

            ! Set the projected point to either the start or end node if
            ! the projection lies outside the node
            if (u .lt. 0) then
              u = 0
            else if (u .gt. 1) then
              u = 1
            end if

            ! Compute the new projected point coordinates in the global frame
            newProjPoint = x1 + u*x21
            vec = newProjPoint - pt
            newDist = dot_product(vec, vec)

            ! If the distance for this new projected point is less than the
            ! provided distance, save the new point
            if (newDist .lt. dist) then
              dist = newDist
              projPoint = newProjPoint
              jj = j
            end if

          end do

          ! Store the saved projected point in the full array
          allDist(i) = dist
          allProjPoints(:, i) = projPoint
          tangents(:, i) = nodalTangents(:, jj)

        end do

        end subroutine minDistanceCurve

      end module curveProj
