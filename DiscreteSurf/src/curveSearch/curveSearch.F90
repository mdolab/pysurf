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

        subroutine minDistanceCurve(nCoor, nBars, coor, barNodes, allProjPoints, allDist)
!
!       ****************************************************************
!       *                                                              *
!       *                                                              *
!       ****************************************************************
!
        implicit none

        !f2py intent(in) nCoor, nBars, coor, barNodes
        !f2py intent(out) projPoints

!
!       Subroutine arguments.
!
        ! Input
        integer(kind=intType), intent(in) :: nCoor, nBars
        real(kind=realType), dimension(3,nCoor), intent(in) :: coor
        real(kind=realType), dimension(3,nBars+1), intent(in) :: barNodes

        ! Output
        real(kind=realType), dimension(3,nCoor), intent(inout) :: allProjPoints
        real(kind=realType), dimension(nCoor), intent(inout) :: allDist

        ! Working
        real(kind=realType), dimension(nCoor) :: connect_count
        real(kind=realType) :: pt(3), dist, ppt(3), x21(3), vec(3)
        real(kind=realType) :: x1(3), x2(3), p1(3), u, mag2, d(3)
        real(kind=realType) :: newProjPoint(3), newDist, projPoint(3)
        integer(kind=intType) :: i, j

        !===============================================================

        ! Loop over all points to be projected
        do i = 1,nCoor
          pt = coor(:, i)
          dist = allDist(i)

          ! Loop over all bars
          do j = 1,nBars
            x1 = barNodes(:, j)
            x2 = barNodes(:, j+1)

            x21 = x2 - x1
            mag2 = dot_product(x21, x21)
            p1 =  pt - x1

            u = dot_product(x21, p1) / mag2

            if (u .lt. 0) then
              u = 0
            else if (u .gt. 1) then
              u = 1
            end if

            newProjPoint = x1 + u*x21
            vec = newProjPoint - pt
            newDist = sqrt(dot_product(vec, vec))

            if (newDist .lt. dist) then
              dist = newDist
              projPoint = newProjPoint
            end if

          end do

          allDist(i) = dist
          allProjPoints(:, i) = projPoint

        end do

        end subroutine minDistanceCurve


      end module curveProj
