module curveUtils

    use precision
    use constants
    implicit none

contains

    subroutine barProjection(x1, x2, x, xf, u)

        ! This subroutine projects the point x onto the bar element defined by points
        ! x1 and x2, to obtain the projected point xf.
        !
        ! Ney Secco 2016-11

        use utilities ! This will bring dot
        implicit none

        ! DECLARATIONS

        ! Input variables
        real(kind=realType), dimension(3), intent(in) :: x1, x2, x

        ! Output variables
        real(kind=realType), dimension(3), intent(out) :: xf
        real(kind=realType), intent(out) :: u

        ! Working variables
        real(kind=realType), dimension(3) :: x21, p1, vec, dummyVec
        real(kind=realType) :: mag2, dotResult

        ! EXECUTION

        ! Get the relative vectors for the bar element and the point
        x21 = x2 - x1
        dummyVec = x21
        call dot(x21, dummyVec, mag2)
        p1 = x - x1

        ! Compute the amount of the point that projects onto the element
        call dot(x21, p1, dotResult)
        u = dotResult / mag2

        ! Set the projected point to either the start or end node if
        ! the projection lies outside the node
        if (u .lt. 0) then
            u = 0
        else if (u .gt. 1) then
            u = 1
        end if

        ! Compute the new projected point coordinates in the global frame
        xf = x1 + u * x21

    end subroutine barProjection

    subroutine computeTangent(x1, x2, tangent)

        ! This subroutine computes a normalized vector pointing from x1 to x2
        !
        ! Ney Secco 2016-11

        use utilities
        implicit none

        ! DECLARATIONS

        ! Input variables
        real(kind=realType), dimension(3), intent(in) :: x1, x2

        ! Output variables
        real(kind=realType), dimension(3), intent(out) :: tangent

        ! Working variables
        real(kind=realType), dimension(3) :: x21, dummyVec
        real(kind=realType) :: dotResult

        ! EXECUTION

        ! Get the relative vectors for the bar element
        x21 = x2 - x1

        ! Normalize vector (dot defined in utilities.F90)
        dummyVec = x21
        call dot(x21, dummyVec, dotResult)
        tangent = x21 / sqrt(dotResult)

    end subroutine computeTangent

end module curveUtils
