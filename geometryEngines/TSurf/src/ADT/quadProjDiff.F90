! This module only contains hand-differentiated quad projection subroutines
! I had to place it here so that it could use differentiated versions of
! subroutines present in adtProjections.F90, such as quadProjSubIter.
!
! Ney Secco 2016-11

module quadProjDiff

use precision
use constants
implicit none

contains

  subroutine quadProjection_d(x1, x1d, x2, x2d, x3, x3d, x4, x4d, x, xd, &
                              xfd, u, ud, v, vd)

    ! This subroutine computes the projection (xf) of a given point (x) into the
    ! quad defined by four nodes (x1, x2, x3, and x4).
    ! This is differentiated by hand.
    ! 
    ! INPUTS:
    !
    ! x1: real(3) -> Coordinates (X,Y,Z) of the first quad node.
    !
    ! x2: real(3) -> Coordinates (X,Y,Z) of the second quad node.
    !
    ! x3: real(3) -> Coordinates (X,Y,Z) of the third quad node.
    !
    ! x4: real(3) -> Coordinates (X,Y,Z) of the fourth quad node.
    !
    ! x: real(3) -> Coordinates (X,Y,Z) of the point that should be projected.
    !
    ! OUTPUTS:
    !
    ! xf: real(3) -> Coordinates (X,Y,Z) of the projected point.
    !
    ! u: real -> Parametric coordinate of the projected point on the quad element.
    !
    ! v: real -> Parametric coordinate of the projected point on the quad element.
    !
    ! val: real -> Distance**2 between the point (x) and its projection (xf).
    !
    ! Ney Secco - 2016-10

    use adtProjections_d
    implicit none

    ! DECLARATIONS

    ! Input variables
    real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
    real(kind=realType), dimension(3), intent(in) :: x1d, x2d, x3d, x4d
    real(kind=realType), dimension(3), intent(in) :: x
    real(kind=realType), dimension(3), intent(in) :: xd
    real(kind=realType), intent(in) :: u, v

    ! Output variables
    real(kind=realType), dimension(3), intent(out) :: xfd
    real(kind=realType), intent(out) :: ud, vd

    ! Working variables
    real(kind=realType), dimension(2) :: residual, residuald, aux2
    real(kind=realType), dimension(2,2) :: invJac
    real(kind=realType), dimension(3) :: xf

    ! EXECUTION

    ! Here we need to solve an adjoint for the projection function
    ! to avoid using the Newton iteration all over again.
    !
    ! The Newton iteration on the original function are used to solve
    ! grad(dist2) = 0, which is a non-linear system with two equations
    ! and two unknowns (u and v). We can create adjoint systems for
    ! this problem to propagate derivatives.

    ! Call projection function at the final u,v to get inverse Jacobian
    ! The inverse Jacobian will be used to solve the adjoint system and
    ! get the correct seed of ud and vd.
    ! We use the same call to get the partial derivative of the constraint
    ! function with respect to the input variables only (errord). Thus we can set the
    ! seeds of the state variables to zero (ud = vd = 0)
    ud = 0.0
    vd = 0.0
    call quadProjResidual_d(x1, x1d, x2, x2d, x3, x3d, x4, x4d, x, xd, &
                            u, ud, v, vd, &
                            residual, residuald, invJac)

    ! Now we can compute the values of the seed ud and vd by
    ! solving the "adjoint" system. The inverse Jacobian that
    ! we compute inside the original function already represents
    ! the inversion of the right hand side of the system.
    aux2 = matmul(invJac, residuald)
    ud = -aux2(1)
    vd = -aux2(2)

    ! Compute derivatives of the outputs
    call quadProjOutput_d(x1, x1d, x2, x2d, x3, x3d, x4, x4d, &
                          u, ud, v, vd, &
                          xf, xfd)

  end subroutine quadProjection_d

  !=========================================================
  !=========================================================
  !=========================================================

  subroutine quadProjection_b(x1, x1b, x2, x2b, x3, x3b, x4, x4b, x, xb, &
                              xf, xfb, u, ub, v, vb)

    ! This subroutine computes the projection (xf) of a given point (x) into the
    ! quad defined by four nodes (x1, x2, x3, and x4).
    ! This is differentiated by hand.
    ! 
    ! INPUTS:
    !
    ! x1: real(3) -> Coordinates (X,Y,Z) of the first quad node.
    !
    ! x2: real(3) -> Coordinates (X,Y,Z) of the second quad node.
    !
    ! x3: real(3) -> Coordinates (X,Y,Z) of the third quad node.
    !
    ! x4: real(3) -> Coordinates (X,Y,Z) of the fourth quad node.
    !
    ! x: real(3) -> Coordinates (X,Y,Z) of the point that should be projected.
    !
    ! OUTPUTS:
    !
    ! xf: real(3) -> Coordinates (X,Y,Z) of the projected point.
    !
    ! u: real -> Parametric coordinate of the projected point on the quad element.
    !
    ! v: real -> Parametric coordinate of the projected point on the quad element.
    !
    ! val: real -> Distance**2 between the point (x) and its projection (xf).
    !
    ! Ney Secco - 2016-10

    use adtProjections_b
    implicit none

    ! DECLARATIONS

    ! Input variables
    real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
    real(kind=realType), dimension(3), intent(in) :: x
    real(kind=realType), intent(in) :: u, v
    real(kind=realType), intent(in) :: ub, vb
    real(kind=realType), dimension(3), intent(in) :: xf
    real(kind=realType), dimension(3), intent(in) :: xfb

    ! Output variables
    real(kind=realType), dimension(3), intent(out) :: x1b, x2b, x3b, x4b
    real(kind=realType), dimension(3), intent(out) :: xb

    ! Working variables
    real(kind=realType), dimension(3) :: xf_temp
    real(kind=realType), dimension(3) :: xfb_temp
    real(kind=realType), dimension(3) :: x1b_partial, x2b_partial, x3b_partial, x4b_partial

    real(kind=realType) :: ub_temp, vb_temp, ub_partial, vb_partial
    real(kind=realType), dimension(2) :: residual, residualb
    real(kind=realType), dimension(2,2) :: Jac, invJac
    real(kind=realType), dimension(2,15) :: dgdx
    real(kind=realType), dimension(15) :: xb2full
    real(kind=realType) :: dotResult

    ! EXECUTION

    ! Here we need to solve an adjoint for the projection function
    ! to avoid using the Newton iteration all over again.
    !
    ! The Newton iteration on the original function are used to solve
    ! grad(dist2) = 0, which is a non-linear system with two equations
    ! and two unknowns (u and v). We can create adjoint systems for
    ! this problem to propagate derivatives.

    ! We want the error to remain at 0, even when changing inputs/outputs.
    ! Therefore, we can set its seed to zero

    ! Get partial derivatives of the outputs with respect to input and state variables.
    ! We need a temporary copy of the output derivatives since the code modifies it.
    xf_temp = xf
    xfb_temp = xfb

    call quadProjOutput_b(x1, x1b_partial, x2, x2b_partial, x3, x3b_partial, x4, x4b_partial, &
                          u, ub_partial, v, vb_partial, &
                          xf_temp, xfb_temp)
    ub_partial = ub_partial + ub
    vb_partial = vb_partial + vb

    ! Now we need to assemble the Jacobian (dg/du) and the RHS of the adjoint equation (dg/dx)
    residualb = [1.0, 0.0]
    call quadProjResidual_b(x1, x1b, x2, x2b, x3, x3b, x4, x4b, &
                            x, xb, u, ub_temp, v, vb_temp, &
                            residual, residualb, invJac)
    dgdx(1,1:3) = x1b
    dgdx(1,4:6) = x2b
    dgdx(1,7:9) = x3b
    dgdx(1,10:12) = x4b
    dgdx(1,13:15) = xb
    Jac(1,:) = [ub_temp, vb_temp]

    residualb = [0.0, 1.0]
    call quadProjResidual_b(x1, x1b, x2, x2b, x3, x3b, x4, x4b, &
                            x, xb, u, ub_temp, v, vb_temp, &
                            residual, residualb, invJac)
    dgdx(2,1:3) = x1b
    dgdx(2,4:6) = x2b
    dgdx(2,7:9) = x3b
    dgdx(2,10:12) = x4b
    dgdx(2,13:15) = xb
    Jac(2,:) = [ub_temp, vb_temp]

    ! Invert the Jacobian
    call invert2x2(Jac, invJac)

    ! Use the adjoint to compute total derivatives
    xb2full = matmul([ub_partial, vb_partial], matmul(invJac,dgdx))

    x1b = x1b_partial - xb2full(1:3)
    x2b = x2b_partial - xb2full(4:6)
    x3b = x3b_partial - xb2full(7:9)
    x4b = x4b_partial - xb2full(10:12)
    xb = -xb2full(13:15)

  end subroutine quadProjection_b



end module quadProjDiff
