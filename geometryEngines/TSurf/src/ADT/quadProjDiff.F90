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

  subroutine quadProjection_d(x1, x1d, x2, x2d, x3, x3d, x4, x4d, &
                              x, xd, &
                              xf, xfd, u, ud, v, vd, val)

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
    real(kind=realType), intent(in) :: u, v, val
    real(kind=realType), dimension(3), intent(in) :: xf

    ! Output variables
    real(kind=realType), dimension(3), intent(out) :: xfd
    real(kind=realType), intent(out) :: ud, vd

    ! Working variables
    real(kind=realType), dimension(3) :: x21, x41, x3142
    real(kind=realType), dimension(3) :: x21d, x41d, x3142d, xftemp
    real(kind=realType), dimension(2) :: error, errord, aux2
    real(kind=realType), dimension(2,2) :: invJac
    real(kind=realType) :: du, dv

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
    xftemp = xf

    call quadProjSubIter_d(x1, x1d, x2, x2d, x3, x3d, x4, x4d, &
                           x, xd, u, ud, v, vd, &
                           du, dv, error, errord, invJac, xftemp, xfd)

    print *,'x1'
    print *,x1
    print *,'x2'
    print *,x2
    print *,'x3'
    print *,x3
    print *,'x4'
    print *,x4
    print *,'x'
    print *,x
    print *,'x1d'
    print *,x1d
    print *,'x2d'
    print *,x2d
    print *,'x3d'
    print *,x3d
    print *,'x4d'
    print *,x4d
    print *,'xd'
    print *,xd
    print *,'invJac'
    print *,invJac(1,:)
    print *,invJac(2,:)
    print *,'error'
    print *,error
    print *,'errord'
    print *,errord
    print *,'u'
    print *,u
    print *,'v'
    print *,v
    print *,'du'
    print *,du
    print *,'dv'
    print *,dv

    ! Now we can compute the values of the seed ud and vd by
    ! solving the "adjoint" system. The inverse Jacobian that
    ! we compute inside the original function already represents
    ! the inversion of the right hand side of the system.
    aux2 = matmul(invJac, errord)
    ud = -aux2(1)
    vd = -aux2(2)

    ! Determine auxiliary vectors and its derivatives
    x21 = x2 - x1
    x21d = x2d - x1d
    x41 = x4 - x1
    x41d = x4d - x1d
    x3142 = x3 - x1 - x21 - x41
    x3142d = x3d - x1d - x21d - x41d

    ! Derivatives with respect to xf
    ! Remember that xf = x1 + u*x21 + v*x41 + u*v*x3142
    xfd = x1d + &
          ud*x21 + u*x21d + &
          vd*x41 + v*x41d + &
          ud*v*x3142 + u*vd*x3142 + u*v*x3142d

  end subroutine quadProjection_d

  !=========================================================
  !=========================================================
  !=========================================================

  subroutine quadProjection_b(x1, x1b, x2, x2b, x3, x3b, x4, x4b, &
                              x, xb, &
                              xf, xfb, u, ub, v, vb, val)

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
    real(kind=realType), intent(in) :: u, v, val
    real(kind=realType), intent(in) :: ub, vb
    real(kind=realType), dimension(3), intent(in) :: xf
    real(kind=realType), dimension(3), intent(in) :: xfb

    ! Output variables
    real(kind=realType), dimension(3), intent(out) :: x1b, x2b, x3b, x4b
    real(kind=realType), dimension(3), intent(out) :: xb

    ! Working variables
    real(kind=realType), dimension(3) :: xf_temp
    real(kind=realType), dimension(3) :: xfb_temp
    real(kind=realType) :: ub_temp, vb_temp
    real(kind=realType), dimension(2) :: error, errorb, aux2
    real(kind=realType), dimension(2,2) :: invJac
    real(kind=realType) :: du, dv

    ! EXECUTION

    ! Here we need to solve an adjoint for the projection function
    ! to avoid using the Newton iteration all over again.
    !
    ! The Newton iteration on the original function are used to solve
    ! grad(dist2) = 0, which is a non-linear system with two equations
    ! and two unknowns (u and v). We can create adjoint systems for
    ! this problem to propagate derivatives.

    ! Set temporary variable to xf to avoid overwriting
    xf_temp = xf
    xfb_temp = xfb
    ub_temp = ub
    vb_temp = vb

    ! We want the error to remain at 0, even when changing inputs/outputs.
    ! Therefore, we can set its seed to zero
    errorb = 0.0

    print *,'RUNNING REVERSE MODE'
    print *,'xfb'
    print *,xfb
    print *,'ub'
    print *,ub
    print *,'vb'
    print *,vb
    print *,'x1'
    print *,x1
    print *,'x2'
    print *,x2
    print *,'x3'
    print *,x3
    print *,'x4'
    print *,x4
    print *,'x'
    print *,x
    print *,'u'
    print *,u
    print *,'v'
    print *,v

    ! Check new error derivative
    call quadProjSubIter_b(x1, x1b, x2, x2b, x3, x3b, x4, x4b, &
                           x, xb, u, ub_temp, v, vb_temp, &
                           du, dv, error, errorb, invJac, xf_temp, xfb_temp)

    print *,'xfb'
    print *,xfb
    print *,'ub'
    print *,ub
    print *,'vb'
    print *,vb
    print *,'x1'
    print *,x1
    print *,'x2'
    print *,x2
    print *,'x3'
    print *,x3
    print *,'x4'
    print *,x4
    print *,'x'
    print *,x
    print *,'x1b'
    print *,x1b
    print *,'x2b'
    print *,x2b
    print *,'x3b'
    print *,x3b
    print *,'x4b'
    print *,x4b
    print *,'xb'
    print *,xb
    print *,'invJac'
    print *,invJac(1,:)
    print *,invJac(2,:)
    print *,'error'
    print *,error
    print *,'errorb'
    print *,errorb
    print *,'du'
    print *,du
    print *,'dv'
    print *,dv

  end subroutine quadProjection_b



end module quadProjDiff
