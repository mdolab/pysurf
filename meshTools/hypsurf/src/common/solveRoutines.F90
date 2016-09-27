module solveRoutines


use precision
implicit none
contains

#ifndef USE_TAPENADE

  ! TODO : Refactor the solve routines to take in "m" such that we can
  ! have more RHS.
  subroutine solve(A, y, b, n, ipiv)

    implicit none
    real(kind=realType), intent(inout) :: A(n,n), y(n), b(n)
    integer, intent(in) :: n
    integer, intent(inout) :: ipiv(n)
    integer :: info, nrhs
    character :: trans
    real(kind=realType) :: Acopy(n,n)
    Acopy = A
    nrhs = 1
    call DGETRF(N, N, Acopy, N, ipiv, info)
    y = b
    call dgetrs('N', N, nrhs, ACopy, N, ipiv, y, N, info)

  end subroutine solve


  subroutine solve_b(A, Ab, y, yb, b, bb, n, ipiv)
    implicit none
    real(kind=realType):: A(n,n), Ab(n, n), y(n), yb(n)
    real(kind=realType) :: b(n), bb(n)
    integer, intent(inout) :: ipiv(n)
    integer, intent(in) :: n
    integer :: info, i, j
    real(kind=realType) :: incrbb(n), ACopy(n, n), RHS(n), zero

    zero = 0.

    Acopy = A
    call DGETRF(N, N, ACopy, N, ipiv, info)

    ! Transpose solve for the bb increment.
    incrbb = yb

    call dgetrs('T', N, 1, ACopy, N, ipiv, incrbb, N, info)

    bb = bb + incrbb

    RHS = b
    call dgetrs('N', N, 1, ACopy, N, ipiv, RHS, N, info)

    ! RHS now contains y
    y = RHS

    ! Accumulate into Ab
    do j=1,n
       do i=1,n
          Ab(i,j) = Ab(i,j) - RHS(j)*incrbb(i)
       end do
    end do
    yb = zero
  end subroutine solve_b


  subroutine solve_d(A, Ad, y, yd, b, bd, n, ipiv)
    implicit none
    real(kind=realType), intent(inout) :: A(n,n), Ad(n,n), y(n), yd(n)
    real(kind=realType) :: b(n), bd(n)
    integer, intent(inout) :: ipiv(n)
    integer, intent(in) :: n
    ! Working
    integer :: info, i,j
    real(kind=realType) :: RHSd(n), ACopy(n, n)

    ! The forward AD for Ax = b is
    !    Ad*y + A*yd = bd
    ! => yd = A^{-1}*(bd - Ad*y)

    ! Calculate y since we need to use that on the RHS, call regular solve
    call solve(A, y, b, n, ipiv)

    ! Now find (bd - Ad*y)
    do i=1,n
      RHSd(i) = bd(i)
      do j=1,n
        RHSd(i) = RHSd(i) - Ad(i,j)*y(j)
      end do
    end do

    ! Now solve yd = A^{-1}*(bd - Ad*y) with regular solve
    call solve(A, yd, RHSd, n, ipiv)

  end subroutine solve_d



  ! ----------------------------------------------------------------------
  ! These are the dummy routines that tapenade sees just to get
  ! the dependencies correct
  ! ----------------------------------------------------------------------
#else
  subroutine solve(A, y, b, n, ipiv)
    use precision
    implicit none
    real(kind=realType), intent(inout) :: A(n,n), y(n), b(n)
    integer, intent(in) :: n, ipiv(n)
    integer(kind=intType) :: info

    ! THIS CODE LITERALLY DOES NOT MATTER! AS LONG AS 'y' DEPENDS ON 'A'
    ! and 'b' IT IS FINE.
    y(1) = A(1,1)*b(1)

  end subroutine solve


#endif

end module
