!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
!
MODULE SOLVEROUTINES_B
  USE PRECISION_B
  IMPLICIT NONE

CONTAINS
!  Differentiation of solve in reverse (adjoint) mode (with options i4 dr8 r8):
!   gradient     of useful results: y
!   with respect to varying inputs: a b
  SUBROUTINE SOLVE_B(a, ab, y, yb, b, bb, n, ipiv)
    USE PRECISION_B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=realtype), INTENT(INOUT) :: a(n, n), y(n), b(n)
    REAL(kind=realtype) :: ab(n, n), yb(n), bb(n)
    INTEGER, INTENT(INOUT) :: ipiv(n)
    INTEGER(kind=inttype) :: ii, jj
    ab = 0.0_8
    bb = 0.0_8
    DO ii=n,1,-1
      DO jj=n,1,-1
        ab(ii, jj) = ab(ii, jj) + b(ii)*yb(ii)
        bb(ii) = bb(ii) + a(ii, jj)*yb(ii)
        yb(ii) = 0.0_8
      END DO
    END DO
  END SUBROUTINE SOLVE_B
  SUBROUTINE SOLVE(a, y, b, n, ipiv)
    USE PRECISION_B
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=realtype), INTENT(INOUT) :: a(n, n), y(n), b(n)
    INTEGER, INTENT(INOUT) :: ipiv(n)
    INTEGER(kind=inttype) :: ii, jj
! THIS CODE LITERALLY DOES NOT MATTER! AS LONG AS 'y' DEPENDS ON 'A'
! and 'b' IT IS FINE.
    DO ii=1,n
      DO jj=1,n
        y(ii) = a(ii, jj)*b(ii)
      END DO
    END DO
  END SUBROUTINE SOLVE
END MODULE SOLVEROUTINES_B