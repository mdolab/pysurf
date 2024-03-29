!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 22 Aug 2023 15:51
!
MODULE CURVEUTILS_D
  USE PRECISION
  USE CONSTANTS
  IMPLICIT NONE

CONTAINS
!  Differentiation of barprojection in forward (tangent) mode:
!   variations   of useful results: xf
!   with respect to varying inputs: x x1 x2
!   RW status of diff variables: x:in xf:out x1:in x2:in
  SUBROUTINE BARPROJECTION_D(x1, x1d, x2, x2d, x, xd, xf, xfd, u)
! This subroutine projects the point x onto the bar element defined by points
! x1 and x2, to obtain the projected point xf.
!
! Ney Secco 2016-11
! This will bring dot
    USE UTILITIES_D
    IMPLICIT NONE
! DECLARATIONS
! Input variables
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1, x2, x
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1d, x2d, xd
! Output variables
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: xf
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: xfd
    REAL(kind=realtype), INTENT(OUT) :: u
    REAL(kind=realtype) :: ud
! Working variables
    REAL(kind=realtype), DIMENSION(3) :: x21, p1, vec, dummyvec
    REAL(kind=realtype), DIMENSION(3) :: x21d, p1d, dummyvecd
    REAL(kind=realtype) :: mag2, dotresult
    REAL(kind=realtype) :: mag2d, dotresultd
! EXECUTION
! Get the relative vectors for the bar element and the point
    x21d = x2d - x1d
    x21 = x2 - x1
    dummyvecd = x21d
    dummyvec = x21
    CALL DOT_D0(x21, x21d, dummyvec, dummyvecd, mag2, mag2d)
    p1d = xd - x1d
    p1 = x - x1
! Compute the amount of the point that projects onto the element
    CALL DOT_D0(x21, x21d, p1, p1d, dotresult, dotresultd)
    ud = (dotresultd-dotresult*mag2d/mag2)/mag2
    u = dotresult/mag2
! Set the projected point to either the start or end node if
! the projection lies outside the node
    IF (u .LT. 0) THEN
      u = 0
      ud = 0.0_8
    ELSE IF (u .GT. 1) THEN
      u = 1
      ud = 0.0_8
    END IF
! Compute the new projected point coordinates in the global frame
    xfd = x1d + x21*ud + u*x21d
    xf = x1 + u*x21
  END SUBROUTINE BARPROJECTION_D

  SUBROUTINE BARPROJECTION(x1, x2, x, xf, u)
! This subroutine projects the point x onto the bar element defined by points
! x1 and x2, to obtain the projected point xf.
!
! Ney Secco 2016-11
! This will bring dot
    USE UTILITIES_D
    IMPLICIT NONE
! DECLARATIONS
! Input variables
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1, x2, x
! Output variables
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: xf
    REAL(kind=realtype), INTENT(OUT) :: u
! Working variables
    REAL(kind=realtype), DIMENSION(3) :: x21, p1, vec, dummyvec
    REAL(kind=realtype) :: mag2, dotresult
! EXECUTION
! Get the relative vectors for the bar element and the point
    x21 = x2 - x1
    dummyvec = x21
    CALL DOT(x21, dummyvec, mag2)
    p1 = x - x1
! Compute the amount of the point that projects onto the element
    CALL DOT(x21, p1, dotresult)
    u = dotresult/mag2
! Set the projected point to either the start or end node if
! the projection lies outside the node
    IF (u .LT. 0) THEN
      u = 0
    ELSE IF (u .GT. 1) THEN
      u = 1
    END IF
! Compute the new projected point coordinates in the global frame
    xf = x1 + u*x21
  END SUBROUTINE BARPROJECTION

!  Differentiation of computetangent in forward (tangent) mode:
!   variations   of useful results: tangent
!   with respect to varying inputs: x1 x2
!   RW status of diff variables: tangent:out x1:in x2:in
  SUBROUTINE COMPUTETANGENT_D(x1, x1d, x2, x2d, tangent, tangentd)
! This subroutine computes a normalized vector pointing from x1 to x2
!
! Ney Secco 2016-11
    USE UTILITIES_D
    IMPLICIT NONE
! DECLARATIONS
! Input variables
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1, x2
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1d, x2d
! Output variables
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: tangent
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: tangentd
! Working variables
    REAL(kind=realtype), DIMENSION(3) :: x21, dummyvec
    REAL(kind=realtype), DIMENSION(3) :: x21d, dummyvecd
    REAL(kind=realtype) :: dotresult
    REAL(kind=realtype) :: dotresultd
    INTRINSIC SQRT
    REAL(kind=realtype) :: result1
    REAL(kind=realtype) :: result1d
    REAL(kind=realtype) :: temp
! EXECUTION
! Get the relative vectors for the bar element
    x21d = x2d - x1d
    x21 = x2 - x1
! Normalize vector (dot defined in utilities.F90)
    dummyvecd = x21d
    dummyvec = x21
    CALL DOT_D0(x21, x21d, dummyvec, dummyvecd, dotresult, dotresultd)
    temp = SQRT(dotresult)
    IF (dotresult .EQ. 0.0) THEN
      result1d = 0.0_8
    ELSE
      result1d = dotresultd/(2.0*temp)
    END IF
    result1 = temp
    tangentd = (x21d-x21*result1d/result1)/result1
    tangent = x21/result1
  END SUBROUTINE COMPUTETANGENT_D

  SUBROUTINE COMPUTETANGENT(x1, x2, tangent)
! This subroutine computes a normalized vector pointing from x1 to x2
!
! Ney Secco 2016-11
    USE UTILITIES_D
    IMPLICIT NONE
! DECLARATIONS
! Input variables
    REAL(kind=realtype), DIMENSION(3), INTENT(IN) :: x1, x2
! Output variables
    REAL(kind=realtype), DIMENSION(3), INTENT(OUT) :: tangent
! Working variables
    REAL(kind=realtype), DIMENSION(3) :: x21, dummyvec
    REAL(kind=realtype) :: dotresult
    INTRINSIC SQRT
    REAL(kind=realtype) :: result1
! EXECUTION
! Get the relative vectors for the bar element
    x21 = x2 - x1
! Normalize vector (dot defined in utilities.F90)
    dummyvec = x21
    CALL DOT(x21, dummyvec, dotresult)
    result1 = SQRT(dotresult)
    tangent = x21/result1
  END SUBROUTINE COMPUTETANGENT

END MODULE CURVEUTILS_D

