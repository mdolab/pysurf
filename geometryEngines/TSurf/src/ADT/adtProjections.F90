!
!     ******************************************************************
!     *                                                                *
!     * File:          adtProjections.F90                              *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 10-27-2016                                      *
!     * Last modified: 10-27-2016                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtProjections
!
!     ******************************************************************
!     *                                                                *
!     * Module, which contains small subroutines which perform         *
!     * projection tasks. They are grouped here so they could be       *
!     * easily differentiated.                                         *
!     *                                                                *
!     ******************************************************************
!
      use precision
      use constants
      implicit none

      contains

        !===============================================================


        subroutine triaProjection(x1, x2, x3, &
                                  x, &
                                  xf, u, v, val)

          ! This subroutine computes the projection (xf) of a given point (x) into the
          ! triangle defined by three nodes (x1, x2, and x3).
          ! 
          ! INPUTS:
          !
          ! x1: real(3) -> Coordinates (X,Y,Z) of the first triangle node.
          !
          ! x2: real(3) -> Coordinates (X,Y,Z) of the second triangle node.
          !
          ! x3: real(3) -> Coordinates (X,Y,Z) of the third triangle node.
          !
          ! x: real(3) -> Coordinates (X,Y,Z) of the point that should be projected.
          !
          ! OUTPUTS:
          !
          ! xf: real(3) -> Coordinates (X,Y,Z) of the projected point.
          !
          ! u: real -> Parametric coordinate of the projected point on the triangle element.
          !
          ! v: real -> Parametric coordinate of the projected point on the triangle element.
          !
          ! val: real -> Distance**2 between the point (x) and its projection (xf).
          !
          ! Ney Secco - 2016-10
          
          implicit none

          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(3), intent(in) :: x1, x2, x3
          real(kind=realType), dimension(3), intent(in) :: x

          ! Output variables
          real(kind=realType), dimension(3), intent(out) :: xf
          real(kind=realType), intent(out) :: u, v, val

          ! Working variables
          real(kind=realType), dimension(3) :: a, b, vf, vt
          real(kind=realType), dimension(3) :: an, bn, norm
          real(kind=realType) :: invLen, vn, dx, dy, dz, uv

          ! EXECUTION

          ! Determine the tangent vectors in u- and v-direction.
          ! Store these in a and b respectively.
          a = x2 - x1
          b = x3 - x1
 
          ! Determine the normal vector of the face by taking the
          ! cross product of a and b. Afterwards this vector will
          ! be scaled to a unit vector.

          norm(1) = a(2)*b(3) - a(3)*b(2)
          norm(2) = a(3)*b(1) - a(1)*b(3)
          norm(3) = a(1)*b(2) - a(2)*b(1)

          invLen = one/max(eps,sqrt(norm(1)*norm(1) &
               +                        norm(2)*norm(2) &
               +                        norm(3)*norm(3)))

          norm(1) = norm(1)*invLen
          norm(2) = norm(2)*invLen
          norm(3) = norm(3)*invLen

          ! Determine the vector vf from xf to given coordinate.

          vf(1) = x(1) - x1(1)
          vf(2) = x(2) - x1(2)
          vf(3) = x(3) - x1(3)

          ! Determine the projection of the vector vf onto
          ! the face.

          vn = vf(1)*norm(1) + vf(2)*norm(2) + vf(3)*norm(3)
          vt(1) = vf(1) - vn*norm(1)
          vt(2) = vf(2) - vn*norm(2)
          vt(3) = vf(3) - vn*norm(3)

          ! The vector vt points from the current point on the
          ! face to the new point. However this new point lies on
          ! the plane determined by the vectors a and b, but not
          ! necessarily on the face itself. The new point on the
          ! face is obtained by projecting the point in the a-b
          ! plane onto the face. this can be done by determining
          ! the coefficients du and dv, such that vt = du*a + dv*b.
          ! To solve du and dv the vectors normal to a and b
          ! inside the plane ab are needed.

          an(1) = a(2)*norm(3) - a(3)*norm(2)
          an(2) = a(3)*norm(1) - a(1)*norm(3)
          an(3) = a(1)*norm(2) - a(2)*norm(1)

          bn(1) = b(2)*norm(3) - b(3)*norm(2)
          bn(2) = b(3)*norm(1) - b(1)*norm(3)
          bn(3) = b(1)*norm(2) - b(2)*norm(1)

          ! Solve parametric coordinates u and v.
          ! The clipping of vn should not be
          ! active, as this would mean that the vectors a and b
          ! are parallel. This corresponds to a tria degenerated
          ! to a line, which should not occur in the surface mesh.

          vn = a(1)*bn(1) + a(2)*bn(2) + a(3)*bn(3)
          vn = sign(max(eps,abs(vn)),vn)
          u = (vt(1)*bn(1) + vt(2)*bn(2) + vt(3)*bn(3))/vn

          vn = b(1)*an(1) + b(2)*an(2) + b(3)*an(3)
          vn = sign(max(eps,abs(vn)),vn)
          v = (vt(1)*an(1) + vt(2)*an(2) + vt(3)*an(3))/vn

          ! Triangles should be bounded by the line u + v = 1

          uv = u + v
          if (uv > one) then
             u = u / uv
             v = v / uv
          end if

          ! Determine the new parameter values u and v. These
          ! are limited to 0 <= (u,v) <= 1.

          u = min(one,max(zero,u))
          v = min(one,max(zero,v))

          ! Determine the new coordinates of the point xf.
          xf = x1 + u*a + v*b
          
          ! Compute the distance squared between the given
          ! coordinate and the point xf.

          dx = x(1) - xf(1)
          dy = x(2) - xf(2)
          dz = x(3) - xf(3)

          val = dx*dx + dy*dy + dz*dz

        end subroutine triaProjection

        !***************************************************************
        !***************************************************************

        subroutine quadProjection(x1, x2, x3, x4, &
                                  x, &
                                  xf, u, v, val)

          ! This subroutine computes the projection (xf) of a given point (x) into the
          ! quad defined by four nodes (x1, x2, x3, and x4).
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
          
          implicit none

          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
          real(kind=realType), dimension(3), intent(in) :: x

          ! Output variables
          real(kind=realType), dimension(3), intent(out) :: xf
          real(kind=realType), intent(out) :: u, v, val

          ! Working variables
          real(kind=realType), dimension(2) :: residual
          real(kind=realType), dimension(2,2) :: invJac
          real(kind=realType) :: u_old, v_old, du, dv
          real(kind=realType) :: dx, dy, dz, update
          integer(kind=intType) :: ll

          ! Local parameters used in the Newton algorithm.
          integer(kind=intType), parameter :: iterMax   = 15
          real(kind=realType),   parameter :: thresConv = 1.e-12_realType

          ! EXECUTION
          
          ! Initialize u and v to 0.5 and determine the
          ! corresponding coordinates on the face, which is the
          ! centroid.

          u  = half
          v  = half

          ! Newton loop to determine the point on the surface,
          ! which minimizes the distance to the given coordinate.

          NewtonQuads: do ll=1,iterMax

             ! Store previous parametric coordinates
             u_old = u
             v_old = v

             call quadProjResidual(x1, x2, x3, x4, x, &
                                   u, v, &
                                   residual, invJac)

             ! Compute the updates (remember to flip signs)
             call dotProd(invJac(1,:), residual, du)
             call dotProd(invJac(2,:), residual, dv)
             du = -du
             dv = -dv

             ! Determine the new parameter values uu and vv. These
             ! are limited to 0 <= (uu,vv) <= 1.

             u = u + du
             u = min(one,max(zero,u))

             v = v + dv
             v = min(one,max(zero,v))

             ! Update error metric (after the cropping)
             du = u - u_old
             dv = v - v_old
             update = du*du + dv*dv

             ! Exit the loop if the update of the parametric
             ! weights is below the threshold
             if(sqrt(update) <= thresConv) exit NewtonQuads

          enddo NewtonQuads

          ! Call projection one more time for the updated value of u and v
          call quadProjOutput(x1, x2, x3, x4, &
                              u, v, &
                              xf)

          ! Compute the distance squared between the given
          ! coordinate and the point xf.

          dx = x(1) - xf(1)
          dy = x(2) - xf(2)
          dz = x(3) - xf(3)

          val = dx*dx + dy*dy + dz*dz

        end subroutine quadProjection

        subroutine quadProjResidual(x1, x2, x3, x4, x, &
                                    u, v, &
                                    residual, invJac)

          ! This subroutine computes the error that should be reduced
          ! by the iterative process in the subroutine quadProjection.
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
          ! u: real -> Parametric coordinate of the projected point on the quad element.
          !
          ! v: real -> Parametric coordinate of the projected point on the quad element.
          !
          ! OUTPUTS:
          !
          ! du: real -> Change in parametric coordinate due to solving process.
          !
          ! dv: real -> Change in parametric coordinate due to solving process.
          !
          ! error: real(2) -> gradient of the distance squared between the point and its
          !                   projection, which should be zero for minimum distance.
          !                   This is the residual of the equations we want to solve.
          !
          ! xf: real(3) -> Coordinates (X,Y,Z) of the projected point.
          !
          ! Ney Secco - 2016-10
          !
          ! If you write d2 (the distance between the projected and the original
          ! point squared) in terms of u and v you will get:
          ! dist2 = A(T)*UV
          ! where A and UV are given in the code below. The u and v terms are only in UV.
          ! This function will give the updates for u and v based on Broyden's method to
          ! solve grad(dist2) = 0. In other words, we will compute:
          ! delta(u,v) = -inv(Jac(grad(dist2)))*grad(dist2).
          
          implicit none

          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
          real(kind=realType), dimension(3), intent(in) :: x
          real(kind=realType), intent(in) :: u, v

          ! Output variables
          real(kind=realType), dimension(2), intent(out) :: residual
          real(kind=realType), dimension(2,2), intent(out) :: invJac

          ! Working variables
          real(kind=realType), dimension(3) :: x10, x21, x41, x3142
          real(kind=realType), dimension(9) :: A
          real(kind=realType), dimension(2,9) :: gradUV
          real(kind=realType), dimension(2) :: gradDist2
          real(kind=realType), dimension(2,2) :: Jac
          real(kind=realType) :: dotResult

          ! EXECUTION

          ! Determine auxiliary vectors
          x10 = x1 - x
          x21 = x2 - x1
          x41 = x4 - x1
          x3142 = x3 - x1 - x21 - x41

          ! Determine vector of coefficients (A)

          call dotProd(x10, x10, dotResult)
          A(1) = dotResult

          call dotProd(x10, x21, dotResult)
          A(2) = 2*dotResult

          call dotProd(x10, x41, dotResult)
          A(3) = 2*dotResult

          call dotProd(x21, x21, dotResult)
          A(4) = dotResult

          call dotProd(x10, x3142, dotResult)
          A(5) = 2*dotResult
          call dotProd(x21, x41, dotResult)
          A(5) = A(5) + 2*dotResult

          call dotProd(x41, x41, dotResult)
          A(6) = dotResult

          call dotProd(x21, x3142, dotResult)
          A(7) = 2*dotResult

          call dotProd(x41, x3142, dotResult)
          A(8) = 2*dotResult

          call dotProd(x3142, x3142, dotResult)
          A(9) = dotResult

          ! Assemble the gradient of the vector of independent variables (grad(X))
          ! gradUV(1,:) are derivatives with respect to u
          ! gradUV(2,:) are derivatives with respect to v

          gradUV(1,1) = 0.0
          gradUV(1,2) = 1.0
          gradUV(1,3) = 0.0
          gradUV(1,4) = 2.0*u
          gradUV(1,5) = v
          gradUV(1,6) = 0.0
          gradUV(1,7) = 2.0*u*v
          gradUV(1,8) = v*v
          gradUV(1,9) = 2.0*u*v*v

          gradUV(2,1) = 0.0
          gradUV(2,2) = 0.0
          gradUV(2,3) = 1.0
          gradUV(2,4) = 0.0
          gradUV(2,5) = u
          gradUV(2,6) = 2.0*v
          gradUV(2,7) = u*u
          gradUV(2,8) = 2.0*u*v
          gradUV(2,9) = 2.0*u*u*v

          ! Now compute the gradient of the objective function (grad(dist2))
          call dotProd(A, gradUV(1,:), gradDist2(1))
          call dotProd(A, gradUV(2,:), gradDist2(2))

          ! The residual is the objective function gradient
          residual = gradDist2

          ! Now assemble the Jacobian of the function we are solving (Jac(grad(dist2)))
          Jac(1,1) = 2.0*(A(4) + v*A(7) + v*v*A(9))
          Jac(1,2) = A(5) + 2.0*(u*A(7) + v*A(8) + 2*u*v*A(9))
          Jac(2,1) = Jac(1,2)
          Jac(2,2) = 2.0*(A(6) + u*A(8) + u*u*A(9))

          ! Invert the Jacobian
          call invert2x2(Jac,invJac)

        end subroutine quadProjResidual

        !===============================================================

        subroutine quadProjOutput(x1, x2, x3, x4, &
                                  u, v, &
                                  xf)

          ! This subroutine computes the error that should be reduced
          ! by the iterative process in the subroutine quadProjection.
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
          ! u: real -> Parametric coordinate of the projected point on the quad element.
          !
          ! v: real -> Parametric coordinate of the projected point on the quad element.
          !
          ! OUTPUTS:
          !
          ! du: real -> Change in parametric coordinate due to solving process.
          !
          ! dv: real -> Change in parametric coordinate due to solving process.
          !
          ! error: real(2) -> gradient of the distance squared between the point and its
          !                   projection, which should be zero for minimum distance.
          !                   This is the residual of the equations we want to solve.
          !
          ! xf: real(3) -> Coordinates (X,Y,Z) of the projected point.
          !
          ! Ney Secco - 2016-10
          !
          ! If you write d2 (the distance between the projected and the original
          ! point squared) in terms of u and v you will get:
          ! dist2 = A(T)*UV
          ! where A and UV are given in the code below. The u and v terms are only in UV.
          ! This function will give the updates for u and v based on Broyden's method to
          ! solve grad(dist2) = 0. In other words, we will compute:
          ! delta(u,v) = -inv(Jac(grad(dist2)))*grad(dist2).
          
          implicit none

          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, x4
          real(kind=realType), intent(in) :: u, v

          ! Output variables
          real(kind=realType), dimension(3), intent(out) :: xf

          ! Working variables
          real(kind=realType), dimension(3) :: x21, x41, x3142

          ! EXECUTION

          ! Determine auxiliary vectors
          x21 = x2 - x1
          x41 = x4 - x1
          x3142 = x3 - x1 - x21 - x41

          ! Compute guess for projection point
          xf = x1 + u*x21 + v*x41 + u*v*x3142

        end subroutine quadProjOutput

        !===============================================================
        !===============================================================
        !===============================================================
        ! WEIGHT FUNCTIONS


        subroutine triaWeights(u,v,weight)

          ! This subroutine compute the interpolation weights based
          ! on parametric coordinates.

          ! DECLARATIONS

          ! Inputs variables
          real(kind=realType), intent(in) :: u,v

          ! Output variables
          real(kind=realType), dimension(8), intent(out) :: weight          

          ! EXECUTION
          
          ! Initialize weight vector
          weight = 0.0

          ! Update values according to the current element type
          weight(1) = one - u - v
          weight(2) =       u 
          weight(3) =           v

        end subroutine triaWeights

        !===============================================================

        subroutine quadWeights(u,v,weight)

          ! This subroutine compute the interpolation weights based
          ! on parametric coordinates.

          ! DECLARATIONS

          ! Inputs variables
          real(kind=realType), intent(in) :: u,v

          ! Output variables
          real(kind=realType), dimension(8), intent(out) :: weight          

          ! EXECUTION
          
          ! Initialize weight vector
          weight = 0.0

          ! Update values according to the current element type
          weight(1) = (one - u)*(one - v)
          weight(2) =        u *(one - v)
          weight(3) =        u *       v
          weight(4) = (one - u)*       v

        end subroutine quadWeights

        !===============================================================
        !===============================================================
        !===============================================================

        subroutine computeNodalNormals(nCoor, nTria, nQuads, coor, triaConn, &
                                       quadsConn, nodalNormals)
!
!       ****************************************************************
!       *                                                              *
!       * This routine computes the nodal normals for each node by     *
!       * obtaining the normals for each panel face then interpolating *
!       * the data to the nodes.                                       *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of points for which the element must be        *
!       *        determined.                                           *
!       * nTria: Number of triangle elements.                          *
!       * nQuads: Number of quad elements.                             *
!       * coor:  The coordinates of these points.                      *
!       * triaConn: Connectivity information for the triangle elments. *
!       * quadsConn: Connectivity information for the quad elments.    *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * nodalNormals: The interpolated normals at each coordinate    *
!       *               node.                                          *
!       *                                                              *
!       ****************************************************************
!
!       John Jasa

        implicit none

!
!       Subroutine arguments.
!
        ! Input
        integer(kind=intType), intent(in) :: nCoor, nTria, nQuads
        real(kind=realType), dimension(3,nCoor), intent(in) :: coor
        integer(kind=intType), dimension(3,nTria), intent(in) :: triaConn
        integer(kind=intType), dimension(4,nQuads), intent(in) :: quadsConn

        ! Output
        real(kind=realType), dimension(3,nCoor), intent(out) :: nodalNormals

        ! Working
        integer(kind=intType), dimension(nCoor) :: connect_count
        real(kind=realType) :: normal1(3), normal2(3), normal3(3), normal4(3)
        integer(kind=intType) :: i, ind1, ind2, ind3, ind4
        real(kind=realType) :: x1(3), x2(3), x3(3), x4(3)
        real(kind=realType) :: x12(3), x23(3), x13(3), x24(3)
        real(kind=realType) :: dotResult


        !===============================================================

        ! Initialize cumulative variables
        nodalNormals = 0.0
        connect_count = 0

        ! Loop over triangle connectivities
        do i = 1,nTria

          ! Get the indices for each node of the triangle element
          ind1 = triaConn(1, i)
          ind2 = triaConn(2, i)
          ind3 = triaConn(3, i)

          ! Get the coordinates for each node of the triangle element
          x1 = coor(:, ind1)
          x2 = coor(:, ind2)
          x3 = coor(:, ind3)

          ! Compute relative vectors for the triangle sides
          x12 = x1 - x2
          x23 = x2 - x3

          ! Take the cross-product of these vectors to obtain the normal vector
          call crossProd(x12, x23, normal1)

          ! Normalize this normal vector
          call dotProd(normal1, normal1, dotResult)
          normal1 = normal1 / sqrt(dotResult)

          ! Add the contribution of this normal to the nodalNormals array
          nodalNormals(:, ind1) = nodalNormals(:, ind1) + normal1
          nodalNormals(:, ind2) = nodalNormals(:, ind2) + normal1
          nodalNormals(:, ind3) = nodalNormals(:, ind3) + normal1

          ! Add connectivity information to the connect_count array so we
          ! know how many edges are connected to a node.
          ! We divide the nodalNormals vector by this number to obtain the
          ! averaged nodal normal.
          connect_count(ind1) = connect_count(ind1) + 1
          connect_count(ind2) = connect_count(ind2) + 1
          connect_count(ind3) = connect_count(ind3) + 1

        end do

        ! Loop over quad connectivities
        do i = 1,nQuads

          ! Get the indices for each node of the quad element
          ind1 = quadsConn(1, i)
          ind2 = quadsConn(2, i)
          ind3 = quadsConn(3, i)
          ind4 = quadsConn(4, i)

          ! Get the coordinates for each node of the quad element
          x1 = coor(:, ind1)
          x2 = coor(:, ind2)
          x3 = coor(:, ind3)
          x4 = coor(:, ind4)

          ! Compute relative vectors for the quad sides
          x13 = x3 - x1
          x24 = x4 - x2

          ! Take the cross-product of these vectors to obtain the normal vectors
          ! Normalize these normal vectors
          call crossProd(x13, x24, normal1)
          call dotProd(normal1, normal1, dotResult)
          normal1 = normal1 / sqrt(dotResult)

          ! Add the contribution of this normal to the nodalNormals array
          nodalNormals(:, ind1) = nodalNormals(:, ind1) + normal1
          nodalNormals(:, ind2) = nodalNormals(:, ind2) + normal1
          nodalNormals(:, ind3) = nodalNormals(:, ind3) + normal1
          nodalNormals(:, ind4) = nodalNormals(:, ind4) + normal1

          ! Add connectivity information to the connect_count array so we
          ! know how many edges are connected to a node.
          ! We divide the nodalNormals vector by this number to obtain the
          ! averaged nodal normal.
          connect_count(ind1) = connect_count(ind1) + 1
          connect_count(ind2) = connect_count(ind2) + 1
          connect_count(ind3) = connect_count(ind3) + 1
          connect_count(ind4) = connect_count(ind4) + 1

        end do

        do i = 1,3
          ! Divide the nodal normals by the number of edges that each node has
          nodalNormals(i, :) = nodalNormals(i, :) / connect_count
        end do

        do i=1,nCoor
          ! Normalize these new averaged nodal normals
          normal1 = nodalNormals(:, i)
          call dotProd(normal1, normal1, dotResult)
          nodalNormals(:, i) = normal1 / sqrt(dotResult)
        end do

        end subroutine computeNodalNormals

        !===============================================================
        !===============================================================
        !===============================================================

        subroutine crossProd(A, B, C)
!
!       ****************************************************************
!       *       Obtain the cross-product of two vectors, A and B.      *
!       ****************************************************************
!
          implicit none

          real(kind=realType), intent(in) :: A(3), B(3)
          real(kind=realType), intent(out) :: C(3)

          C(1) = A(2) * B(3) - A(3) * B(2)
          C(2) = A(3) * B(1) - A(1) * B(3)
          C(3) = A(1) * B(2) - A(2) * B(1)

        end subroutine crossProd

        !===============================================================

        subroutine dotProd(a,b,c)

          ! This subroutine simply computes a dot product c = a.dot.b
          ! This function is redundant with the one defined in
          ! Utilities.F90, but I decided to do it in order to
          ! keep ADT independent of other pySurf modules.
          !
          ! Ney Secco 2016-10

          implicit none
          
          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(:), intent(in) :: a,b

          ! Output variables
          real(kind=realType), intent(out) :: c

          ! Working variables
          integer(kind=intType) :: ii

          ! EXECUTION
          c = 0.0
          do ii=1,size(a)
             c = c + a(ii)*b(ii)
          end do

        end subroutine dotProd

        !===============================================================

        subroutine invert2x2(A,Ainv)

          ! This subroutine inverts a 2x2 matrix
          !
          ! Ney Secco 2016-11

          implicit none

          ! DECLARATIONS

          ! Input variables
          real(kind=realType), dimension(2,2), intent(in) :: A

          ! Output variables
          real(kind=realType), dimension(2,2), intent(out) :: Ainv

          ! Working variables
          real(kind=realType) :: detInv

          ! EXECUTION

          ! Compute the determinant
          detInv = 1.0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

          ! Compute inverted matrix
          Ainv(1,1) = A(2,2)*detInv
          Ainv(1,2) = -A(1,2)*detInv
          Ainv(2,1) = -A(2,1)*detInv
          Ainv(2,2) = A(1,1)*detInv

        end subroutine invert2x2


      end module adtProjections
