!
!     ******************************************************************
!     *                                                                *
!     * File:          hypsurf.F90                                     *
!     * Author:        John Jasa and Ney Secco                         *
!     * Starting date: 08-11-2016                                      *
!     * Last modified: 08-11-2016                                      *
!     *                                                                *
!     ******************************************************************
!
!
!     ******************************************************************
!     *                                                                *
!     * Contains subroutines for hyperbolic surface mesh generation.   *
!     *                                                                *
!     ******************************************************************
!
      module hypsurf

      use precision
      implicit none

      contains

        !=================================================================


        !=================================================================

        subroutine computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, K, f, numNodes)

          implicit none

          !f2py intent(in) numNodes, r0, N0, S0, rm1, Sm1, theta, sigmaSplay, numLayers, epsE0
          !f2py intent(in) bc1, bc2, layerIndex
          !f2py intent(out) K, f
          !f2py depends(numNodes) r0, N0, S0, rm1, Sm1, K, f


          integer, intent(in) :: layerIndex, numNodes, numLayers
          real*8, intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
          real*8, intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
          real*8, intent(in) :: sigmaSplay, epsE0
          character*32, intent(in) :: bc1, bc2
          real*8, intent(out) :: K(3*numNodes, 3*numNodes), f(3*numNodes)

          real*8 :: r_curr(3), r_next(3), r_prev(3), d_vec(3), d_vec_rot(3), eye(3, 3)
          integer :: index, i

          ! Initialize arrays
          K(:, :) = 0._8
          f(:) = 0._8
          eye(:, :) = 0._8
          do i=1,3
            eye(i, i) = 1.
          end do

          ! Now loop over each node

          index = 1

          if (bc1 .eq. 'splay') then

            ! Get coordinates

            r_curr = r0(:3)
            r_next = r0(4:6)

            ! Get vector that connects r_next to r_curr
            d_vec = r_next - r_curr

            ! Get marching direction vector (orthogonal to the curve and to the surface normal)
            call cross(N0(:, 1), d_vec, d_vec_rot)

            ! Populate matrix
            K(1, :3) = d_vec_rot
            K(2, :3) = N0(:,index)
            K(3, :3) = d_vec

            f(:3) = [S0(1) * (1-sigmaSplay), 0._8, 0._8]

          else if (bc1 .eq. 'constX') then

            ! Populate matrix
            K(2, 4:6) = [0, -1, 0]
            K(3, 4:6) = [0, 0, -1]
            do i=1,3
              K(i, i) = 1.
            end do

          else if (bc1 .eq. 'constY') then
            ! Populate matrix
            K(1, 4:6) = [-1, 0, 0]
            K(3, 4:6) = [0, 0, -1]
            do i=1,3
              K(i, i) = 1.
            end do

          else if (bc1 .eq. 'constZ') then
            ! Populate matrix
            K(1, 4:6) = [-1, 0, 0]
            K(2, 4:6) = [0, -1, 0]
            do i=1,3
              K(i, i) = 1.
            end do

          else if (bc1(:5) .eq. 'curve') then

            ! Populate matrix
            do i=1,3
              K(i, i) = 1.
            end do
            f(:3) = S0(1) * N0(:,1)

          else

            ! Call assembly routine
            call matrixBuilder(index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, numNodes, K, f)

          end if

          do index=2,3!numNodes-1

            ! Call assembly routine
            call matrixBuilder(index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, numNodes, K, f)

          end do

          index = numNodes

          if (bc2 .eq. 'continuous') then
            ! Populate matrix (use same displacements of first node)
            K(3*(index-1)+1:, 3*(index-1)+1:) = eye
            K(3*(index-1)+1:, :3) = -eye

          else if (bc2 .eq. 'splay') then

            ! Get coordinates

            r_curr = r0(3*(index-1)+1:)
            r_prev = r0(3*(index-2)+1:3*(index-2)+3)

            ! Get vector that connects r_next to r_curr
            d_vec = r_curr - r_prev

            ! Get marching direction vector (orthogonal to the curve and to the surface normal)
            call cross(N0(:, index), d_vec, d_vec_rot)

            ! Populate matrix
            K(3*index-2, 3*index-2:) = d_vec_rot
            K(3*index-1, 3*index-2:) = N0(:,index)
            K(3*index-0, 3*index-2:) = d_vec

            f(3*(index-1)+1:3*index) = [S0(1) * (1-sigmaSplay), 0._8, 0._8]

          else if (bc2 .eq. 'constX') then

            ! Populate matrix
            K(3*index-1, 3*(index-2)+1:3*(index-1)+1) = [0, -1, 0]
            K(3*index-0, 3*(index-2)+1:3*(index-1)+1) = [0, 0, -1]
            do i=3*index-2, 3*index
              K(i, i) = 1.
            end do

          else if (bc2 .eq. 'constY') then
            ! Populate matrix
            K(3*index-2, 3*(index-2)+1:3*(index-1)+1) = [-1, 0, 0]
            K(3*index-0, 3*(index-2)+1:3*(index-1)+1) = [0, 0, -1]
            do i=3*index-2, 3*index
              K(i, i) = 1.
            end do

          else if (bc2 .eq. 'constZ') then
            ! Populate matrix
            K(3*index-2, 3*(index-2)+1:3*(index-1)+1) = [-1, 0, 0]
            K(3*index-1, 3*(index-2)+1:3*(index-1)+1) = [0, -1, 0]
            do i=3*index-2, 3*index
              K(i, i) = 1.
            end do

          else if (bc2(:5) .eq. 'curve') then

            ! Populate matrix
            do i=3*index-2, 3*index
              K(i, i) = 1.
            end do
            f(3*index-2:) = S0(index) * N0(:, index)

          else

            ! Call assembly routine
            call matrixBuilder(index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, index, K, f)

          end if


          end subroutine computeMatrices



          subroutine matrixBuilder(curr_index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, numNodes, K, f)

        implicit none

        integer, intent(in) :: curr_index, numNodes, numLayers, layerIndex
        real*8, intent(in) :: r0(3*numNodes), rm1(3*numNodes), N0(3, numNodes), Sm1(numNodes), epsE0, S0(numNodes), theta
        character*32, intent(in) :: bc1, bc2
        real*8, intent(out) :: K(3*numNodes, 3*numNodes), f(3*numNodes)

        real*8 :: r_curr(3), r_next(3), d_vec(3), d_vec_rot(3)
        real*8 :: r0_xi(3), pi, angle, giveAngle, B0(3, 3), B0inv(3, 3)
        real*8 :: r0_xi_n(3), point(3), neigh1_point(3), neigh2_point(3)
        real*8 :: r0_eta_n(3), A0(3, 3), r0_eta(3), norm
        real*8 :: dSensor, dnum, dden, eye(3, 3), B0invg(3), C0(3, 3)
        real*8 :: epsE, epsI, De(3)
        integer :: index, i, neighbor1_index, neighbor2_index


        pi = 3.1415926535897932384626
        do i=1,3
          eye(i, i) = 1.
        end do

        if (curr_index .eq. 1) then  ! forward case

          if ('bc1' .ne. 'continuous') then

            neighbor1_index = 2
            neighbor2_index = 3

            ! Using forward differencing for xi = 1
            r0_xi = 0.5 * (-3 * r0(:3) + 4 * r0(4:6) - r0(7:9))

            angle = pi

          else

            neighbor1_index = numNodes - 1
            neighbor2_index = 2

            ! Using central differencing for zeta = 2:numNodes-1
            r0_xi = 0.5*(r0(4:6) - r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3))

            ! Compute the local grid angle based on the neighbors
            angle = giveAngle(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3), &
              r0(:3), r0(4:6), N0(:, 1))

          end if


        else if (curr_index == numNodes) then  ! backward case
          neighbor1_index = curr_index - 1
          neighbor2_index = curr_index - 2

          ! Using backward differencing for xi = numNodes
          r0_xi = 0.5*(3*r0(3*(curr_index-1)+1:3*(curr_index)+3) -&
          4*r0(3*(neighbor1_index-1)+1:3*(neighbor1_index)+3) +&
          r0(3*(neighbor2_index-1)+1:3*(neighbor2_index)+3))

          angle = pi

        else ! central case
          neighbor1_index = curr_index - 1
          neighbor2_index = curr_index + 1

          neigh2_point = r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)
          neigh1_point = r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)


          ! Using central differencing for zeta = 2:numNodes-1
          r0_xi = 0.5*(neigh2_point - neigh1_point)

          ! Compute the local grid angle based on the neighbors
          angle = giveAngle(neigh1_point, &
                            r0(3*(curr_index-1)+1:3*(curr_index-1)+3), &
                            neigh2_point, &
                            N0(:,curr_index))

        end if


        call cross(N0(:,curr_index), r0_xi, r0_xi_n)

        ! Assemble B0 matrix
        B0(1, :) = r0_xi
        B0(2, :) = r0_xi_n
        B0(3, :) = N0(:, curr_index)

        ! Invert B0
        call matinv3(B0, B0inv)

        ! Compute eta derivatives
        r0_eta = matmul(B0inv, [0._8, Sm1(curr_index), 0._8])
        call cross(N0(:,curr_index), r0_eta, r0_eta_n)

        ! Assemble A0 matrix
        A0(1, :) = r0_eta
        A0(2, :) = r0_eta_n
        A0(3, :) = 0.


        ! Compute grid distribution sensor (Eq. 6.8a)
        dnum = norm(rm1(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3)) + &
        norm(rm1(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3))
        dden = norm(r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3)) + &
        norm(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3))
        dSensor = dnum / dden

        ! Sharp convex corner detection
        if (angle .lt. 70.*pi/180.) then ! Corner detected

          ! Populate matrix with Eq 8.3
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = -eye
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = 2*eye
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = -eye
          f(3*(curr_index-1)+1:3*(curr_index-1)+3) = [0., 0., 0.]

        else

          ! Compute C0 = B0inv*A0
          C0 = matmul(B0inv, A0)

          ! Compute smoothing coefficients
          call dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0, epsE, epsI)

          ! Compute RHS components
          B0invg = matmul(B0inv, [0._8, S0(curr_index), 0._8])

          if (curr_index .eq. 1) then

            if (bc1 .ne. 'continuous') then ! forwards
              De = epsE*(r0(3*(curr_index-1)+1:3*(curr_index-1)+3) - &
              2*r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) + &
              r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3))

              ! Populate matrix
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = &
              -0.5*(1+theta)*C0 - epsI*eye
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = &
              2*(1+theta)*C0 + 2*epsI*eye
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = &
              -1.5*(1+theta)*C0 + (1-epsI)*eye
              f(3*(curr_index-1)+1:3*(curr_index-1)+3) = B0invg + De

            else

              De = epsE*(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) - &
              2*r0(3*(curr_index-1)+1:3*(curr_index-1)+3) + &
              r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3))

              ! Populate matrix
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = &
              -0.5*(1+theta)*C0 - epsI*eye
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = &
              (1 + 2*epsI)*eye
              K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = &
              0.5*(1+theta)*C0 - epsI*eye
              f(3*(curr_index-1)+1:3*(curr_index-1)+3) = B0invg + De
            end if

          else if (curr_index == numNodes) then  ! backwards
            De = epsE*(r0(3*(curr_index-1)+1:3*(curr_index-1)+3) - &
            2*r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) + &
            r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3))

            ! Populate matrix
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = &
            0.5*(1+theta)*C0 - epsI*eye
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = &
            -2*(1+theta)*C0 + 2*epsI*eye
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = &
            1.5*(1+theta)*C0 + (1-epsI)*eye
            f(3*(curr_index-1)+1:3*(curr_index-1)+3) = B0invg + De

          else  ! central
            De = epsE*(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) - &
            2*r0(3*(curr_index-1)+1:3*(curr_index-1)+3) + &
            r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3))

            ! Populate matrix
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = &
            -0.5*(1+theta)*C0 - epsI*eye
            ! K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = &
            ! (1 + 2*epsI)*eye
            ! K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = &
            ! 0.5*(1+theta)*C0 - epsI*eye
            f(3*(curr_index-1)+1:3*(curr_index-1)+3) = B0invg + De
          end if

        end if

        end subroutine matrixBuilder

        subroutine dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0, epsE, epsI)

        implicit none

        integer, intent(in) :: layerIndex, numLayers
        real*8, intent(in) :: dSensor, angle, epsE0, r0_xi(3), r0_eta(3)
        real*8, intent(out) :: epsE, epsI


        real*8 :: norm, Sl, dbar, a, pi, N, R
        integer :: l, ltrans

        pi = 3.14159265358979323846264338

        ! Compute N (Eq. 6.3)
        N = norm(r0_eta) / norm(r0_xi)

        ! Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layerIndex+2

        ltrans = int(3. / 4. * numLayers)
        if (l .le. ltrans) then
          Sl = sqrt(float(l-1)/float(ltrans-1))
        else
          Sl = sqrt(float(ltrans-1)/float(numLayers-1))
        end if

        ! Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max(dSensor**(2/Sl), 0.1_8)

        ! Compute a (Eq 6.12 adjusted for entire angle (angle=2*alpha))
        if (angle .le. pi) then ! Convex corner
          a = 1.0
        else
          a = 1.0 / (1.0 - cos(angle/2)*cos(angle/2))
        end if

        ! Compute auxiliary variable R (Eq. 6.4)
        R = Sl*dbar*a

        ! Compute the dissipation coefficients
        epsE = epsE0*R*N
        epsI = 2*epsE


        end subroutine dissipationCoefficients

        subroutine matinv3(A, B)

        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real*8, intent(in)  :: A(3,3)   !! Matrix
        real*8, intent(out) :: B(3,3)   !! Inverse matrix
        real*8              :: detinv

        !f2py intent(in) A
        !f2py intent(out) B

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

        end subroutine matinv3

        real*8 function giveAngle(r0, r1, r2, N1) result(angle)

        real*8, intent(in), dimension(3) :: r0, r1, r2, N1

        !f2py intent(in) r0, r1, r2, N1
        !f2py intent(out) angle

        real*8, dimension(3) :: dr1, dr2, dr1crossdr2
        real*8 :: dr1dotdr2, dot, norm, arccos_inside, pi, one

        pi = 3.14159265358979323846264338
        one = 1.0

        dr1 = r1 - r0
        dr2 = r2 - r1

        dr1dotdr2 = dot(dr1, dr2) ! dot product
        call cross(dr1, dr2, dr1crossdr2) ! cross product

        ! Compute acute angle and ensure it's <= 1.0
        arccos_inside = dr1dotdr2 / norm(dr1) / norm(dr2)
        angle = dacos(min(arccos_inside, one))

        ! If the cross product points in the same direction of the surface
        ! normal, we have an acute corner
        if (dot(dr1crossdr2, N1) .gt. 0.) then
          angle = pi + angle
        else
          angle = pi - angle
        end if

        return

        end function giveAngle



        real*8 function norm(v)

        implicit none

        real*8, intent(in) :: v(3)
        real*8 :: dot

        norm = dot(v, v) ** 0.5

        return

        end function norm



        real*8 function dot(a, b)

        implicit none

        real*8, intent(in) :: a(3), b(3)

        dot = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

        return

        end function dot



        subroutine cross(A, B, C)

        implicit none

        real*8, intent(in) :: A(3), B(3)
        real*8, intent(out) :: C(3)

        C(1) = A(2) * B(3) - A(3) * B(2)
        C(2) = A(3) * B(1) - A(1) * B(3)
        C(3) = A(1) * B(2) - A(2) * B(1)

        end subroutine cross

      end module
