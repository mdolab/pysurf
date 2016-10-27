!
!     ******************************************************************
!     *                                                                *
!     * File:          hypsurfMain.F90                                 *
!     * Authors:       John Jasa and Ney Secco                         *
!     * Starting date: 08-11-2016                                      *
!     * Last modified: 09-25-2016                                      *
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
      module hypsurfMain

      use precision
      implicit none

      contains

        !=================================================================


        !=================================================================

        subroutine computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex, theta,&
        sigmaSplay, bc1, bc2, numLayers, epsE0, f, numNodes)


        use solveRoutines, only: solve
        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(out) :: f(3*numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), d_vec(3), d_vec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero, rhs(3*numNodes)

        one = 1.
        zero = 0.

        ! Initialize arrays
        K(:, :) = zero
        f(:) = zero
        eye(:, :) = zero
        do i=1,3
          eye(i, i) = one
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

          f(:3) = zero
          f(1) = S0(1) * (1-sigmaSplay)

        else if (bc1 .eq. 'constx') then

          ! Populate matrix
          K(2, 4) = zero
          K(2, 5) = -one
          K(2, 6) = zero
          K(3, 4) = zero
          K(3, 5) = zero
          K(3, 6) = -one
          do i=1,3
            K(i, i) = one
          end do

        else if (bc1 .eq. 'consty') then
          ! Populate matrix
          K(1, 4) = -one
          K(1, 5) = zero
          K(1, 6) = zero
          K(3, 4) = zero
          K(3, 5) = zero
          K(3, 6) = -one
          do i=1,3
            K(i, i) = one
          end do

        else if (bc1 .eq. 'constz') then
          ! Populate matrix
          K(1, 4) = -one
          K(1, 5) = zero
          K(1, 6) = zero
          K(2, 4) = zero
          K(2, 5) = -one
          K(2, 6) = zero
          do i=1,3
            K(i, i) = one
          end do

        else if (bc1(:5) .eq. 'curve') then

          ! Populate matrix
          do i=1,3
            K(i, i) = one
          end do
          f(:3) = S0(1) * N0(:,1)

        else

          ! Call assembly routine
          call matrixBuilder(index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, numNodes, K, f)

        end if

        do index=2,numNodes-1

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

          f(3*(index-1)+1:3*index) = zero
          f(3*(index-1)+1) = S0(index) * (1-sigmaSplay)

        else if (bc2 .eq. 'constx') then

          ! Populate matrix
          K(3*index-0, 3*(index-2)+1) = zero
          K(3*index-0, 3*(index-2)+2) = zero
          K(3*index-0, 3*(index-2)+3) = -one
          K(3*index-1, 3*(index-2)+1) = zero
          K(3*index-1, 3*(index-2)+2) = -one
          K(3*index-1, 3*(index-2)+3) = zero
          do i=3*index-2, 3*index
            K(i, i) = one
          end do

        else if (bc2 .eq. 'consty') then
          ! Populate matrix
          K(3*index-2, 3*(index-2)+1) = -one
          K(3*index-2, 3*(index-2)+2) = zero
          K(3*index-2, 3*(index-2)+3) = zero
          K(3*index-0, 3*(index-2)+1) = zero
          K(3*index-0, 3*(index-2)+2) = zero
          K(3*index-0, 3*(index-2)+3) = -one
          do i=3*index-2, 3*index
            K(i, i) = one
          end do

        else if (bc2 .eq. 'constz') then
          ! Populate matrix
          K(3*index-2, 3*(index-2)+1) = -one
          K(3*index-2, 3*(index-2)+2) = zero
          K(3*index-2, 3*(index-2)+3) = zero
          K(3*index-1, 3*(index-2)+1) = zero
          K(3*index-1, 3*(index-2)+2) = -one
          K(3*index-1, 3*(index-2)+3) = zero
          do i=3*index-2, 3*index
            K(i, i) = one
          end do

        else if (bc2(:5) .eq. 'curve') then

          ! Populate matrix
          do i=3*index-2, 3*index
            K(i, i) = one
          end do
          f(3*index-2:) = S0(index) * N0(:, index)

        else

          ! Call assembly routine
          call matrixBuilder(index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, index, K, f)

        end if

        ! Set other parameters
        n = 3*numNodes ! Problem size
        nrhs = 1 ! number of right hand sides in f
        ldK = n  ! leading dimension of K (should be = n unless we work with submatrices)
        ldf = n  ! leading dimension of f (should be = n unless we work with submatrices)

        ! call dgesv(n, nrhs, K, ldK, ipiv, f, ldf, info)
        rhs = f

        call solve(K, f, rhs, n, ipiv)

        ! Note that this f is rNext when outputted from computeMatrices_main
        f = r0 + f

        end subroutine computeMatrices_main



        subroutine matrixBuilder(curr_index, bc1, bc2, r0, rm1, N0, S0, Sm1, numLayers, epsE0, layerIndex, theta, numNodes, K, f)

        implicit none

        integer(kind=intType), intent(in) :: curr_index, numNodes, numLayers, layerIndex
        real(kind=realType), intent(in) :: r0(3*numNodes), rm1(3*numNodes), N0(3, numNodes)
        real(kind=realType), intent(in) :: Sm1(numNodes), epsE0, S0(numNodes), theta
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(inout) :: K(3*numNodes, 3*numNodes), f(3*numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), d_vec(3), d_vec_rot(3)
        real(kind=realType) :: r0_xi(3), pi, angle, B0(3, 3), B0inv(3, 3)
        real(kind=realType) :: r0_xi_n(3), point(3), neigh1_point(3), neigh2_point(3)
        real(kind=realType) :: r0_eta_n(3), A0(3, 3), r0_eta(3)
        real(kind=realType) :: dSensor, dnum, dden, eye(3, 3), B0invg(3), C0(3, 3)
        real(kind=realType) :: epsE, epsI, De(3), numnorm1, numnorm2, numnorm3, numnorm4
        integer(kind=intType) :: index, i, j, neighbor1_index, neighbor2_index
        real(kind=realType) :: one, zero

        one = 1.
        zero = 0.


        pi = 3.1415926535897932384626
        eye(:, :) = zero
        do i=1,3
          eye(i, i) = one
        end do

        if (curr_index .eq. 1) then  ! forward case

          if (bc1 .ne. 'continuous') then

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
            call giveAngle(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3), &
              r0(:3), r0(4:6), N0(:, 1), angle)

          end if


        else if (curr_index == numNodes) then  ! backward case
          neighbor1_index = curr_index - 1
          neighbor2_index = curr_index - 2

          ! Using backward differencing for xi = numNodes
          r0_xi = 0.5*(3*r0(3*(curr_index-1)+1:3*(curr_index-1)+3) -&
          4*r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) +&
          r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3))

          angle = pi

        else ! central case
          neighbor1_index = curr_index - 1
          neighbor2_index = curr_index + 1

          neigh2_point = r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)
          neigh1_point = r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)


          ! Using central differencing for zeta = 2:numNodes-1
          r0_xi = 0.5*(neigh2_point - neigh1_point)

          ! Compute the local grid angle based on the neighbors
          call giveAngle(neigh1_point, &
                            r0(3*(curr_index-1)+1:3*(curr_index-1)+3), &
                            neigh2_point, &
                            N0(:,curr_index), angle)

        end if


        call cross(N0(:,curr_index), r0_xi, r0_xi_n)

        ! Assemble B0 matrix
        B0(1, :) = r0_xi
        B0(2, :) = r0_xi_n
        B0(3, :) = N0(:, curr_index)

        ! Invert B0
        call matinv3(B0, B0inv)

        ! Compute eta derivatives
        r0_eta = B0inv(:, 2) * Sm1(curr_index)
        call cross(N0(:,curr_index), r0_eta, r0_eta_n)

        ! Assemble A0 matrix
        A0(1, :) = r0_eta
        A0(2, :) = r0_eta_n
        A0(3, :) = zero

        call norm(rm1(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3), numnorm1)
        call norm(rm1(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3), numnorm2)
        call norm(r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3), numnorm3)
        call norm(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3), numnorm4)
        ! Compute grid distribution sensor (Eq. 6.8a)
        dnum = numnorm1 + numnorm2
        dden = numnorm3 + numnorm4
        dSensor = dnum / dden

        ! Sharp convex corner detection
        if (angle .lt. 70.*pi/180.) then ! Corner detected

          ! Populate matrix with Eq 8.3
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = -eye
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = 2*eye
          K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = -eye
          f(3*(curr_index-1)+1:3*(curr_index-1)+3) = zero

        else

          ! Compute C0 = B0inv*A0
          do i=1,3
            do j=1,3
              call dot(B0inv(i, :), A0(:, j), C0(i, j))
            end do
          end do

          ! Compute smoothing coefficients
          call dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0, epsE, epsI)

          ! Compute RHS components
          B0invg = B0inv(:, 2) * S0(curr_index)

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

            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3) = &
            -0.5*(1+theta)*C0 - epsI*eye
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(curr_index-1)+1:3*(curr_index-1)+3) = &
            (1 + 2*epsI)*eye
            K(3*(curr_index-1)+1:3*(curr_index-1)+3,3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3) = &
            0.5*(1+theta)*C0 - epsI*eye
            f(3*(curr_index-1)+1:3*(curr_index-1)+3) = B0invg + De
          end if

        end if

        end subroutine matrixBuilder

        subroutine dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0, epsE, epsI)

        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numLayers
        real(kind=realType), intent(in) :: dSensor, angle, epsE0, r0_xi(3), r0_eta(3)
        real(kind=realType), intent(out) :: epsE, epsI

        real(kind=realType) :: Sl, dbar, a, pi, N, R, normeta, normxi
        integer(kind=intType) :: l, ltrans

        pi = 3.14159265358979323846264338

        ! Compute N (Eq. 6.3)
        call norm(r0_eta, normeta)
        call norm (r0_xi, normxi)

        N = normeta / normxi

        ! Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layerIndex+2

        ltrans = int(3. / 4. * numLayers)
        if (l .le. ltrans) then
          Sl = dsqrt(float(l-1)/float(ltrans-1))
        else
          Sl = dsqrt(float(ltrans-1)/float(numLayers-1))
        end if

        ! Compute adjusted grid distribution sensor (Eq. 6.7)
        dbar = max(dSensor**(2./Sl), 0.1)

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

        subroutine areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, n, S, maxStretch)

        implicit none

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: r0(3*n), d, nuArea
        integer(kind=intType), intent(in) :: numAreaPasses
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(out) :: S(n), maxStretch

        real(kind=realType) :: r0_extrap(3*(2+n))
        real(kind=realType) :: neighborDist(n), norm_1(n), norm_2(n)
        real(kind=realType) :: Sminus, Splus, stretchRatio(n)

        integer(kind=intType) :: index

        ! Extrapolate the end points and copy starting curve
        r0_extrap(:3) = 2*r0(:3) - r0(4:6)
        r0_extrap(4:3*(n+1)) = r0
        r0_extrap(3*(n+1)+1:) = 2*r0(3*(n-1)+1:) - r0(3*(n-2)+1:3*(n-1))

        ! Compute the distance of each node to its neighbors
        do index=1,n
          call norm(r0_extrap(3*index+1:3*index+3) - r0_extrap(3*index-2:3*index), norm_1(index))
          call norm(r0_extrap(3*index+4:3*index+6) - r0_extrap(3*index+1:3*index+3), norm_2(index))
        end do
        neighborDist = 0.5 * (norm_1 + norm_2)

        ! Multiply distances by the step size to get the areas
        S = d * neighborDist

        ! Divide the marching distance and the neighbor distance to get the stretch ratios
        stretchRatio = d / neighborDist

        ! Get the maximum stretch ratio
        maxStretch = -1.e20
        do index=1, n
          if (stretchRatio(index) .gt. maxStretch) then
            maxStretch = stretchRatio(index)
          end if
        end do

        ! Do the requested number of averagings
        do index=1,numAreaPasses
          ! Store previous values
          Splus = S(2)
          Sminus = S(n-1)
          ! Do the averaging for the central nodes
          S(2:n-1) = (1-nuArea)*S(2:n-1) + nuArea/2*(S(:n-2) + S(3:))
          ! Average for the extremum nodes
          S(1) = (1-nuArea)*S(1) + nuArea*Splus
          S(n) = (1-nuArea)*S(n) + nuArea*Sminus
        end do

        ! If we use curve boundary conditions, we need just the marching distance, and not area, for the end nodes
        if (bc1(:5) .eq. 'curve') then
          S(1) = d
        end if
        if (bc2(:5) .eq. 'curve') then
          S(n) = d
        end if

        end subroutine areaFactor


        subroutine smoothing_main(r, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut)

        implicit none

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses, numLayers
        real(kind=realType), intent(in) :: r(3*n)

        real(kind=realType), intent(out) :: rOut(3*n)

        real(kind=realType) :: r_next(3), r_curr(3), r_prev(3), lp, lm, alphaP
        real(kind=realType) :: r_smooth(3*n)

        integer(kind=intType) :: index, index_pass

        ! Compute alpha'
        alphaP = min(alphaP0, alphaP0 * (eta - 3) / numLayers)

        rOut = r

        ! This function does the grid smoothing

        ! Loop over the desired number of smoothing passes
        do index_pass=1,numSmoothingPasses

          ! Copy nodes
          r_smooth = rOut

          ! Smooth every node
          do index=2,n-1

            ! Get coordinates
            r_curr = rOut(3*(index-1)+1:3*(index-1)+3)
            r_next = rOut(3*(index)+1:3*(index)+3)
            r_prev = rOut(3*(index-2)+1:3*(index-2)+3)

            ! Compute distances
            call norm(r_next - r_curr, lp)
            call norm(r_curr - r_prev, lm)

            ! Compute smoothed coordinates
            r_smooth(3*(index-1)+1:3*(index-1)+3) = (1.-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)

          end do

          ! Copy coordinates to allow next pass
          rOut = r_smooth

        end do

        end subroutine smoothing_main

        subroutine qualityCheck(R, layerIndex, numLayers, numNodes, fail, ratios)

        implicit none

        integer(kind=intType), intent(in) :: numLayers, numNodes
        integer(kind=intType), intent(in), optional :: layerIndex
        real(kind=realType), intent(in) :: R(numLayers, 3*numNodes)
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1)
        integer(kind=intType), intent(out) :: fail

        real(kind=realType) :: XYZ(3, numLayers, numNodes), nodalNormals(3, numLayers, numNodes)
        real(kind=realType) :: panelNormals(3, numLayers-1, numNodes-1), norm_vec(3)
        real(kind=realType) :: vec1(3, numLayers-1, numNodes-1), vec2(3, numLayers-1, numNodes-1)
        real(kind=realType) :: vec3(3, numLayers-2, numNodes-2), vec4(3, numLayers-2, numNodes-2)
        real(kind=realType) :: nodalDerivs(3, 2, numLayers, numNodes), det(numLayers, numNodes)
        real(kind=realType) :: normals(3, numLayers-2, numNodes-2), nodalJacs(3, 3, numLayers, numNodes), norm_val
        integer(kind=intType) :: i, j
        real(kind=realType) :: zero

        zero = 0.

        ! Convert the flattened array R into a 3 x numNodes x numLayers array.
        ! numLayers -> number of layers in the marching direction
        ! numNodes -> number of nodes in direction of curve
        do i=1,numNodes
          XYZ(1, :, i) = R(:, 3*(i-1)+1)
          XYZ(2, :, i) = R(:, 3*(i-1)+2)
          XYZ(3, :, i) = R(:, 3*(i-1)+3)
        end do

        ! Setup nodal normals
        nodalNormals(:, :, :) = zero

        ! Get the panel normals from the interior points of the mesh.
        ! Here we take the cross product of the diagonals of each face
        vec1 = XYZ(:, 2:, 2:) - XYZ(:, :numLayers-1, :numNodes-1)
        vec2 = XYZ(:, 2:, :numNodes-1) - XYZ(:, :numLayers-1, 2:)
        do i=1,numNodes-1
          do j=1,numLayers-1
            call cross(vec2(:, j, i), vec1(:, j, i), norm_vec)
            call norm(norm_vec, norm_val)
            panelNormals(:, j, i) = norm_vec / norm_val
          end do
        end do

        ! Set the interior normals using an average of the panel normals
        vec3 = panelNormals(:, 2:, 2:) + panelNormals(:, :numLayers-2, :numNodes-2)
        vec4 = panelNormals(:, 2:, :numNodes-2) + panelNormals(:, :numLayers-2, 2:)
        normals = vec3 + vec4
        do i=2,numNodes-1
          do j=2,numLayers-1
            call norm(normals(:, j-1, i-1), norm_val)
            nodalNormals(:, j, i) = normals(:, j-1, i-1) / norm_val
          end do
        end do

        ! Set the boundary normals
        nodalNormals(:, 2:, 1) = panelNormals(:, :, 1)
        nodalNormals(:, 1, :numNodes-1) = panelNormals(:, 1, :)
        nodalNormals(:, :numLayers-1, numNodes) = panelNormals(:, :, numNodes-1)
        nodalNormals(:, numLayers, 2:) = panelNormals(:, numLayers-1, :)

        ! Setup nodal derivatives
        nodalDerivs(:, :, :, :) = zero

        ! Compute interior derivatives using 2nd order central differencing
        nodalDerivs(:, 1, 2:numLayers-1, 2:numNodes-1) = (XYZ(:, 3:, 2:numNodes-1) - XYZ(:, :numLayers-2, 2:numNodes-1)) / 2.
        nodalDerivs(:, 2, 2:numLayers-1, 2:numNodes-1) = (XYZ(:, 2:numLayers-1, 3:) - XYZ(:, 2:numLayers-1, :numNodes-2)) / 2.

        ! Compute i derivatives using 1st order differencing
        nodalDerivs(:, 1, 1, :) = XYZ(:, 2, :) - XYZ(:, 1, :)
        nodalDerivs(:, 1, numLayers, :) = XYZ(:, numLayers, :) - XYZ(:, numLayers-1, :)

        nodalDerivs(:, 1, 2:numLayers-1, 1) = (XYZ(:, 3:, 1) - XYZ(:, :numLayers-2, 1)) / 2.
        nodalDerivs(:, 1, 2:numLayers-1, numNodes) = (XYZ(:, 3:, numNodes) - XYZ(:, :numLayers-2, numNodes)) / 2.

        ! Compute j derivatives using 1st order differencing
        nodalDerivs(:, 2, :, 1) = XYZ(:, :, 2) - XYZ(:, :, 1)
        nodalDerivs(:, 2, :, numNodes) = XYZ(:, :, numNodes) - XYZ(:, :, numNodes-1)

        nodalDerivs(:, 2, 1, 2:numNodes-1) = (XYZ(:, 1, 3:) - XYZ(:, 1, :numNodes-2)) / 2.
        nodalDerivs(:, 2, numLayers, 2:numNodes-1) = (XYZ(:, numLayers, 3:) - XYZ(:, numLayers, :numNodes-2)) / 2.


        ! Assemble nodal Jacobians
        nodalJacs(:, :, :, :) = zero
        nodalJacs(1, :, :, :) = nodalDerivs(:, 1, :, :)
        nodalJacs(2, :, :, :) = nodalDerivs(:, 2, :, :)
        nodalJacs(3, :, :, :) = nodalNormals(:, :, :)

        ! Compute determinants of Jacobians and find ratio of min to max per face
        ratios(:, :) = zero

        ! Compute the determinants of each nodal Jacobian
        do i=1,numNodes
          do j=1,numLayers
            call m33det(nodalJacs(:, :, j, i), det(j, i))
          end do
        end do

        ! Find the ratio of the minimum valued determinant to the maximum
        ! valued determinant.
        ! This is a measure of quality, with 1 being desirable and anything
        ! less than 0 meaning the mesh is no longer valid.
        do i=1,numNodes-1
          do j=1,numLayers-1
            ratios(j, i) = minval(det(j:j+1, i:i+1)) / maxval(det(j:j+1, i:i+1))
          end do
        end do

        fail = 0
        ! Throw an error and set the failure flag if the mesh is not valid
        do i=1,numNodes-1
          do j=1,numLayers-1
            if (((ratios(j, i) .ne. ratios(j, i)) .or. ratios(j, i) .le. zero) .and. (layerIndex .ge. 1)) then
              print *,'========= FAILURE DETECTED ============'
              fail = 1
            end if
          end do
        end do
        if (fail .eq. 1) then
          print *,"The mesh is not valid after step", layerIndex+1
        end if

        ! Throw a warning if the mesh is low quality
        if ((minval(ratios) .le. .2) .and. (layerIndex .ge. 1)) then
          print *,"The mesh may be low quality after step", layerIndex+1
        end if

        end subroutine qualityCheck


        subroutine findRadius(r, numNodes, radius)

        implicit none

        real(kind=realType), intent(in) :: r(3*numNodes)
        integer(kind=intType), intent(in) :: numNodes

        real(kind=realType), intent(out) :: radius

        real(kind=realType) :: x, y, z
        real(kind=realType) :: minX, maxX, minY, maxY, minZ, maxZ
        integer(kind=intType) :: i

        minX = 1.e20
        minY = 1.e20
        minZ = 1.e20
        maxX = -1.e20
        maxY = -1.e20
        maxZ = -1.e20

        ! Split coordinates and find max and min values
        do i=1,numNodes
          x = r(3*(i-1)+1)
          if (x .gt. maxX) then
            maxX = x
          end if
          if (x .lt. minX) then
            minX = x
          end if

          y = r(3*(i-1)+2)
          if (y .gt. maxY) then
            maxY = y
          end if
          if (y .lt. minY) then
            minY = y
          end if

          z = r(3*(i-1)+3)
          if (z .gt. maxZ) then
            maxZ = z
          end if
          if (z .lt. minZ) then
            minZ = z
          end if
        end do

        ! Find largest radius (we give only half of the largest side to be considered as radius)
        radius = -1.e20
        if (maxX-minX .gt. radius) then
          radius = maxX-minX
        end if
        if (maxY-minY .gt. radius) then
          radius = maxY-minY
        end if
        if (maxZ-minZ .gt. radius) then
          radius = maxZ-minZ
        end if
        radius = radius / 2.

        end subroutine findRadius


        subroutine findRatio(dMax, d0, numLayers, ratioGuess, q)

        implicit none

        real(kind=realType), intent(in) :: dMax, d0, ratioGuess
        integer(kind=intType), intent(in) :: numLayers

        real(kind=realType), intent(out) :: q

        real(kind=realType) :: Rdot, R
        integer(kind=intType) :: nIters

        ! Note that this counter is not the intType that we use for all other integers.
        ! This is done so that Tapenade correctly backwards differentiates this subroutine.
        integer :: i

        ! Extra parameters
        nIters = 200 ! Maximum number of iterations for Newton search

        ! Initialize ratio
        q = ratioGuess

        ! Newton search loop
        do i=1,nIters
          ! Residual function
          R = d0 * (1. - q**(numLayers-1)) - dMax * (1. - q)

          ! Residual derivative
          Rdot = -(numLayers-1)*d0*q**(numLayers-2) + dMax
          ! Update ratio with Newton search
          q = q - R/Rdot
        end do

        ! Check if we got a reasonable value
        if ((q .le. 1) .or. (q .ge. ratioGuess)) then
          print *, ''
          print *, ''
          stop 'Ratio may be too large... Increase number of cells or reduce extension'
        end if

        end subroutine findRatio

        subroutine march_main(rStart, dStart, theta, sigmaSplay, bc1, bc2,&
        epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, R, fail, ratios, majorIndices)

        implicit none

        real(kind=realType), intent(in) :: rStart(3*numNodes),dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, extension, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses

        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1), R(numLayers, 3*numNodes)
        integer(kind=intType), intent(out) :: majorIndices(numLayers)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)
        real(kind=realType) :: rNext_in(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta, min_ratio, ratios_small(3, numNodes-1)
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor
        integer(kind=intType) :: arraySize, nAllocations

        rNext = rStart
        NNext = 0.
        do layerIndex=1,numNodes
          NNext(3, layerIndex) = -1.
        end do

        ! Initialize step size and total marched distance
        d = dStart

        ! Find the characteristic radius of the mesh
        call findRadius(rNext, numNodes, radius)

        ! Find the desired marching distance
        dMax = radius * (extension-1.)

        ! Compute the growth ratio necessary to match this distance
        call findRatio(dMax, dStart, numLayers, ratioGuess, dGrowth)

        ! We need a guess for the first-before-last curve in order to compute the grid distribution sensor
        ! As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
        rm1 = rNext
        N0 = NNext

        !===========================================================

        ! Some functions require the area factors of the first-before-last curve
        ! We will repeat the first curve areas for simplicity.
        ! rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
        ! But we still need to find Sm1

        call areaFactor(rNext, d, nuArea, numAreaPasses, bc1, bc2, numNodes, Sm1, maxStretch)

        R(1, :) = rNext

        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0 = rNext

          ! Compute the new area factor for the desired marching distance
          call areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, maxStretch)

          cfactor = int(maxstretch/cmax) + 1

          ! Constrain the marching distance if the stretching ratio is too high
          dPseudo = d / cFactor

          ! Subiteration
          ! The number of subiterations is the one required to meet the desired marching distance
          do indexSubIter=1,cFactor

            ! Recompute areas with the pseudo-step
            call areaFactor(r0, dPseudo, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, maxStretch)

            ! March using the pseudo-marching distance
            eta = layerIndex+2

            ! Generate matrices of the linear system
            call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex-1, theta,&
            sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, numNodes)

            ! Smooth coordinates
            call smoothing_main(rNext, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed)

            ! Placeholder for projection
            rNext = rSmoothed
            NNext = N0

            Sm1 = S0
            rm1 = r0

            r0 = rNext
            N0 = NNext

          end do

          ! Store grid points
          R(layerIndex+1, :) = rNext

          ! Update step size
          d = d*dGrowth

        end do

        end subroutine march_main





        subroutine matinv3(A, B)

        implicit none

        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real(kind=realType), intent(in)  :: A(3,3)   !! Matrix
        real(kind=realType), intent(out) :: B(3,3)   !! Inverse matrix
        real(kind=realType)              :: detinv

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

        subroutine giveAngle(r0, r1, r2, N1, angle)

        implicit none

        real(kind=realType), intent(in), dimension(3) :: r0, r1, r2, N1

        !f2py intent(in) r0, r1, r2, N1
        !f2py intent(out) angle

        real(kind=realType), dimension(3) :: dr1, dr2, dr1crossdr2
        real(kind=realType) :: dr1dotdr2, arccos_inside, pi, one, angle, tmp
        real(kind=realType) :: normdr1, normdr2, dr1crossdr2dotN1

        pi = 3.14159265358979323846264338
        one = 1.0

        dr1 = r1 - r0
        dr2 = r2 - r1

        call dot(dr1, dr2, dr1dotdr2) ! dot product
        call cross(dr1, dr2, dr1crossdr2) ! cross product

        ! Compute acute angle and ensure it's <= 1.0
        call norm(dr1, normdr1)
        call norm(dr2, normdr2)

        arccos_inside = dr1dotdr2 / normdr1 / normdr2
        angle = dacos(min(arccos_inside, one))

        ! If the cross product points in the same direction of the surface
        ! normal, we have an acute corner
        call dot(dr1crossdr2, N1, dr1crossdr2dotN1)
        if (dr1crossdr2dotN1 .gt. 0.) then
          angle = pi + angle
        else
          angle = pi - angle
        end if

        end subroutine giveAngle

        !============================================================

        subroutine dot(A, B, dot_)

          ! John Jasa - 2016-08

          implicit none

          real(kind=realType), intent(in) :: A(3), B(3)
          real(kind=realType), intent(out) :: dot_

          dot_ = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)

        end subroutine dot

        !============================================================

        subroutine norm(A, norm_)

          ! John Jasa - 2016-08

          implicit none

          real(kind=realType), intent(in) :: A(3)
          real(kind=realType), intent(out) :: norm_

          norm_ = dsqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))

        end subroutine norm


        subroutine cross(A, B, C)

        implicit none

        real(kind=realType), intent(in) :: A(3), B(3)
        real(kind=realType), intent(out) :: C(3)

        C(1) = A(2) * B(3) - A(3) * B(2)
        C(2) = A(3) * B(1) - A(1) * B(3)
        C(3) = A(1) * B(2) - A(2) * B(1)

        end subroutine cross


        subroutine m33det(A, det)

        implicit none

        real(kind=realType), dimension(3, 3), intent(in) :: A
        real(kind=realType), intent(out) :: det

        det =   A(1,1)*A(2,2)*A(3,3)  &
              - A(1,1)*A(2,3)*A(3,2)  &
              - A(1,2)*A(2,1)*A(3,3)  &
              + A(1,2)*A(2,3)*A(3,1)  &
              + A(1,3)*A(2,1)*A(3,2)  &
              - A(1,3)*A(2,2)*A(3,1)


        end subroutine m33det

      end module
