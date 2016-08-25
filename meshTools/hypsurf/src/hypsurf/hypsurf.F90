!
!     ******************************************************************
!     *                                                                *
!     * File:          hypsurf.F90                                     *
!     * Authors:       John Jasa and Ney Secco                         *
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

        subroutine computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, f, numNodes)

        implicit none

        !f2py intent(in) numNodes, r0, N0, S0, rm1, Sm1, theta, sigmaSplay, numLayers, epsE0
        !f2py intent(in) bc1, bc2, layerIndex
        !f2py intent(out) f
        !f2py depends(numNodes) r0, N0, S0, rm1, Sm1, f


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

        ! Initialize arrays
        K(:, :) = 0.
        f(:) = 0.
        eye(:, :) = 0.
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

          f(:3) = [S0(1) * (1-sigmaSplay), 0., 0.]

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

          f(3*(index-1)+1:3*index) = [S0(index) * (1-sigmaSplay), 0., 0.]

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

        ! Set other parameters
        n = 3*numNodes ! Problem size
        nrhs = 1 ! number of right hand sides in f
        ldK = n  ! leading dimension of K (should be = n unless we work with submatrices)
        ldf = n  ! leading dimension of f (should be = n unless we work with submatrices)

        call dgesv(n, nrhs, K, ldK, ipiv, f, ldf, info)

        end subroutine computeMatrices



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
        integer(kind=intType) :: index, i, neighbor1_index, neighbor2_index


        pi = 3.1415926535897932384626
        eye(:, :) = 0.
        do i=1,3
          eye(i, i) = 1.
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
        r0_eta = matmul(B0inv, [0., Sm1(curr_index), 0.])
        call cross(N0(:,curr_index), r0_eta, r0_eta_n)

        ! Assemble A0 matrix
        A0(1, :) = r0_eta
        A0(2, :) = r0_eta_n
        A0(3, :) = 0.

        numnorm1 = norm2(rm1(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3))
        numnorm2 = norm2(rm1(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-rm1(3*(curr_index-1)+1:3*(curr_index-1)+3))
        numnorm3 = norm2(r0(3*(neighbor2_index-1)+1:3*(neighbor2_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3))
        numnorm4 = norm2(r0(3*(neighbor1_index-1)+1:3*(neighbor1_index-1)+3)-r0(3*(curr_index-1)+1:3*(curr_index-1)+3))
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
          f(3*(curr_index-1)+1:3*(curr_index-1)+3) = [0., 0., 0.]

        else

          ! Compute C0 = B0inv*A0
          C0 = matmul(B0inv, A0)

          ! Compute smoothing coefficients
          call dissipationCoefficients(layerIndex, r0_xi, r0_eta, dSensor, angle, numLayers, epsE0, epsE, epsI)

          ! Compute RHS components
          B0invg = matmul(B0inv, [0., S0(curr_index), 0.])

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
            !   PROBLEM IS HERE
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
        normeta = norm2(r0_eta)
        normxi = norm2(r0_xi)

        N = normeta / normxi

        ! Compute Sl (Eq. 6.5) based on a transition l of 3/4 of max
        l = layerIndex+2

        ltrans = int(3. / 4. * numLayers)
        if (l .le. ltrans) then
          Sl = sqrt(float(l-1)/float(ltrans-1))
        else
          Sl = sqrt(float(ltrans-1)/float(numLayers-1))
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

        real(kind=realType) :: r0minus1(3), r0plus1(3), r0_extrap(3*(2+n))
        real(kind=realType) :: R0_extrap_3d(3, 2+n), neighborDist(n), norm_1(n), norm_2(n)
        real(kind=realType) :: Sminus, Splus, stretchRatio(n)

        integer(kind=intType) :: index

        !f2py intent(in) n, r0, d, nuArea, numAreaPasses, bc1, bc2
        !f2py intent(out) S, maxStretch
        !f2py depends(n) r0, S

        ! Extrapolate the end points
        r0minus1 = 2*r0(:3) - r0(4:6)
        r0plus1 = 2*r0(3*(n-1)+1:) - r0(3*(n-2)+1:3*(n-1))

        r0_extrap(:3) = r0minus1
        r0_extrap(4:3*(n+1)) = r0
        r0_extrap(3*(n+1)+1:) = r0plus1

        ! Reshape so we have 3D array
        R0_extrap_3d = RESHAPE(r0_extrap, (/3, 2+n/))

        ! Compute the distance of each node to its neighbors
        norm_1 = norm2(R0_extrap_3d(:, 2:n+1) - R0_extrap_3d(:, :n), 1)
        norm_2 = norm2(R0_extrap_3d(:, 3:) - R0_extrap_3d(:, 2:n+1), 1)
        neighborDist = 0.5 * (norm_1 + norm_2)

        ! Multiply distances by the step size to get the areas
        S = d * neighborDist

        ! Divide the marching distance and the neighbor distance to get the stretch ratios
        stretchRatio = d / neighborDist

        ! Get the maximum stretch ratio
        maxStretch = maxval(stretchRatio)

        ! Do the requested number of averagings
        do index=1,numAreaPasses
          ! Store previous values
          Splus = S(2)
          Sminus = S(n-2)
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


        subroutine smoothing(r, eta, alphaP0, numSmoothingPasses, numLayers, n)

        implicit none

        !f2py intent(in) eta, alphaP0, numSmoothingPasses, numLayers, n
        !f2py intent(inout) r
        !f2py depends(n) r

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses, numLayers
        real(kind=realType), intent(inout) :: r(3*n)

        real(kind=realType) :: r_next(3), r_curr(3), r_prev(3), lp, lm, alphaP
        real(kind=realType) :: r_smooth(3*n)

        integer(kind=intType) :: index, index_pass

        ! This function does the grid smoothing

            ! Loop over the desired number of smoothing passes
            do index_pass=1,numSmoothingPasses

              ! Copy nodes
              r_smooth = r

              ! Smooth every node
              do index=2,n-1

                ! Get coordinates
                r_curr = r(3*(index-1)+1:3*(index-1)+3)
                r_next = r(3*(index)+1:3*(index)+3)
                r_prev = r(3*(index-2)+1:3*(index-2)+3)

                ! Compute distances
                lp = norm2(r_next - r_curr)
                lm = norm2(r_curr - r_prev)

                ! Compute alpha'
                alphaP = min(alphaP0, alphaP0 * (eta - 2) / numLayers)

                ! Compute smoothed coordinates
                r_smooth(3*(index-1)+1:3*(index-1)+3) = (1-alphaP)*r_curr + alphaP*(lm*r_next + lp*r_prev)/(lp + lm)
              end do

              ! Copy coordinates to allow next pass
              r = r_smooth
            end do


        end subroutine smoothing

        subroutine qualityCheck(R, layerIndex, nw1, nw2, fail, ratios)

        implicit none

        integer(kind=intType), intent(in) :: nw1, nw2
        integer(kind=intType), intent(in), optional :: layerIndex
        real(kind=realType), intent(in) :: R(3*nw2, nw1)
        real(kind=realType), intent(out) :: ratios(nw2-1, nw1-1)
        integer(kind=intType), intent(out) :: fail


        real(kind=realType) :: XYZ(3, nw2, nw1), nodalNormals(3, nw2, nw1)
        real(kind=realType) :: panelNormals(3, nw2-1, nw1-1), norm_vec(3)
        real(kind=realType) :: vec1(3, nw2-1, nw1-1), vec2(3, nw2-1, nw1-1)
        real(kind=realType) :: vec3(3, nw2-2, nw1-2), vec4(3, nw2-2, nw1-2)
        real(kind=realType) :: nodalDerivs(3, 2, nw2, nw1), det(nw2, nw1)
        real(kind=realType) :: normals(3, nw2-2, nw1-2), nodalJacs(3, 3, nw2, nw1)
        integer(kind=intType) :: i, j

        ! Convert the flattened array R into a 3 x nw2 x nw1 array.
        ! nw1 -> number of nodes in the direction of marching
        ! nw2 -> number of nodes in direction of curve
        do i=1,nw2
          XYZ(1, i, :) = R(3*(i-1)+1, :)
          XYZ(2, i, :) = R(3*(i-1)+2, :)
          XYZ(3, i, :) = R(3*(i-1)+3, :)
        end do

        ! Setup nodal normals
        nodalNormals(:, :, :) = 0.

        ! Get the panel normals from the interior points of the mesh.
        ! Here we take the cross product of the diagonals of each face
        vec1 = XYZ(:, 2:, 2:) - XYZ(:, :nw2-1, :nw1-1)
        vec2 = XYZ(:, :nw2-1, 2:) - XYZ(:, 2:, :nw1-1)
        do i=1,nw2-1
          do j=1,nw1-1
            call cross(vec2(:, i, j), vec1(:, i, j), norm_vec)
            panelNormals(:, i, j) = norm_vec / norm2(norm_vec)
          end do
        end do

        ! Set the interior normals using an average of the panel normals
        vec3 = panelNormals(:, 2:, 2:) + panelNormals(:, :nw2-2, :nw1-2)
        vec4 = panelNormals(:, :nw2-2, 2:) + panelNormals(:, 2:, :nw1-2)
        normals = vec3 + vec4
        do i=2,nw2-1
          do j=2,nw1-1
            nodalNormals(:, i, j) = normals(:, i-1, j-1) / norm2(normals(:, i-1, j-1))
          end do
        end do

        ! Set the boundary normals
        nodalNormals(:, 1, 2:) = panelNormals(:, 1, :)
        nodalNormals(:, :nw2-1, 1) = panelNormals(:, :, 1)
        nodalNormals(:, nw2, :nw1-1) = panelNormals(:, nw2-1, :)
        nodalNormals(:, 2:, nw1) = panelNormals(:, :, nw1-1)

        ! Setup nodal derivatives
        nodalDerivs(:, :, :, :) = 0.

        ! Compute interior derivatives using 2nd order central differencing
        nodalDerivs(:, 1, 2:nw2-1, 2:nw1-1) = (XYZ(:, 2:nw2-1, 3:) - XYZ(:, 2:nw2-1, :nw1-2)) / 2.
        nodalDerivs(:, 2, 2:nw2-1, 2:nw1-1) = (XYZ(:, 3:, 2:nw1-1) - XYZ(:, :nw2-2, 2:nw1-1)) / 2.

        ! Compute i derivatives using 1st order differencing
        nodalDerivs(:, 1, :, 1) = XYZ(:, :, 2) - XYZ(:, :, 1)
        nodalDerivs(:, 1, :, nw1) = XYZ(:, :, nw1) - XYZ(:, :, nw1-1)

        nodalDerivs(:, 1, 1, 2:nw1-1) = (XYZ(:, 1, 3:) - XYZ(:, 1, :nw1-2)) / 2.
        nodalDerivs(:, 1, nw2, 2:nw1-1) = (XYZ(:, nw2, 3:) - XYZ(:, nw2, :nw1-2)) / 2.

        ! Compute j derivatives using 1st order differencing
        nodalDerivs(:, 2, 1, :) = XYZ(:, 2, :) - XYZ(:, 1, :)
        nodalDerivs(:, 2, nw2, :) = XYZ(:, nw2, :) - XYZ(:, nw2-1, :)

        nodalDerivs(:, 2, 2:nw2-1, 1) = (XYZ(:, 3:, 1) - XYZ(:, :nw2-2, 1)) / 2.
        nodalDerivs(:, 2, 2:nw2-1, nw1) = (XYZ(:, 3:, nw1) - XYZ(:, :nw2-2, nw1)) / 2.


        ! Assemble nodal Jacobians
        nodalJacs(:, :, :, :) = 0.
        nodalJacs(1, :, :, :) = nodalDerivs(:, 1, :, :)
        nodalJacs(2, :, :, :) = nodalDerivs(:, 2, :, :)
        nodalJacs(3, :, :, :) = nodalNormals(:, :, :)

        ! Compute determinants of Jacobians and find ratio of min to max per face
        ratios(:, :) = 0.

        ! Compute the determinants of each nodal Jacobian
        do i=1,nw2
          do j=1,nw1
            call m33det(nodalJacs(:, :, i, j), det(i, j))
          end do
        end do

        ! Find the ratio of the minimum valued determinant to the maximum
        ! valued determinant.
        ! This is a measure of quality, with 1 being desirable and anything
        ! less than 0 meaning the mesh is no longer valid.
        do i=1,nw2-1
          do j=1,nw1-1
            ratios(i, j) = minval(det(i:i+1, j:j+1)) / maxval(det(i:i+1, j:j+1))
          end do
        end do

        fail = 0
        ! Throw an error and set the failure flag if the mesh is not valid
        do i=1,nw2-1
          do j=1,nw1-1
            if ((isnan(ratios(i, j)) .or. ratios(i, j) .le. 0.) .and. (layerIndex .ge. 1)) then
              print *,"The mesh is not valid after step", layerIndex+1
              fail = 1
            end if
          end do
        end do

        ! Throw a warning if the mesh is low quality
        if ((minval(ratios) .le. .2) .and. (layerIndex .ge. 1)) then
          print *,"The mesh may be low quality after step", layerIndex+1
        end if

        end subroutine qualityCheck


        subroutine subIteration(py_projection, r0, N0, S0, rm1, Sm1, layerIndex, theta, sigmaSplay, bc1, bc2,&
        numLayers, epsE0, eta, alphaP0, numSmoothingPasses, numNodes, rNext, NNext)

        implicit none

        external py_projection
        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses

        real(kind=realType), intent(out) :: rNext(3*numNodes), NNext(3, numNodes)


        real(kind=realType) :: rNext_in(3*numNodes)
        real(kind=realType) :: dr(3*numNodes)

        ! Generate matrices of the linear system
        call computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, dr, numNodes)

        ! Update r
        rNext = r0 + dr

        ! Smooth coordinates
        call smoothing(rNext, eta, alphaP0, numSmoothingPasses, numLayers, numNodes)

        rNext_in = rNext

        call py_projection(rNext_in, rNext, NNext, numNodes)

        end subroutine subIteration



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
        real(kind=realType) :: normdr1, normdr2

        pi = 3.14159265358979323846264338
        one = 1.0

        dr1 = r1 - r0
        dr2 = r2 - r1

        dr1dotdr2 = dot_product(dr1, dr2) ! dot product
        call cross(dr1, dr2, dr1crossdr2) ! cross product

        ! Compute acute angle and ensure it's <= 1.0
        normdr1 = norm2(dr1)
        normdr2 = norm2(dr2)

        arccos_inside = dr1dotdr2 / normdr1 / normdr2
        angle = dacos(min(arccos_inside, one))

        ! If the cross product points in the same direction of the surface
        ! normal, we have an acute corner
        if (dot_product(dr1crossdr2, N1) .gt. 0.) then
          angle = pi + angle
        else
          angle = pi - angle
        end if

        end subroutine giveAngle


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
