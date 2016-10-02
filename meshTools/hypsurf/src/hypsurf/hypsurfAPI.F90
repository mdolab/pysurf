!
!     ******************************************************************
!     *                                                                *
!     * File:          hypsurfAPI.F90                                  *
!     * Authors:       John Jasa and Ney Secco                         *
!     * Starting date: 09-25-2016                                      *
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
      module hypsurfapi

      use hypsurfMain
      use precision
      implicit none

      real(kind=realType), dimension(:,:), allocatable :: R_initial_march
      real(kind=realType), dimension(:,:), allocatable :: R_smoothed
      real(kind=realType), dimension(:,:), allocatable :: R_final
      real(kind=realType), dimension(:,:), allocatable :: S
      real(kind=realType), dimension(:,:, :), allocatable :: N


      contains

        !=================================================================

        !=======================================================================

        subroutine march(py_projection, rStart, dStart, theta, sigmaSplay, bc1, bc2,&
        epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, R, fail, ratios, majorIndices)

        implicit none

        external py_projection
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
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor, nSubIters
        integer(kind=intType) :: arraySize, nAllocations

        real(kind=realType), dimension(:,:), allocatable :: ext_temp_real
        real(kind=realType), dimension(:,:), allocatable :: ext_temp_real_1dim
        real(kind=realType), dimension(:,:, :), allocatable :: ext_temp_real_N
        real(kind=realType), dimension(:,:), allocatable :: ext_R_initial_march
        real(kind=realType), dimension(:,:), allocatable :: ext_R_smoothed
        real(kind=realType), dimension(:,:), allocatable :: ext_R_final
        real(kind=realType), dimension(:,:), allocatable :: ext_S
        real(kind=realType), dimension(:,:, :), allocatable :: ext_N

        ! Project onto the surface or curve (if applicable)
        call py_projection(rStart, rNext, NNext, numNodes)

        ! Initialize step size and total marched distance
        d = dStart
        dTot = 0

        ! Find the characteristic radius of the mesh
        call findRadius(rNext, numNodes, radius)

        ! Find the desired marching distance
        dMax = radius * (extension-1.)

        ! Compute the growth ratio necessary to match this distance
        call findRatio(dMax, dStart, numLayers, ratioGuess, dGrowth)

        ! Print growth ratio
        print *, 'Growth ratio: ',dGrowth


        ! ! Issue a warning message if the projected points are far from the given points
        ! if max(abs(rNext - rStart)) > 1.0e-5:
        !     warn('The given points (rStart) might not belong to the given surface (surf)\nThese points were projected for the surface mesh generation')


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

        fail = 0
        nSubIters = 0

        ! Initialize extended connectivity arrays for the intersection curves.
        ! These arrays will be cropped to get the actual outputs.
        ! For now we will allocate arrays of size arraySize. We will increase them
        ! if necessary
        arraySize = numLayers
        allocate(ext_R_initial_march(arraySize, 3*numNodes))
        allocate(ext_R_smoothed(arraySize, 3*numNodes))
        allocate(ext_R_final(arraySize, 3*numNodes))
        allocate(ext_S(arraySize, numNodes))
        allocate(ext_N(arraySize, 3, numNodes))

        ! Initialize values
        ext_R_initial_march = 0.
        ext_R_smoothed = 0.
        ext_R_final = 0.
        ext_S = 0.
        ext_N = 0.

        ! Store the initial curves
        ext_R_initial_march(1, :) = rNext
        ext_R_smoothed(1, :) = rNext
        ext_R_final(1, :) = rNext
        ext_S(1, :) = Sm1
        ext_N(1, :, :) = N0
        R(1, :) = rNext
        majorIndices(1) = 1

        ! For now, we have just one allocation
        nAllocations = 1

        print *, ''
        print *, '=================='
        print *, "Marching Progress:"
        print *, '=================='
        print *, ''
        print *, 'Major iteration | Subiters for this step | Lowest mesh quality value | Marching distance'
        print *, '========================================================================================'
        min_ratio = 1.

        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0 = rNext

          ! Compute the new area factor for the desired marching distance
          call areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, maxStretch)

          ! The subiterations will use pseudo marching steps.
          ! If the required marching step is too large, the hyperbolic marching might become
          ! unstable. So we subdivide the current marching step into multiple smaller pseudo-steps

          ! Compute the factor between the current stretching ratio and the allowed one.
          ! If the current stretching ratio is smaller than cMax, the cFactor will be 1.0, and
          ! The pseudo-step will be the same as the desired step.

          cFactor = ceiling(maxStretch/cMax)

          ! Constrain the marching distance if the stretching ratio is too high
          dPseudo = d / cFactor

          print *, layerIndex, '        ', cFactor, '          ', min_ratio, d

          ! Subiteration
          ! The number of subiterations is the one required to meet the desired marching distance
          do indexSubIter=1,cFactor

            ! Increase the number of intersections elements known so far
            nSubIters = nSubIters + 1

            ! Check if we already exceeded the memory allocated so far
            if (nSubIters+1 .gt. size(ext_R_initial_march,1)) then

              ! We need to allocate more memory
              nAllocations = nAllocations + 1

              ! REALLOCATING ext_R_initial_march

              ! Create temporary array and initialize
              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_initial_march

              ! Now move the new allocation back to ext_R_initial_march.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_initial_march)

              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_smoothed

              ! Now move the new allocation back to ext_R_smoothed.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_smoothed)

              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_final

              ! Now move the new allocation back to ext_R_final.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_final)

              allocate(ext_temp_real_1dim(nAllocations*arraySize, numNodes))
              ext_temp_real_1dim = 0.0

              ! Tranfer data from original array
              ext_temp_real_1dim(:(nAllocations-1)*arraySize, :) = ext_S

              ! Now move the new allocation back to ext_S.
              ! ext_temp_real_1dim is deallocated in this process.
              call move_alloc(ext_temp_real_1dim, ext_S)

              allocate(ext_temp_real_N(nAllocations*arraySize, 3, numNodes))
              ext_temp_real_N = 0.0

              ! Tranfer data from original array
              ext_temp_real_N(:(nAllocations-1)*arraySize, :, :) = ext_N

              ! Now move the new allocation back to ext_N.
              ! ext_temp_real_N is deallocated in this process.
              call move_alloc(ext_temp_real_N, ext_N)

            end if

            ! Recompute areas with the pseudo-step
            call areaFactor(r0, dPseudo, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, maxStretch)

            ! March using the pseudo-marching distance
            eta = layerIndex+2

            ! Record the index if this is the coordinate set for a major step
            if (cFactor .eq. indexSubIter) then
              ! Account for python indexing; start at 0
              majorIndices(layerIndex+1) = nSubIters+1
            end if

            ! Generate matrices of the linear system
            call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex-1, theta,&
            sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, numNodes)

            ! Save this initially marched row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_initial_march(nSubIters+1, :) = rNext

            ! Smooth coordinates
            call smoothing_main(rNext, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed)

            ! Save this smoothed row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_smoothed(nSubIters+1, :) = rSmoothed

            call py_projection(rSmoothed, rNext, NNext, numNodes)

            ! Save this final projected row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_final(nSubIters+1, :) = rNext

            ! Save S0 into S for backwards derivative computation
            ext_S(nSubIters+1, :) = S0

            ! Save N0 into N for backwards derivative computation
            ext_N(nSubIters+1, :, :) = N0

            ! Update Sm1 (Store the previous area factors)
            Sm1 = S0

            ! Update rm1
            rm1 = r0

            ! Update r0
            r0 = rNext

            ! Update the Normals with the values computed in the last iteration.
            ! We do this because, the last iteration already projected the new
            ! points to the surface and also computed the normals. So we don't
            ! have to repeat the projection step
            N0 = NNext

          end do

          ! Store grid points
          R(layerIndex+1, :) = rNext

          ! Check quality of the mesh
          if (layerIndex .gt. 2) then
            if (layerIndex .gt. numLayers - 1) then
              call qualityCheck(R(layerIndex-3:layerIndex, :), layerIndex, 4, numNodes, fail, ratios_small)
            else
              call qualityCheck(R(layerIndex-2:layerIndex+1, :), layerIndex, 4, numNodes, fail, ratios_small)
            end if
            min_ratio = minval(ratios_small)
          end if

          ! Compute the total marched distance so far
          dTot = dTot + d

          ! Update step size
          d = d*dGrowth
        end do

        ! Crop extended array
        allocate(R_initial_march(nSubIters+1, 3*numNodes))
        R_initial_march(:,:) = ext_R_initial_march(:nSubIters+1, :)

        ! Crop extended array
        allocate(R_smoothed(nSubIters+1, 3*numNodes))
        R_smoothed(:,:) = ext_R_smoothed(:nSubIters+1, :)

        ! Crop extended array
        allocate(R_final(nSubIters+1, 3*numNodes))
        R_final(:,:) = ext_R_final(:nSubIters+1, :)

        ! Crop extended array
        allocate(S(nSubIters+1, numNodes))
        S(:,:) = ext_S(:nSubIters+1, :)

        ! Crop extended array
        allocate(N(nSubIters+1, 3, numNodes))
        N(:,:, :) = ext_N(:nSubIters+1, :, :)

        call qualityCheck(R, 0, numLayers, numNodes, fail, ratios)

        print *, ''
        print *, 'Final mesh quality minimum:', minval(ratios)
        print *, ''


        end subroutine march

        !=======================================================================

        subroutine march_d(py_projection, rStart, rStartd, dStart, theta, sigmaSplay, bc1, bc2,&
        epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, R, Rd, fail, ratios, majorIndices)

        use hypsurfMain_d, only: computeMatrices_main_d, smoothing_main_d, areafactor_d, findRadius_d, findRatio_d
        implicit none

        external py_projection
        real(kind=realType), intent(in) :: rStart(3*numNodes), rStartd(3*numNodes)
        real(kind=realType), intent(in) :: dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, extension, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses

        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1)
        real(kind=realType), intent(out) :: R(numLayers, 3*numNodes), Rd(numLayers, 3*numNodes)
        integer(kind=intType), intent(out) :: majorIndices(numLayers)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)
        real(kind=realType) :: rNext_in(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor, nSubIters
        integer(kind=intType) :: arraySize, nAllocations

        ! Derivative variables
        real(kind=realType):: r0d(3*numNodes), N0d(3, numNodes), S0d(numNodes)
        real(kind=realType) :: rm1d(3*numNodes), Sm1d(numNodes), rNextd(3*numNodes)
        real(kind=realType) :: rSmoothedd(3*numNodes), NNextd(3, numNodes)
        real(kind=realType) :: r0d_copy(3*numNodes), dd, dGrowthd, dMaxd
        real(kind=realType) :: dStartd, radiusd, dPseudod

        call py_projection(rStart, rNext, NNext, numNodes)

        ! Need py_projection_d here
        rNextd = rStartd
        rNext = rStart
        NNext = 0.
        do layerindex=1,numnodes
          NNext(3, layerindex) = -1.
        end do

        ! Initialize step size and total marched distance
        d = dStart
        dTot = 0

        ! Find the characteristic radius of the mesh
        call findRadius_d(rNext, rNextd, numNodes, radius, radiusd)

        ! Find the desired marching distance
        dMax = radius * (extension-1.)
        dMaxd = radiusd * (extension-1.)

        dStartd = 0.

        ! Compute the growth ratio necessary to match this distance
        call findratio_d(dMax, dMaxd, dStart, dStartd, numLayers, ratioGuess, dGrowth, dGrowthd)

        ! Print growth ratio
        print *, 'Growth ratio: ',dGrowth

        ! We need a guess for the first-before-last curve in order to compute the grid distribution sensor
        ! As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
        rm1 = rNext
        rm1d = rNextd
        N0 = NNext

        !===========================================================

        ! Some functions require the area factors of the first-before-last curve
        ! We will repeat the first curve areas for simplicity.
        ! rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
        ! But we still need to find Sm1
        Sm1d = 0.
        dd = 0.
        call areaFactor_d(rNext, rNextd, d, dd, nuArea, numAreaPasses, bc1, bc2, numNodes, Sm1, Sm1d, maxStretch)

        fail = 0
        nSubIters = 0

        Rd = 0.
        R(1, :) = rNext
        Rd(1, :) = rNextd

        dd = 0.
        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0 = rNext
          r0d = rNextd

          S0d = 0.
          ! Compute the new area factor for the desired marching distance
          call areaFactor_d(r0, r0d, d, dd, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, S0d, maxStretch)

          ! The subiterations will use pseudo marching steps.
          ! If the required marching step is too large, the hyperbolic marching might become
          ! unstable. So we subdivide the current marching step into multiple smaller pseudo-steps

          ! Compute the factor between the current stretching ratio and the allowed one.
          ! If the current stretching ratio is smaller than cMax, the cFactor will be 1.0, and
          ! The pseudo-step will be the same as the desired step.
          cFactor = ceiling(maxStretch/cMax)

          ! Constrain the marching distance if the stretching ratio is too high
          dPseudo = d / cFactor
          dPseudod = dd / cFactor

          ! Subiteration
          ! The number of subiterations is the one required to meet the desired marching distance
          do indexSubIter=1,cFactor

            S0d = 0.
            ! Recompute areas with the pseudo-step
            call areaFactor_d(r0, r0d, dPseudo, dPseudod, nuArea, numAreaPasses, bc1, bc2, numNodes, S0, S0d, maxStretch)

            ! March using the pseudo-marching distance
            eta = layerIndex+2

            rNextd = 0.
            N0d = 0.
            call computeMatrices_main_d(r0, r0d, N0, N0d, S0, S0d, rm1, rm1d, &
            Sm1, Sm1d, layerIndex-1, theta, sigmaSplay, bc1, bc2, &
            numLayers, epsE0, rNext, rNextd, numNodes)

            rSmoothedd = 0.
            call smoothing_main_d(rNext, rNextd, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed, rSmoothedd)

            ! call py_projection(rSmoothed, rNext, NNext, numNodes)
            ! rNext = rSmoothed
            ! rNextd = rSmoothedd
            ! NNext = N0
            ! NNextd = N0d

            ! Update Sm1 (Store the previous area factors)
            Sm1 = S0
            Sm1d = S0d

            ! Update rm1
            rm1 = r0
            rm1d = r0d

            ! Update r0
            r0 = rSmoothed
            r0d = rSmoothedd

            ! Update the normals with the values computed in the last iteration.
            ! We do this because the last iteration already projected the new
            ! points to the surface and also computed the normals, so we don't
            ! have to repeat the projection step
            N0 = NNext
            N0d = 0.

          end do

          ! Store grid points
          R(layerIndex+1, :) = r0
          Rd(layerIndex+1, :) = r0d

          ! Compute the total marched distance so far
          dTot = dTot + d

          ! Update step size
          dd = dd*dGrowth + d*dGrowthd
          d = d*dGrowth
        end do

        end subroutine march_d

        !=======================================================================

        ! subroutine march_b(py_projection, rStart, rStartb, R_initial_march, R_smoothed, &
        ! R_final, N, majorIndices, dStart, theta, sigmaSplay, bc1, bc2, &
        ! epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, &
        ! numSmoothingPasses, numAreaPasses, &
        ! numLayers, numNodes, nSubIters, R, Rb, fail)
        !
        ! use hypsurfMain_b, only: computeMatrices_main_b, smoothing_main_b, areafactor_b
        ! implicit none
        !
        ! external py_projection
        ! real(kind=realType), intent(in) :: rStart(3*numNodes), Rb(numLayers, 3*numNodes)
        ! real(kind=realType), dimension(nSubIters, 3*numNodes), intent(in) :: R_initial_march, R_smoothed, R_final
        ! real(kind=realType), intent(in) :: N(nSubIters, 3, numNodes)
        ! real(kind=realType), intent(in) :: dStart, theta, sigmaSplay
        ! integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses, nSubIters
        ! real(kind=realType), intent(in) :: epsE0, extension, nuArea
        ! character*32, intent(in) :: bc1, bc2
        ! real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        ! integer(kind=intType), intent(in) :: numSmoothingPasses
        ! integer(kind=intType), intent(in) :: majorIndices(numLayers)
        !
        ! integer(kind=intType), intent(out) :: fail
        ! real(kind=realType), intent(out) :: R(numLayers, 3*numNodes), rStartb(3*numNodes)
        !
        ! real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        ! real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        ! real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes), NNextb(3, numNodes)
        ! real(kind=realType) :: rNext_in(3*numNodes)
        ! real(kind=realType) :: R_fullb(nSubIters, 3*numNodes), S_fullb(nSubIters, numNodes)
        ! real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        ! real(kind=realType) :: dPseudo, maxStretch, radius, eta
        ! integer(kind=intType) :: layerIndex, indexSubIter, cFactor
        ! integer(kind=intType) :: arraySize, nAllocations, i, j, majorIndexArray(nSubIters), majorIndex
        ! logical :: isMajor
        !
        ! ! Derivative variables
        ! real(kind=realType):: r0b(3*numNodes), N0b(3, numNodes), S0b(numNodes)
        ! real(kind=realType) :: rm1b(3*numNodes), Sm1b(numNodes), rNextb(3*numNodes)
        ! real(kind=realType) :: rSmoothedb(3*numNodes)
        !
        ! ! Pseudocode outline of backwards code necessary
        !
        ! ! loop backwards over the subiterations:
        ! ! backwards project
        ! ! backwards smooth
        ! ! backwards computematrices
        ! ! this is the seed for the next row
        ! ! if this row is a major iteration, add the weights from Rb
        !
        ! ! Project onto the surface or curve (if applicable)
        ! call py_projection(rStart, rNext, NNext, numNodes)
        !
        ! ! Initialize step size and total marched distance
        ! d = dStart
        ! dTot = 0
        !
        ! ! Find the characteristic radius of the mesh
        ! call findRadius(rNext, numNodes, radius)
        !
        ! ! Find the desired marching distance
        ! dMax = radius * (extension-1.)
        !
        ! ! Compute the growth ratio necessary to match this distance
        ! call findRatio(dMax, dStart, numLayers, ratioGuess, dGrowth)
        !
        ! ! We need a guess for the first-before-last curve in order to compute the grid distribution sensor
        ! ! As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
        ! rm1 = rNext
        ! N0 = NNext
        !
        ! ! Some functions require the area factors of the first-before-last curve
        ! ! We will repeat the first curve areas for simplicity.
        ! ! rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
        ! ! But we still need to find Sm1
        !
        ! call areaFactor(rNext, d, nuArea, numAreaPasses, bc1, bc2, numNodes, Sm1, maxStretch)
        !
        ! ! Start at the last marched layer
        ! majorIndex = numLayers
        !
        ! ! Initial seed for normals is 0. because no external functions have
        ! ! dependencies on the normal values
        ! NNextb = 0.
        ! rNextb = 0.
        !
        ! R_fullb = 0.
        ! S_fullb = 0.
        !
        ! do indexSubIter=nSubIters, 2, -1
        !
        !   rNextb = rNextb + R_fullb(indexSubIter, :)
        !
        !   ! Check to see if this subiteration is a major iteration; if so, set a flag
        !   isMajor = .false.
        !   do i=1,numLayers
        !     if (indexSubIter .eq. majorIndices(i)) then
        !       majorIndex = i
        !       rNextb = rNextb + Rb(majorIndex, :)
        !       isMajor = .true.
        !     end if
        !   end do
        !
        !   print *,indexSubIter, isMajor
        !
        !   ! Obtain the final projected coordinates of this subiteration
        !   rSmoothed = R_smoothed(indexSubIter, :)
        !
        !   ! Backwards project here
        !   ! call py_projection_b(rSmoothed, rSmoothedb, rNext, rNextb, NNext, NNextb, numNodes)
        !
        !   ! Placeholder until we get the projection differentiation
        !   rSmoothedb = rNextb
        !
        !   rNext = R_initial_march(indexSubIter, :)
        !
        !   ! Backwards smooth
        !   call smoothing_main_b(rNext, rNextb, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed, rSmoothedb)
        !
        !   ! Obtain the correct values to use in computeMatrices
        !   r0 = R_final(indexSubIter-1, :)
        !   N0 = N(indexSubIter-1, :, :)
        !   S0 = S(indexSubIter-1, :)
        !   if (indexSubIter .eq. 2) then
        !     rm1 = R_final(indexSubIter-1, :)
        !     Sm1 = S(indexSubIter-1, :)
        !   else
        !     rm1 = R_final(indexSubIter-2, :)
        !     Sm1 = S(indexSubIter-2, :)
        !   end if
        !
        !   ! Backwards computematrices
        !   call computeMatrices_main_b(r0, r0b, N0, N0b, S0, S0b, rm1, rm1b, &
        !   Sm1, Sm1b, layerIndex-1, theta, sigmaSplay, bc1, bc2, &
        !   numLayers, epsE0, rNext, rNextb, numNodes)
        !
        !   R_fullb(indexSubIter-1, :) = R_fullb(indexSubIter-1, :) + r0b
        !   S_fullb(indexSubIter-1, :) = S_fullb(indexSubIter-1, :) + S0b
        !   if (indexSubIter .eq. 2) then
        !     R_fullb(indexSubIter-1, :) = R_fullb(indexSubIter-1, :) + rm1b
        !     S_fullb(indexSubIter-1, :) = S_fullb(indexSubIter-1, :) + Sm1b
        !   else
        !     R_fullb(indexSubIter-2, :) = R_fullb(indexSubIter-2, :) + rm1b
        !     S_fullb(indexSubIter-2, :) = S_fullb(indexSubIter-2, :) + Sm1b
        !   end if
        !
        !   r0 = R_final(indexSubIter-1, :)
        !
        !   ! Compute reverse areaFactor
        !   call areaFactor_b(r0, r0b, d, db, nuarea, numareapasses, bc1, bc2, numNodes, &
        !   S0, S0b, maxstretch)
        !
        !   R_fullb(indexSubIter-1, :) = R_fullb(indexSubIter-1, :) + r0b
        !
        !
        ! end do
        !
        ! end subroutine march_b

        !=======================================================================

        subroutine releaseMemory()

          ! This subroutine just deallocates memory used by the intersection code.
          ! Remember to call this function after you copy the outputs in Python.

          implicit none

          ! Deallocate output variables
          if (allocated(R_initial_march)) then
            deallocate(R_initial_march, R_smoothed, R_final, S, N)
          end if

        end subroutine releaseMemory

        !=================================================================

        subroutine computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, numNodes)

        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(out) :: rNext(3*numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), d_vec(3), d_vec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex-1, theta,&
        sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, numNodes)

        end subroutine computeMatrices

        !=======================================================================

        subroutine computeMatrices_b(r0, r0_b, N0, N0_b, S0, S0_b, rm1, rm1_b, Sm1,&
        Sm1_b, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, rNext_b, numNodes)

        use hypsurfmain_b, only: computematrices_main_b
        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0, rNext_b(3*numNodes)
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(out) :: rNext(3*numNodes)
        real(kind=realType), intent(out) :: r0_b(3*numNodes), N0_b(3, numNodes), S0_b(numNodes)
        real(kind=realType), intent(out) :: rm1_b(3*numNodes), Sm1_b(numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), d_vec(3), d_vec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computematrices_main_b(r0, r0_b, N0, N0_b, S0, S0_b, rm1, rm1_b, Sm1, Sm1_b,&
        layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, rNext_b, numNodes)

        end subroutine computeMatrices_b

        !=================================================================

        subroutine computeMatrices_d(r0, r0_d, N0, N0_d, S0, S0_d, rm1, rm1_d, Sm1,&
        Sm1_d, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, rNext_d, numNodes)

        use hypsurfmain_d, only: computematrices_main_d
        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(out) :: rNext(3*numNodes), rNext_d(3*numNodes)
        real(kind=realType), intent(in) :: r0_d(3*numNodes), N0_d(3, numNodes), S0_d(numNodes)
        real(kind=realType), intent(in) :: rm1_d(3*numNodes), Sm1_d(numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), d_vec(3), d_vec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computematrices_main_d(r0, r0_d, N0, N0_d, S0, S0_d, rm1, rm1_d, Sm1, Sm1_d,&
        layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, rNext, rNext_d, numNodes)

        end subroutine computeMatrices_d

        !=================================================================

        subroutine smoothing(r, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut)

        implicit none

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses, numLayers
        real(kind=realType), intent(in) :: r(3*n)
        real(kind=realType), intent(out) :: rOut(3*n)

        call smoothing_main(r, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut)

        end subroutine smoothing

        !=================================================================

        subroutine smoothing_b(r, rb, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut, rOutb)

        use hypsurfmain_b, only: smoothing_main_b
        implicit none

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses, numLayers
        real(kind=realType), intent(in) :: r(3*n), rOutb(3*n)
        real(kind=realType), intent(out) :: rOut(3*n), rb(3*n)

        call smoothing_main_b(r, rb, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut, rOutb)

        end subroutine smoothing_b

        !=================================================================

        subroutine smoothing_d(r, rd, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut, rOutd)

        use hypsurfmain_d, only: smoothing_main_d
        implicit none

        integer(kind=intType), intent(in) :: n
        real(kind=realType), intent(in) :: eta, alphaP0
        integer(kind=intType), intent(in) :: numSmoothingPasses, numLayers
        real(kind=realType), intent(in) :: r(3*n), rd(3*n)
        real(kind=realType), intent(out) :: rOut(3*n), rOutd(3*n)

        call smoothing_main_d(r, rd, eta, alphaP0, numSmoothingPasses, numLayers, n, rOut, rOutd)

        end subroutine smoothing_d


      end module
