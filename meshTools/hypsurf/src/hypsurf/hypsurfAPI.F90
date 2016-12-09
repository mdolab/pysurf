!
!     ******************************************************************
!     *                                                                *
!     * File:          hypsurfAPI.F90                                  *
!     * Authors:       John Jasa and Ney Secco                         *
!     * Starting date: 09-25-2016                                      *
!     * Last modified: 10-26-2016                                      *
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
        epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, extension_given,&
        guideIndices, retainSpacing, numSmoothingPasses, numAreaPasses, numGuides,&
        numLayers, numNodes, R, fail, ratios, majorIndices)

        implicit none

        external py_projection
        real(kind=realType), intent(in) :: rStart(3*numNodes),dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, marchParameter, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses
        logical, intent(in) :: extension_given
        integer(kind=intType), intent(in) :: numGuides
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing

        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1), R(numLayers, 3*numNodes)
        integer(kind=intType), intent(out) :: majorIndices(numLayers)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)
        real(kind=realType) :: rNext_in(3*numNodes), rNext_(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta, min_ratio, ratios_small(3, numNodes-1)
        real(kind=realType) :: extension, growthRatio
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor, nSubIters
        integer(kind=intType) :: arraySize, nAllocations

        integer(kind=intType) :: breakNodes(numGuides+2), intervalID, node1, node2
        integer(kind=intType) :: numInInterval, numIntervals

        real(kind=realType), dimension(:), allocatable :: rInterval, arcLength, startArcLength, rNew

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

        ! Really only do this if we're remeshing
        breakNodes(1) = 1
        breaknodes(2:numGuides+1) = guideIndices
        breakNodes(numGuides+2) = numNodes

        !!! Still need to do this!
        ! Remove duplicate entries from breakNodes if they exist.
        ! This will only occur at the beginning or end because we've already
        ! sorted the internal part of breakNodes by sorting guideIndices.
        ! if (breakNodes(1) .eq. breakNodes(2))

        numIntervals = numGuides + 1

        ! Initialize step size and total marched distance
        d = dStart
        dTot = 0

        ! Find the characteristic radius of the mesh
        call findRadius(rNext, numNodes, radius)

        if (extension_given) then
          ! Find the desired marching distance
          ! Here marchParameter = extension
          dMax = radius * (marchParameter-1.)

          ! Compute the growth ratio necessary to match this distance
          call findRatio(dMax, dStart, numLayers, ratioGuess, dGrowth)
        else
          ! Here marchParameter = growthRatio
          dGrowth = marchParameter
        end if

        ! Print growth ratio
        print *, 'Growth ratio: ',dGrowth

        ! Issue a warning message if the projected points are far from the given points
        if (maxval(abs(rNext - rStart)) .gt. 1.0e-5) then
          print *, ''
          print *, 'The given points (rStart) might not belong to the given surface (surf)'
          print *, 'These points were projected for the surface mesh generation'
          print *, ''
        end if


        ! We need a guess for the first-before-last curve in order to compute the grid distribution sensor
        ! As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
        rm1 = rNext
        N0 = NNext

        !===========================================================

        ! Some functions require the area factors of the first-before-last curve
        ! We will repeat the first curve areas for simplicity.
        ! rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
        ! But we still need to find Sm1

        call areaFactor(rNext, d, nuArea, numAreaPasses, bc1, bc2, guideIndices, retainSpacing,&
        numGuides, numNodes, Sm1, maxStretch)

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
        print *, '========================================================================================'
        print *, "                                Marching Progress:"
        print *, ''
        print *, 'Major iteration | Subiters for this step | Lowest mesh quality value | Marching distance'
        print *, '========================================================================================'
        min_ratio = 1.

        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0 = rNext

          ! Compute the new area factor for the desired marching distance
          call areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, guideIndices, retainSpacing,&
          numGuides, numNodes, S0, maxStretch)

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

            ! Recompute areas with the pseudo-step
            call areaFactor(r0, dPseudo, nuArea, numAreaPasses, bc1, bc2, guideIndices, retainSpacing,&
            numGuides, numNodes, S0, maxStretch)

            ! Save S0 into S for backwards derivative computation
            ext_S(nSubIters+1, :) = S0

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

            ! March using the pseudo-marching distance
            eta = layerIndex+2

            ! Record the index if this is the coordinate set for a major step
            if (cFactor .eq. indexSubIter) then
              ! Account for python indexing; start at 0
              majorIndices(layerIndex+1) = nSubIters+1
            end if

            ! Generate matrices of the linear system
            call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex-1, theta,&
            sigmaSplay, bc1, bc2, numLayers, epsE0, guideIndices, retainSpacing, rNext, numNodes, numGuides)

            ! Save this initially marched row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_initial_march(nSubIters+1, :) = rNext

            ! Smooth coordinates
            call smoothing_main(rNext, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed)

            ! Save this smoothed row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_smoothed(nSubIters+1, :) = rSmoothed

            call py_projection(rSmoothed, rNext, NNext, numNodes)

            rNext_(:) = 0.

            if (retainSpacing) then

              ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
              do intervalID=1, numGuides+1

                ! Get indices of the first and last nodes of the interval
                node1 = breakNodes(intervalID)
                node2 = breakNodes(intervalID+1)

                numInInterval = node2 - node1 + 1

                allocate(rInterval(3*numInInterval))
                allocate(arcLength(numInInterval))
                allocate(startArcLength(numInInterval))
                allocate(rNew(3*numInInterval))

                ! Take the slice of the nodal coordinates that corresponds to the
                ! current subinterval for the start curve
                rInterval = rStart(3*(node1-1)+1:3*node2)

                ! Compute the normalized arc-lengths of this interval
                call compute_arc_length(rInterval, numInInterval, startArcLength)

                ! Get the current nodal coordinates
                rInterval = rNext(3*(node1-1)+1:3*node2)

                ! Compute the normalized arc-lengths of this interval
                call compute_arc_length(rInterval, numInInterval, arcLength)

                ! Redistribute the nodes
                call redistribute_nodes_by_arc_length(rInterval, arcLength, startArcLength, numInInterval, rNew)

                rNext_(3*(node1-1)+1:3*node2) = rNew

                deallocate(rInterval, arcLength, startArcLength, rNew)

              end do

              ! Project the remeshed curve back onto the surface
              call py_projection(rNext_, rNext, NNext, numNodes)

            end if

            ! Save this final projected row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_final(nSubIters+1, :) = rNext

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

        subroutine march_d(py_projection, py_projection_d, rStart, rStartd, dStart, theta, sigmaSplay, bc1, bc2,&
        epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, guideIndices, retainSpacing,&
        extension_given, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, numGuides, R, Rd, fail, ratios, majorIndices)

        use hypsurfMain_d, only: computeMatrices_main_d, smoothing_main_d, areafactor_d, findRadius_d, findRatio_d
        implicit none

        external py_projection
        external py_projection_d
        real(kind=realType), intent(in) :: rStart(3*numNodes), rStartd(3*numNodes)
        real(kind=realType), intent(in) :: dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, marchParameter, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses
        logical, intent(in) :: extension_given
        integer(kind=intType), intent(in) :: numGuides
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing

        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1)
        real(kind=realType), intent(out) :: R(numLayers, 3*numNodes), Rd(numLayers, 3*numNodes)
        integer(kind=intType), intent(out) :: majorIndices(numLayers)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)
        real(kind=realType) :: rNext_in(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta, extension
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
        ! call py_projection_d(rStart, rStartd, rNext, rNextd, NNext, NNextd, numNodes)
        rNextd = rStartd
        rNext = rStart
        NNextd(:, :) = 0.

        ! Initialize step size and total marched distance
        d = dStart
        dTot = 0

        if (extension_given) then
          extension = marchParameter
          ! Find the characteristic radius of the mesh
          call findRadius_d(rNext, rNextd, numNodes, radius, radiusd)

          ! Find the desired marching distance
          dMax = radius * (extension-1.)
          dMaxd = radiusd * (extension-1.)

          dStartd = 0.

          ! Compute the growth ratio necessary to match this distance
          call findratio_d(dMax, dMaxd, dStart, dStartd, numLayers, ratioGuess, dGrowth, dGrowthd)
        else
          dGrowth = marchParameter
          dGrowthd = 0.
        end if


        ! We need a guess for the first-before-last curve in order to compute the grid distribution sensor
        ! As we still don't have a "first-before-last curve" yet, we will just repeat the coordinates
        rm1 = rNext
        rm1d = rNextd
        N0 = NNext
        N0d = NNextd

        !===========================================================

        ! Some functions require the area factors of the first-before-last curve
        ! We will repeat the first curve areas for simplicity.
        ! rNext, NNext, rm1 for the first iteration are computed at the beginning of the function.
        ! But we still need to find Sm1
        Sm1d = 0.
        dd = 0.
        call areaFactor_d(rNext, rNextd, d, dd, nuArea, numAreaPasses,&
        bc1, bc2, guideIndices, retainSpacing, numguides, numNodes, Sm1, Sm1d, maxStretch)

        fail = 0
        nSubIters = 0

        Rd = 0.
        Rd(1, :) = rNextd
        R(1, :) = rNext

        dd = 0.
        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0d = rNextd
          r0 = rNext

          S0d = 0.
          ! Compute the new area factor for the desired marching distance
          call areaFactor_d(r0, r0d, d, dd, nuArea, numAreaPasses, bc1, bc2,&
          guideIndices, retainSpacing, numguides, numNodes, S0, S0d, maxStretch)

          ! The subiterations will use pseudo marching steps.
          ! If the required marching step is too large, the hyperbolic marching might become
          ! unstable. So we subdivide the current marching step into multiple smaller pseudo-steps

          ! Compute the factor between the current stretching ratio and the allowed one.
          ! If the current stretching ratio is smaller than cMax, the cFactor will be 1.0, and
          ! The pseudo-step will be the same as the desired step.
          cFactor = ceiling(maxStretch/cMax)

          ! Constrain the marching distance if the stretching ratio is too high
          dPseudod = dd / cFactor
          dPseudo = d / cFactor

          ! Subiteration
          ! The number of subiterations is the one required to meet the desired marching distance
          do indexSubIter=1,cFactor

            S0d = 0.
            ! Recompute areas with the pseudo-step
            call areaFactor_d(r0, r0d, dPseudo, dPseudod, nuArea, numAreaPasses,&
            bc1, bc2, guideIndices, retainSpacing, numguides, numNodes, S0, S0d, maxStretch)

            ! March using the pseudo-marching distance
            eta = layerIndex+2

            rNextd = 0.
            call computeMatrices_main_d(r0, r0d, N0, N0d, S0, S0d, rm1, rm1d, &
            Sm1, Sm1d, layerIndex-1, theta, sigmaSplay, bc1, bc2, &
            numLayers, epsE0, guideIndices, retainSpacing, rNext, rNextd, numNodes, numGuides)

            rSmoothedd = 0.
            call smoothing_main_d(rNext, rNextd, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed, rSmoothedd)

            call py_projection(rSmoothed, rNext, NNext, numNodes)
            call py_projection_d(rSmoothed, rSmoothedd, rNext, rNextd, NNext, NNextd, numNodes)
            rNext = rSmoothed
            rNextd = rSmoothedd
            NNext = N0
            NNextd = N0d

            ! Update Sm1 (Store the previous area factors)
            Sm1d = S0d
            Sm1 = S0

            ! Update rm1
            rm1d = r0d
            rm1 = r0

            ! Update r0
            r0d = rNextd
            r0 = rNext

            ! Update the normals with the values computed in the last iteration.
            ! We do this because the last iteration already projected the new
            ! points to the surface and also computed the normals, so we don't
            ! have to repeat the projection step
            N0 = NNext
            N0d = NNextd

          end do

          ! Store grid points
          Rd(layerIndex+1, :) = rNextd
          R(layerIndex+1, :) = rNext

          ! Compute the total marched distance so far
          dTot = dTot + d

          ! Update step size
          dd = dd*dGrowth + d*dGrowthd
          d = d*dGrowth
        end do

        end subroutine march_d

        !=======================================================================

        subroutine march_b(py_projection, rStart, rStartb, R_initial_march, R_smoothed, &
        R_final, N, majorIndices, dStart, theta, sigmaSplay, &
        &   bc1, bc2, epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax,&
        guideIndices, retainSpacing, extension_given, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, nSubIters, numGuides, R, Rb, fail&
        &   , ratios)

          use hypsurfMain_b, only: computeMatrices_main_b, smoothing_main_b, areaFactor_b, findRadius_b, findRatio_b
          implicit none

          external py_projection
          integer(kind=inttype), intent(in) :: numNodes, numLayers, &
        &   numAreaPasses, nSubIters
          real(kind=realType), dimension(nSubIters, 3*numNodes), intent(in) :: R_initial_march, R_smoothed, R_final
          real(kind=realType), intent(in) :: N(nSubIters, 3, numNodes)
          real(kind=realtype), intent(in) :: rStart(3*numNodes), dStart, theta&
        &   , sigmaSplay
          real(kind=realtype) :: dStartb
          real(kind=realtype), intent(in) :: epsE0, marchParameter, nuArea
          character(len=32), intent(in) :: bc1, bc2
          real(kind=realtype), intent(in) :: alphaP0, ratioGuess, cMax
          integer(kind=inttype), intent(in) :: numSmoothingPasses
          integer(kind=intType), intent(in) :: majorIndices(numLayers)
          logical, intent(in) :: extension_given
          integer(kind=intType), intent(in) :: numGuides
          integer(kind=intType), intent(in) :: guideIndices(numGuides)
          logical, intent(in) :: retainSpacing

          integer(kind=inttype), intent(out) :: fail
          real(kind=realtype), intent(out) :: rStartb(3*numNodes)

          real(kind=realtype) :: ratios(numLayers-1, numNodes-1), R(numLayers&
        &   , 3*numNodes), extension
          real(kind=realtype) :: Rb(numLayers, 3*numNodes)
          real(kind=realtype) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
          real(kind=realtype) :: r0b(3*numNodes), N0b(3, numNodes), S0b(&
        &   numNodes)
          real(kind=realtype) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*&
        &   numNodes)
          real(kind=realtype) :: rm1b(3*numNodes), Sm1b(numNodes), rSmoothedb(&
        &   3*numNodes)
          real(kind=realtype) :: rNext(3*numNodes), nNext(3, numNodes)
          real(kind=realtype) :: rNextb(3*numNodes)
          real(kind=realtype) :: rNext_in(3*numNodes)
          real(kind=realtype) :: dr(3*numNodes), d, dMax, dGrowth
          real(kind=realtype) :: db, dMaxb, dGrowthb, dVec(numLayers)
          real(kind=realtype) :: dPseudo, maxStretch, radius, eta
          real(kind=realtype) :: dPseudob, radiusb
          integer(kind=inttype) :: layerIndex, indexSubIter, cFactor
          integer(kind=inttype) :: minorIndex, cFactorVec(numLayers)
          intrinsic int
          integer :: arg1
          integer :: ad_to

          ! Placeholder before we get projection_b
          call py_projection(rStart, rNext, NNext, numNodes)

        ! find the characteristic radius of the mesh
          call findradius(rNext, numNodes, radius)

        ! find the desired marching distance
          extension = marchParameter
          dMax = radius*(extension-1.)

        ! compute the growth ratio necessary to match this distance
          call findratio(dMax, dStart, numLayers, ratioGuess, dGrowth)

          d = dStart
          cFactorVec = 0.
          dVec(1) = d

          ! Compute and store d and cFactor values for the backwards loop
          do layerIndex=1,numLayers-1
            call areaFactor(rNext, d, nuArea, numAreaPasses, bc1, bc2, guideIndices, retainSpacing,&
            numGuides, numNodes, S0, maxStretch)
            cFactor = ceiling(maxStretch/cMax)
            cFactorVec(layerIndex) = cFactor
            d = d*dGrowth
            dVec(layerIndex+1) = d
          end do

          ! Initialize the seeds as 0s
          db = 0.0_8
          rNextb = 0.0_8
          Sm1b = 0.0_8
          rm1b = 0.0_8
          dGrowthb = 0.0_8

          ! minorIndex is the index for the pseudomesh, which includes info
          ! from all minor steps, including subiterations
          minorIndex = nSubIters
          S0 = S(minorIndex-1, :)

          ! Loop over the major steps in reverse order
          do layerIndex=numLayers-1,1,-1

            ! call popreal8array(d, realtype/8)
            d = dVec(layerIndex)

            dGrowthb = dGrowthb + d*db
            db = dGrowth*db

            rNextb = rNextb + rb(layerIndex+1, :)
            rb(layerIndex+1, :) = 0.0_8
            dPseudo = d/cFactor
            dPseudob = 0.0_8
            r0b = 0.0_8

            ! call popinteger4(ad_to)
            ad_to = cFactorVec(layerIndex)

            ! Loop over a certain number of subiterations for each major step.
            ! This number, the cFactor, was determined previously
            do indexSubIter=ad_to,1,-1

              ! call popreal8array(N0, realtype*3*numNodes/8)
              N0 = N(minorIndex, :, :)
              ! call popreal8array(r0, realtype*3*numNodes/8)
              r0 = R_final(minorIndex-1, :)

              rNextb = rNextb + r0b
              r0b = 0.0_8
              ! call popreal8array(rm1, realtype*3*numNodes/8)
              if (minorIndex .le. 2) then
                rm1 = R_final(minorIndex-1, :)
                Sm1 = S(minorIndex-1, :)
              else
                rm1 = R_final(minorIndex-2, :)
                Sm1 = S(minorIndex-2, :)
              end if

              r0b = rm1b
              S0b = 0.0_8
              ! S0b = S0b + Sm1b
              ! call popreal8array(Sm1, realtype*numNodes/8)

              S0b = Sm1b
              rSmoothedb = 0.0_8

              ! call popreal8array(rNext, realtype*3*numNodes/8)
              rNext = R_initial_march(minorIndex, :)

              rSmoothedb = rNextb
              eta = layerIndex + 2

              call smoothing_main_b(rNext, rNextb, eta, alphaP0, &
        &                       numSmoothingPasses, numLayers, numNodes, &
        &                       rSmoothed, rSmoothedb)
              arg1 = layerIndex - 1
              !
              ! N0b = 0.
              ! rm1b = 0.
              ! Sm1b = 0.

              ! call popreal8array(rNext, realtype*3*numNodes/8)
              rNext = R_initial_march(minorIndex-1, :)


              call computematrices_main_b(r0, r0b, N0, N0b, S0, S0b, rm1, rm1b&
        &                             , Sm1, Sm1b, arg1, theta, sigmaSplay, bc1&
        &                             , bc2, numLayers, epsE0, guideIndices, retainSpacing,&
                                        rNext, rNextb, numNodes, numGuides)

              ! call popreal8array(S0, realtype*numNodes/8)
              S0 = S(minorIndex, :)


              call areaFactor_b(r0, r0b, dPseudo, dPseudob, nuArea, &
        &                   numAreaPasses, bc1, bc2, guideIndices, retainSpacing, &
                            numGuides, numNodes, S0, S0b, maxStretch)
              rNextb = 0.0_8
              minorIndex = minorIndex - 1
            end do
            db = db + dPseudob/cFactor

            ! call popinteger4array(cFactor, inttype/4)
            cFactor = cFactorVec(layerIndex)

            ! call popreal8array(S0, realtype*numNodes/8)
            if (minorIndex .le. 1) then
              S0 = S(minorIndex, :)
            else
              S0 = S(minorIndex-1, :)
            end if

            S0b = 0.0_8
            call areaFactor_b(r0, r0b, d, db, nuArea, numAreaPasses, bc1, bc2,&
                              guideIndices, retainSpacing, numGuides, numNodes, S0, S0b, maxStretch)
            rNextb = rNextb + r0b
          end do

          rNextb = rNextb + rb(1, :)
          ! rb(1, :) = 0.
          db = 0.0_8

          call areaFactor_b(rNext, rNextb, d, db, nuArea, numAreaPasses, bc1, &
        &               bc2, guideIndices, retainSpacing, numguides, numNodes, Sm1, Sm1b, maxStretch)

          rNextb = rNextb + rm1b

          call findRatio_b(dMax, dMaxb, dStart, dStartb, numLayers, ratioGuess&
        &              , dGrowth, dGrowthb)

          radiusb = (extension-1.)*dMaxb

          call findRadius_b(rNext, rNextb, numNodes, radius, radiusb)

          rStartb = 0.0_8
          rStartb = rNextb

        end subroutine march_b


        !=======================================================================

        subroutine releaseMemory()

          ! This subroutine just deallocates memory used by the marching code.
          ! Remember to call this function after you copy the outputs in Python.

          implicit none

          ! Deallocate output variables
          if (allocated(R_initial_march)) then
            deallocate(R_initial_march, R_smoothed, R_final, S, N)
          end if

        end subroutine releaseMemory

        !=================================================================

        subroutine computeMatrices(r0, N0, S0, rm1, Sm1, layerIndex, theta,&
        sigmaSplay, bc1, bc2, guideIndices, retainSpacing, numLayers, epsE0, rNext, numNodes, numGuides)

        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        integer(kind=intType), intent(in) :: numGuides
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing

        real(kind=realType), intent(out) :: rNext(3*numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), dVec(3), dVec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex-1, theta,&
        sigmaSplay, bc1, bc2, numLayers, epsE0, guideIndices, retainSpacing, rNext, numNodes, numGuides)

        end subroutine computeMatrices

        !=======================================================================

        subroutine computeMatrices_b(r0, r0_b, N0, N0_b, S0, S0_b, rm1, rm1_b, Sm1,&
        Sm1_b, layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, guideIndices, retainSpacing,&
        rNext, rNext_b, numNodes, numGuides)

        use hypsurfmain_b, only: computematrices_main_b
        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0, rNext_b(3*numNodes)
        character*32, intent(in) :: bc1, bc2
        integer(kind=intType), intent(in) :: numGuides
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing

        real(kind=realType), intent(out) :: rNext(3*numNodes)
        real(kind=realType), intent(out) :: r0_b(3*numNodes), N0_b(3, numNodes), S0_b(numNodes)
        real(kind=realType), intent(out) :: rm1_b(3*numNodes), Sm1_b(numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), dVec(3), dVec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computematrices_main_b(r0, r0_b, N0, N0_b, S0, S0_b, rm1, rm1_b, Sm1, Sm1_b,&
        layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, guideIndices, retainSpacing,&
        rNext, rNext_b, numNodes, numGuides)

        end subroutine computeMatrices_b

        !=================================================================

        subroutine computeMatrices_d(r0, r0_d, N0, N0_d, S0, S0_d, rm1, rm1_d, Sm1,&
        Sm1_d, layerIndex, theta, sigmaSplay, bc1, bc2, guideIndices, retainSpacing, numLayers,&
        epsE0, rNext, rNext_d, numNodes, numGuides)

        use hypsurfmain_d, only: computematrices_main_d
        implicit none

        integer(kind=intType), intent(in) :: layerIndex, numNodes, numLayers
        real(kind=realType), intent(in) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType), intent(in) :: rm1(3*numNodes), Sm1(numNodes), theta
        real(kind=realType), intent(in) :: sigmaSplay, epsE0
        character*32, intent(in) :: bc1, bc2
        integer(kind=intType), intent(in) :: numGuides
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing
        real(kind=realType), intent(in) :: r0_d(3*numNodes), N0_d(3, numNodes), S0_d(numNodes)
        real(kind=realType), intent(in) :: rm1_d(3*numNodes), Sm1_d(numNodes)

        real(kind=realType), intent(out) :: rNext(3*numNodes), rNext_d(3*numNodes)

        real(kind=realType) :: r_curr(3), r_next(3), r_prev(3), dVec(3), dVec_rot(3), eye(3, 3)
        real(kind=realType) :: K(3*numNodes, 3*numNodes)
        integer(kind=intType) :: index, i
        integer(kind=intType) :: ipiv(3*numNodes)
        integer(kind=intType) :: n, nrhs, ldK, ldf, info
        real(kind=realType) :: one, zero

        call computematrices_main_d(r0, r0_d, N0, N0_d, S0, S0_d, rm1, rm1_d, Sm1, Sm1_d,&
        layerIndex, theta, sigmaSplay, bc1, bc2, numLayers, epsE0, guideIndices, retainSpacing,&
        rNext, rNext_d, numNodes, numGuides)

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

        !=================================================================

        subroutine areafactor_test_d(r0, r0d, d, dd, nuarea, numareapasses, bc1, &
      &   bc2, guideindices, retainspacing, numguides, n, s, sd, maxstretch)

        use hypsurfmain_d, only: areafactor_d
        implicit none
        integer(kind=inttype), intent(in) :: n
        real(kind=realtype), intent(in) :: r0(3*n), d, nuarea
        real(kind=realtype), intent(in) :: r0d(3*n), dd
        integer(kind=inttype), intent(in) :: numareapasses
        character(len=32), intent(in) :: bc1, bc2
        real(kind=realtype), intent(out) :: s(n), maxstretch
        real(kind=realtype), intent(out) :: sd(n)
        integer(kind=inttype), intent(in) :: numguides
        integer(kind=inttype), intent(in) :: guideindices(numguides)
        logical, intent(in) :: retainspacing

        call areafactor_d(r0, r0d, d, dd, nuarea, numareapasses, bc1, &
      &   bc2, guideindices, retainspacing, numguides, n, s, sd, maxstretch)

        end subroutine

        !=================================================================

        subroutine areafactor_test_b(r0, r0b, d, db, nuarea, numareapasses, bc1, &
      &   bc2, guideindices, retainspacing, numguides, n, s, sb, maxstretch)

        use hypsurfmain_b, only: areafactor_b
        implicit none
        integer(kind=inttype), intent(in) :: n
        real(kind=realtype), intent(in) :: r0(3*n), d, nuarea
        real(kind=realtype), intent(out) :: r0b(3*n), db
        integer(kind=inttype), intent(in) :: numareapasses
        character(len=32), intent(in) :: bc1, bc2
        real(kind=realtype), intent(in) :: s(n), maxstretch
        real(kind=realtype), intent(in) :: sb(n)
        integer(kind=inttype), intent(in) :: numguides
        integer(kind=inttype), intent(in) :: guideindices(numguides)
        logical, intent(in) :: retainspacing

        call areafactor_b(r0, r0b, d, db, nuarea, numareapasses, bc1, &
      &   bc2, guideindices, retainspacing, numguides, n, s, sb, maxstretch)

        end subroutine


      end module
