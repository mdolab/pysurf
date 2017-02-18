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

      ! R_initial_march(i) has the nodes of layer i right after the solution of the linear system of the previous layer
      real(kind=realType), dimension(:,:), allocatable :: R_initial_march
      ! R_smoothed(i) has the nodes of layer i after the smoothing step
      real(kind=realType), dimension(:,:), allocatable :: R_smoothed
      ! R_final(i) has the nodes of layer i after the first projection step
      real(kind=realType), dimension(:,:), allocatable :: R_projected
      ! R_remeshed(i) has the nodes of layer i after the remeshing step
      real(kind=realType), dimension(:,:), allocatable :: R_remeshed
      ! R_final(i) has the nodes of layer i after the final projection step
      real(kind=realType), dimension(:,:), allocatable :: R_final

      real(kind=realType), dimension(:,:), allocatable :: Sm1_hist
      real(kind=realType), dimension(:,:), allocatable :: S0_hist
      real(kind=realType), dimension(:,:, :), allocatable :: N_projected
      real(kind=realType), dimension(:,:, :), allocatable :: N_final


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

        real(kind=realType) :: normArcLength(numNodes)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)
        real(kind=realType) :: rNext_in(3*numNodes), rRemeshed(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta, min_ratio, ratios_small(3, numNodes-1)
        real(kind=realType) :: extension, growthRatio
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor, nSubIters
        integer(kind=intType) :: arraySize, nAllocations

        integer(kind=intType) :: breakNodes(numGuides+2), intervalID, node1, node2
        integer(kind=intType) :: numInInterval, numIntervals

        real(kind=realType), dimension(:,:), allocatable :: ext_temp_real
        real(kind=realType), dimension(:,:), allocatable :: ext_temp_real_1dim
        real(kind=realType), dimension(:,:, :), allocatable :: ext_temp_real_N
        real(kind=realType), dimension(:,:), allocatable :: ext_R_initial_march
        real(kind=realType), dimension(:,:), allocatable :: ext_R_smoothed
        real(kind=realType), dimension(:,:), allocatable :: ext_R_projected
        real(kind=realType), dimension(:,:), allocatable :: ext_R_remeshed
        real(kind=realType), dimension(:,:), allocatable :: ext_R_final
        real(kind=realType), dimension(:,:), allocatable :: ext_S0
        real(kind=realType), dimension(:,:), allocatable :: ext_Sm1
        real(kind=realType), dimension(:,:, :), allocatable :: ext_N_projected
        real(kind=realType), dimension(:,:, :), allocatable :: ext_N_final

        ! Project onto the surface or curve (if applicable)
        call py_projection(rStart, rNext, NNext, numNodes)

        ! Initialize step size and total marched distance
        d = dStart
        dTot = 0

        if (extension_given) then

           ! Find the characteristic radius of the mesh
           call findRadius(rNext, numNodes, radius)

           print *,'radius'
           print *,radius

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

        !===========================================================

        ! Compute arclengths of each subinterval if the user requested the remesh feature

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

        ! Initialize arclenghts
        normArcLength = 0.0

        ! normArcLength is the normalized arclengths of each subinterval. A new subinterval
        ! starts when normArcLength(i) = 0. For instance, If we have two subintervals we
        ! might have:
        ! normArcLength = [0.0, 0.3, 0.7, 0.0, 0.1, 0.4, 0.8, 1.0]

        ! Compute normalized arclengths of each subinterval
        if (retainSpacing) then

           ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
           do intervalID=1, numGuides+1

              ! Get indices of the first and last nodes of the interval
              node1 = breakNodes(intervalID)
              node2 = breakNodes(intervalID+1)

              numInInterval = node2 - node1 + 1

              ! Compute the normalized arc-lengths of this interval
              call compute_arc_length(rStart(3*(node1-1)+1:3*node2), numInInterval, normArcLength(node1:node2))

           end do

        end if

        !===========================================================

        fail = 0
        nSubIters = 0

        ! Initialize extended connectivity arrays for the intersection curves.
        ! These arrays will be cropped to get the actual outputs.
        ! For now we will allocate arrays of size arraySize. We will increase them
        ! if necessary
        arraySize = numLayers
        allocate(ext_R_initial_march(arraySize, 3*numNodes))
        allocate(ext_R_smoothed(arraySize, 3*numNodes))
        allocate(ext_R_projected(arraySize, 3*numNodes))
        allocate(ext_R_remeshed(arraySize, 3*numNodes))
        allocate(ext_R_final(arraySize, 3*numNodes))
        allocate(ext_S0(arraySize, numNodes))
        allocate(ext_Sm1(arraySize, numNodes))
        allocate(ext_N_projected(arraySize, 3, numNodes))
        allocate(ext_N_final(arraySize, 3, numNodes))

        ! Initialize values
        ext_R_initial_march = 0.
        ext_R_smoothed = 0.
        ext_R_projected = 0.
        ext_R_remeshed = 0.
        ext_R_final = 0.
        ext_S0 = 0.
        ext_Sm1 = 0.
        ext_N_projected = 0.
        ext_N_final = 0.

        ! Store the initial curves
        ext_R_initial_march(1, :) = rNext
        ext_R_smoothed(1, :) = rNext
        ext_R_projected(1, :) = rNext
        ext_R_remeshed(1, :) = rNext
        ext_R_final(1, :) = rNext
        ext_N_projected(1, :, :) = NNext
        ext_N_final(1, :, :) = NNext
        ext_S0 = 0.
        ext_Sm1 = 0.

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
          N0 = NNext

          ! Compute the new area factor for the desired marching distance
          call areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, guideIndices, &
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

            ! Recompute areas with the pseudo-step to get S0 and Sm1
            call areaFactor(rm1, dPseudo, nuArea, numAreaPasses, bc1, bc2, guideIndices, &
            numGuides, numNodes, Sm1, maxStretch)
            call areaFactor(r0, dPseudo, nuArea, numAreaPasses, bc1, bc2, guideIndices, &
            numGuides, numNodes, S0, maxStretch)

            ! Increase the number of subiterations done so far
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

              ! REALLOCATING ext_R_smoothed

              ! Create temporary array and initialize
              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_smoothed

              ! Now move the new allocation back to ext_R_smoothed.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_smoothed)

              ! REALLOCATING ext_R_projected

              ! Create temporary array and initialize
              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_projected

              ! Now move the new allocation back to ext_R_projected.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_projected)

              ! REALLOCATING ext_R_remeshed

              ! Create temporary array and initialize
              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_remeshed

              ! Now move the new allocation back to ext_R_remeshed.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_remeshed)

              ! REALLOCATING ext_R_final

              ! Create temporary array and initialize
              allocate(ext_temp_real(nAllocations*arraySize, 3*numNodes))
              ext_temp_real = 0.0

              ! Tranfer data from original array
              ext_temp_real(:(nAllocations-1)*arraySize, :) = ext_R_final

              ! Now move the new allocation back to ext_R_final.
              ! ext_temp_real is deallocated in this process.
              call move_alloc(ext_temp_real, ext_R_final)

              ! REALLOCATING ext_Sm1

              ! Create temporary array and initialize
              allocate(ext_temp_real_1dim(nAllocations*arraySize, numNodes))
              ext_temp_real_1dim = 0.0

              ! Tranfer data from original array
              ext_temp_real_1dim(:(nAllocations-1)*arraySize, :) = ext_Sm1

              ! Now move the new allocation back to ext_Sm1.
              ! ext_temp_real_1dim is deallocated in this process.
              call move_alloc(ext_temp_real_1dim, ext_Sm1)

              ! REALLOCATING ext_S0

              ! Create temporary array and initialize
              allocate(ext_temp_real_1dim(nAllocations*arraySize, numNodes))
              ext_temp_real_1dim = 0.0

              ! Tranfer data from original array
              ext_temp_real_1dim(:(nAllocations-1)*arraySize, :) = ext_S0

              ! Now move the new allocation back to ext_S0.
              ! ext_temp_real_1dim is deallocated in this process.
              call move_alloc(ext_temp_real_1dim, ext_S0)

              ! REALLOCATING ext_N_projected

              ! Create temporary array and initialize
              allocate(ext_temp_real_N(nAllocations*arraySize, 3, numNodes))
              ext_temp_real_N = 0.0

              ! Tranfer data from original array
              ext_temp_real_N(:(nAllocations-1)*arraySize, :, :) = ext_N_projected

              ! Now move the new allocation back to ext_N.
              ! ext_temp_real_N is deallocated in this process.
              call move_alloc(ext_temp_real_N, ext_N_projected)

              ! REALLOCATING ext_N_final

              ! Create temporary array and initialize
              allocate(ext_temp_real_N(nAllocations*arraySize, 3, numNodes))
              ext_temp_real_N = 0.0

              ! Tranfer data from original array
              ext_temp_real_N(:(nAllocations-1)*arraySize, :, :) = ext_N_final

              ! Now move the new allocation back to ext_N.
              ! ext_temp_real_N is deallocated in this process.
              call move_alloc(ext_temp_real_N, ext_N_final)

            end if

            ! Store spacings for the reverse run
            ext_Sm1(nSubIters,:) = Sm1
            ext_S0(nSubIters,:) = S0

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

            ! Save this projected row for later use in the derivatives
            ext_R_projected(nSubIters+1, :) = rNext

            ! Save N0 into N for backwards derivative computation
            ext_N_projected(nSubIters+1, :, :) = NNext

            if (retainSpacing) then

              ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
              do intervalID=1, numGuides+1

                ! Get indices of the first and last nodes of the interval
                node1 = breakNodes(intervalID)
                node2 = breakNodes(intervalID+1)

                numInInterval = node2 - node1 + 1

                ! Temporarily set the first and last arclengths of the normalized interval
                normArcLength(node1) = 0.0
                normArcLength(node2) = 1.0

                ! Redistribute the nodes
                call redistribute_nodes_by_arc_length(rNext(3*(node1-1)+1:3*node2), &
                                                      normArcLength(node1:node2), &
                                                      numInInterval, &
                                                      rRemeshed(3*(node1-1)+1:3*node2))

                ! Save this projected row for later use in the derivatives
                ext_R_remeshed(nSubIters+1, :) = rRemeshed

              end do

              ! Project the remeshed curve back onto the surface
              call py_projection(rRemeshed, rNext, NNext, numNodes)

            end if

            ! Save this final projected row for later use in the backwards
            ! derivative computation if it's the last one of the subiterations
            ext_R_final(nSubIters+1, :) = rNext

            ! Save N0 into N for backwards derivative computation
            ext_N_final(nSubIters+1, :, :) = NNext

            ! Save Sm1 for backwards derivative computation
            !if (nSubiters .eq. 1) then ! Also store Sm1
            !   ext_Sm1(nSubIters, :) = Sm1
            !end if
            !ext_Sm1(nSubIters+1, :) = S0

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
        allocate(R_projected(nSubIters+1, 3*numNodes))
        R_projected(:,:) = ext_R_projected(:nSubIters+1, :)

        ! Crop extended array
        allocate(R_remeshed(nSubIters+1, 3*numNodes))
        R_remeshed(:,:) = ext_R_remeshed(:nSubIters+1, :)

        ! Crop extended array
        allocate(R_final(nSubIters+1, 3*numNodes))
        R_final(:,:) = ext_R_final(:nSubIters+1, :)

        ! Crop extended array
        allocate(Sm1_hist(nSubIters+1, numNodes))
        Sm1_hist(:,:) = ext_Sm1(:nSubIters+1, :)

        ! Crop extended array
        allocate(S0_hist(nSubIters+1, numNodes))
        S0_hist(:,:) = ext_S0(:nSubIters+1, :)

        ! Crop extended array
        allocate(N_projected(nSubIters+1, 3, numNodes))
        N_projected(:,:, :) = ext_N_projected(:nSubIters+1, :, :)

        ! Crop extended array
        allocate(N_final(nSubIters+1, 3, numNodes))
        N_final(:,:, :) = ext_N_final(:nSubIters+1, :, :)

        call qualityCheck(R, 0, numLayers, numNodes, fail, ratios)

        print *, ''
        print *, 'Final mesh quality minimum:', minval(ratios)
        print *, ''


        !do layerIndex=1,size(Sm1_hist,1)
        !   print *,Sm1_hist(layerIndex,:)
        !end do

        end subroutine march

        !=======================================================================

        subroutine march_d(py_projection_d, rStart, rStartd, R_projected_in, R_final_in, &
        N_projected_in, N_final_in, dStart, theta, sigmaSplay, bc1, bc2, &
        epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax, guideIndices, retainSpacing,&
        extension_given, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, numGuides, nSubIters_in, R, Rd, fail, ratios)

        use hypsurfMain_d, only: computeMatrices_main_d, smoothing_main_d, areafactor_d, findRadius_d, findRatio_d, &
                                 compute_arc_length_d, redistribute_nodes_by_arc_length_d
        implicit none

        external py_projection_d
        real(kind=realType), intent(in) :: rStart(3*numNodes), rStartd(3*numNodes)
        real(kind=realType), dimension(nSubIters_in, 3*numNodes), intent(in) :: R_projected_in, R_final_in
        real(kind=realType), dimension(nSubIters_in, 3, numNodes), intent(in) :: N_projected_in, N_final_in
        real(kind=realType), intent(in) :: dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, marchParameter, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses
        logical, intent(in) :: extension_given
        integer(kind=intType), intent(in) :: numGuides, nSubIters_in
        integer(kind=intType), intent(in) :: guideIndices(numGuides)
        logical, intent(in) :: retainSpacing

        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1)
        real(kind=realType), intent(out) :: R(numLayers, 3*numNodes), Rd(numLayers, 3*numNodes)

        real(kind=realType) :: normArcLength(numNodes), normArcLengthd(numNodes)

        real(kind=realType):: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
        real(kind=realType) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes), rRemeshed(3*numNodes)
        real(kind=realType) :: rNext_in(3*numNodes)
        real(kind=realType) :: dr(3*numNodes), d, dTot, dMax, dGrowth
        real(kind=realType) :: dPseudo, maxStretch, radius, eta, extension
        integer(kind=intType) :: layerIndex, indexSubIter, cFactor, nSubIters, nProjs
        integer(kind=intType) :: arraySize, nAllocations, layerID

        integer(kind=intType) :: breakNodes(numGuides+2), intervalID, node1, node2
        integer(kind=intType) :: numInInterval, numIntervals

        ! Derivative variables
        real(kind=realType):: r0d(3*numNodes), N0d(3, numNodes), S0d(numNodes)
        real(kind=realType) :: rm1d(3*numNodes), Sm1d(numNodes), rNextd(3*numNodes)
        real(kind=realType) :: rSmoothedd(3*numNodes), NNextd(3, numNodes), rRemeshedd(3*numNodes)
        real(kind=realType) :: r0d_copy(3*numNodes), dd, dGrowthd, dMaxd
        real(kind=realType) :: dStartd, radiusd, dPseudod


        ! Need py_projection_d here (Remember that rNext and NNext are inputs here)
        rNext = R_projected_in(1,:)
        NNext = N_projected_in(1,:,:)
        layerID = 0
        rNextd = 0.
        NNextd = 0.
        call py_projection_d(rStart, rStartd, rNext, rNextd, NNext, NNextd, layerID, numNodes)
        !call py_projection_d(rStart, rStartd, R_initial_march(1,:), rNextd, N(1,:,:), NNextd, 0, numNodes)
        !rNextd = rStartd
        !NNextd = 0.0

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

        !===========================================================

        ! Really only do this if we're remeshing
        breakNodes(1) = 1
        breakNodes(2:numGuides+1) = guideIndices
        breakNodes(numGuides+2) = numNodes

        !!! Still need to do this!
        ! Remove duplicate entries from breakNodes if they exist.
        ! This will only occur at the beginning or end because we've already
        ! sorted the internal part of breakNodes by sorting guideIndices.
        ! if (breakNodes(1) .eq. breakNodes(2))

        numIntervals = numGuides + 1

        ! Initialize arclenghts
        normArcLength = 0.0

        ! normArcLength is the normalized arclengths of each subinterval. A new subinterval
        ! starts when normArcLength(i) = 0. For instance, If we have two subintervals we
        ! might have:
        ! normArcLength = [0.0, 0.3, 0.7, 0.0, 0.1, 0.4, 0.8, 1.0]

        ! Compute normalized arclengths of each subinterval
        if (retainSpacing) then

           ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
           do intervalID=1, numGuides+1

              ! Get indices of the first and last nodes of the interval
              node1 = breakNodes(intervalID)
              node2 = breakNodes(intervalID+1)

              numInInterval = node2 - node1 + 1

              ! Compute the normalized arc-lengths of this interval
              call compute_arc_length_d(rStart(3*(node1-1)+1:3*node2), rStartd(3*(node1-1)+1:3*node2), &
                                        numInInterval, &
                                        normArcLength(node1:node2), normArcLengthd(node1:node2))

           end do

        end if

        !===========================================================

        fail = 0
        nSubIters = 0

        Rd = 0.
        Rd(1, :) = rNextd
        R(1, :) = rNext

        dd = 0.

        ! Sub iteration counter
        nSubIters = 0

        ! Projection counter
        nProjs = 0

        do layerIndex=1,numLayers-1
          ! Get the coordinates computed by the previous iteration
          r0d = rNextd
          r0 = rNext
          N0d = NNextd
          N0 = NNext

          S0d = 0.
          ! Compute the new area factor for the desired marching distance
          call areaFactor_d(r0, r0d, d, dd, nuArea, numAreaPasses, bc1, bc2,&
          guideIndices, numguides, numNodes, S0, S0d, maxStretch)

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
             Sm1d = 0.
             ! Recompute areas with the pseudo-step
             call areaFactor_d(rm1, rm1d, dPseudo, dPseudod, nuArea, numAreaPasses,&
             bc1, bc2, guideIndices, numguides, numNodes, Sm1, Sm1d, maxStretch)
             call areaFactor_d(r0, r0d, dPseudo, dPseudod, nuArea, numAreaPasses,&
             bc1, bc2, guideIndices, numguides, numNodes, S0, S0d, maxStretch)

             ! Increment subiteration counter
             nSubIters = nSubIters + 1

             ! March using the pseudo-marching distance
             eta = layerIndex+2

             rNextd = 0.
             call computeMatrices_main_d(r0, r0d, N0, N0d, S0, S0d, rm1, rm1d, &
             Sm1, Sm1d, layerIndex-1, theta, sigmaSplay, bc1, bc2, &
             numLayers, epsE0, guideIndices, retainSpacing, rNext, rNextd, numNodes, numGuides)

             rSmoothedd = 0.
             call smoothing_main_d(rNext, rNextd, eta, alphaP0, numSmoothingPasses, numLayers, numNodes, rSmoothed, rSmoothedd)

             ! Get the values of the projection, as py_projection_d uses rNext and NNext as inputs
             rNext = R_projected_in(nSubIters+1,:)
             NNext = N_projected_in(nSubIters+1,:,:)
             rNextd = 0.0
             NNextd = 0.0
             nProjs = nProjs+1

             call py_projection_d(rSmoothed, rSmoothedd, &
                  rNext, rNextd, &
                  NNext, NNextd, &
                  nProjs, numNodes)

             ! Check if we need to remesh the current layer to keep spacing
             if (retainSpacing) then

                ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
                do intervalID=1, numGuides+1

                   ! Get indices of the first and last nodes of the interval
                   node1 = breakNodes(intervalID)
                   node2 = breakNodes(intervalID+1)

                   numInInterval = node2 - node1 + 1

                   ! Temporarily set the first and last arclengths of the normalized interval
                   normArcLength(node1) = 0.0
                   normArcLength(node2) = 1.0

                   ! Redistribute the nodes
                   call redistribute_nodes_by_arc_length_d(rNext(3*(node1-1)+1:3*node2), rNextd(3*(node1-1)+1:3*node2), &
                                                           normArcLength(node1:node2), normArcLengthd(node1:node2), &
                                                           numInInterval, &
                                                           rRemeshed(3*(node1-1)+1:3*node2), rRemeshedd(3*(node1-1)+1:3*node2))

                end do

                ! Project the remeshed curve back onto the surface
                rNext = R_final_in(nSubIters+1,:)
                NNext = N_final_in(nSubIters+1,:,:)
                rNextd = 0.0
                NNextd = 0.0
                nProjs = nProjs + 1 ! This indicates which projection history we should use
                call py_projection_d(rRemeshed, rRemeshedd, rNext, rNextd, NNext, NNextd, nProjs, numNodes)

             end if

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
             N0d = NNextd
             N0 = NNext

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
        !=======================================================================
        !=======================================================================

        subroutine march_b(py_projection_b, rStart, rStartb, R_initial_march_in, R_smoothed_in, &
        R_projected_in, R_remeshed_in, R_final_in, N_projected_in, N_final_in, Sm1_hist_in, S0_hist_in, dStart, &
        theta, sigmaSplay, bc1, bc2, epsE0, alphaP0, marchParameter, nuArea, ratioGuess, cMax,&
        guideIndices, retainSpacing, extension_given, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, nSubIters, numGuides, numProjs, R, Rb)

          ! numProjs: integer -> This is the total number of projections done during
          !                      the forward step.

          use hypsurfMain_b, only: computeMatrices_main_b, smoothing_main_b, areaFactor_b, &
                                   findRadius_b, findRatio_b, compute_arc_length_b, &
                                   redistribute_nodes_by_arc_length_b
          implicit none

          external py_projection_b
          integer(kind=inttype), intent(in) :: numNodes, numLayers, numAreaPasses, nSubIters
          real(kind=realType), dimension(nSubIters, 3*numNodes), intent(in) :: R_initial_march_in, R_smoothed_in
          real(kind=realType), dimension(nSubIters, 3*numNodes), intent(in) :: R_projected_in, R_final_in
          real(kind=realType), dimension(nSubIters, 3*numNodes), intent(in) :: R_remeshed_in
          real(kind=realType), intent(in) :: N_projected_in(nSubIters, 3, numNodes)
          real(kind=realType), intent(in) :: N_final_in(nSubIters, 3, numNodes)
          real(kind=realType), intent(in) :: Sm1_hist_in(nSubIters, numNodes), S0_hist_in(nSubIters, numNodes)
          real(kind=realtype), intent(in) :: rStart(3*numNodes), dStart, theta, sigmaSplay
          real(kind=realtype) :: dStartb, dStartb_dummy
          real(kind=realtype), intent(in) :: epsE0, marchParameter, nuArea
          character(len=32), intent(in) :: bc1, bc2
          real(kind=realtype), intent(in) :: alphaP0, ratioGuess, cMax
          integer(kind=inttype), intent(in) :: numSmoothingPasses
          logical, intent(in) :: extension_given
          integer(kind=intType), intent(in) :: numGuides, numProjs
          integer(kind=intType), intent(in) :: guideIndices(numGuides)
          logical, intent(in) :: retainSpacing
          real(kind=realType), intent(in) :: R(numLayers, 3*numNodes), Rb(numLayers, 3*numNodes)

          real(kind=realtype), intent(out) :: rStartb(3*numNodes)
          real(kind=realtype) :: rStart_dummy(3*numNodes), rStartb_dummy(3*numNodes)

          integer(kind=intType) :: breakNodes(numGuides+2), intervalID, node1, node2
          integer(kind=intType) :: numInInterval, numIntervals

          real(kind=realtype) :: extension
          real(kind=realtype) :: r0(3*numNodes), N0(3, numNodes), S0(numNodes)
          real(kind=realtype) :: r0b(3*numNodes), N0b(3, numNodes), S0b(numNodes)
          real(kind=realtype) :: r0b_dummy(3*numNodes), N0b_dummy(3, numNodes), S0b_dummy(numNodes)
          real(kind=realtype) :: rm1(3*numNodes), Sm1(numNodes), rSmoothed(3*numNodes)
          real(kind=realtype) :: rm1b(3*numNodes), Sm1b(numNodes), rSmoothedb(3*numNodes)
          real(kind=realtype) :: rm1b_dummy(3*numNodes), Sm1b_dummy(numNodes), rSmoothedb_dummy(3*numNodes)
          real(kind=realtype) :: rNext(3*numNodes), NNext(3, numNodes)
          real(kind=realType) :: NNextb(3, numNodes), NNextb_dummy(3, numNodes)
          real(kind=realtype) :: rNextb(3*numNodes), rRemeshed(3*numNodes), rRemeshedb(3*numNodes)
          real(kind=realtype) :: rNextb_dummy(3*numNodes), rRemeshedb_dummy(3*numNodes)
          real(kind=realtype) :: rNext1(3*numNodes), rNext2(3*numNodes)
          real(kind=realType) :: rNext1b(3*numNodes), rNext2b(3*numNodes)
          real(kind=realType) :: rNext1b_dummy(3*numNodes), rNext2b_dummy(3*numNodes)
          real(kind=realType) :: Nnext2(3,numNodes), Nnext2b(3,numNodes), Nnext2b_dummy(3,numNodes)
          real(kind=realType) :: normArcLength(numNodes), normArcLengthb(numNodes), normArcLength_dummy(numNodes)
          real(kind=realType) :: normArcLengthb_dummy(numNodes)
          real(kind=realtype) :: dr(3*numNodes), d, dMax, dGrowth
          real(kind=realtype) :: db, dMaxb, dGrowthb, dVec(numLayers)
          real(kind=realtype) :: db_dummy, dMaxb_dummy, dGrowthb_dummy
          real(kind=realtype) :: dPseudo, maxStretch, maxStretchVec(numLayers), radius, eta
          real(kind=realtype) :: dPseudob, radiusb
          real(kind=realtype) :: dPseudob_dummy, radiusb_dummy
          integer(kind=inttype) :: layerIndex, indexSubIter, cFactor, layerID, projID
          integer(kind=inttype) :: minorIndex, cFactorVec(numLayers), m1Index
          intrinsic int
          integer :: arg10, ii
          integer :: ad_to

          ! ====================================

          ! COMPUTE GROWTH RATIO AND RADIUS

          if (extension_given) then

             ! Find the characteristic radius of the mesh
             rNext = R_final_in(1,:)
             call findRadius(rNext, numNodes, radius)

             ! Find the desired marching distance
             ! Here marchParameter = extension
             dMax = radius * (marchParameter-1.)

             ! Compute the growth ratio necessary to match this distance
             call findRatio(dMax, dStart, numLayers, ratioGuess, dGrowth)
          else
             ! Here marchParameter = growthRatio
             dGrowth = marchParameter
          end if

          ! ====================================

          ! CREATE ARRAY THAT DEFINES REMESH INTERVAL

          ! Compute normalized arclengths of each subinterval
          if (retainSpacing) then

             ! Really only do this if we're remeshing
             breakNodes(1) = 1
             breakNodes(2:numGuides+1) = guideIndices
             breakNodes(numGuides+2) = numNodes

!!! Still need to do this!
             ! Remove duplicate entries from breakNodes if they exist.
             ! This will only occur at the beginning or end because we've already
             ! sorted the internal part of breakNodes by sorting guideIndices.
             ! if (breakNodes(1) .eq. breakNodes(2))

             numIntervals = numGuides + 1

             ! Initialize arclenghts
             normArcLength = 0.0

             ! normArcLength is the normalized arclengths of each subinterval. A new subinterval
             ! starts when normArcLength(i) = 0. For instance, If we have two subintervals we
             ! might have:
             ! normArcLength = [0.0, 0.3, 0.7, 0.0, 0.1, 0.4, 0.8, 1.0]

             ! Do the remesh for each subinterval. We have a self.arcLength entry for each subinterval
             do intervalID=1, numGuides+1

                ! Get indices of the first and last nodes of the interval
                node1 = breakNodes(intervalID)
                node2 = breakNodes(intervalID+1)

                print *,'computing interval'
                print *,node1
                print *,node2

                numInInterval = node2 - node1 + 1

                ! Compute the normalized arc-lengths of this interval
                call compute_arc_length(rStart(3*(node1-1)+1:3*node2), &
                     &                  numInInterval, &
                     &                  normArcLength(node1:node2))

             end do

          end if

          ! ====================================

          ! CREATE ARRAYS OF d AND cFactor

          d = dStart
          cFactorVec = 0.
          maxStretchVec = 0.

          minorIndex = 1
          ! Compute and store d and cFactor values for the backwards loop
          do layerIndex=1,numLayers-1

             ! Assign current d
             dVec(layerIndex) = d

             ! Get nodal coordinates of the current layer
             r0 = R_final_in(minorIndex, :)

             ! Compute the cFactor
             call areaFactor(r0, d, nuArea, numAreaPasses, bc1, bc2, guideIndices, &
                  numGuides, numNodes, S0, maxStretch)

             cFactor = ceiling(maxStretch/cMax)

             ! Store cFactor
             cFactorVec(layerIndex) = cFactor

             ! Store max stretch
             maxStretchVec(layerIndex) = maxStretch

             d = d*dGrowth
             minorIndex = minorIndex + cFactor

          end do

          ! ====================================

          ! Initialize number of projections
          projID = numProjs-1
          minorIndex = minorIndex + 1

          ! ====================================

          ! Initialize derivatives
          db = 0.0
          rnextb = 0.0
          nnextb = 0.0
          normarclengthb = 0.0
          rm1b = 0.0
          sm1b = 0.0
          rremeshedb = 0.0
          dgrowthb = 0.0

          DO layerindex=numlayers-1,1,-1
             !CALL POPREAL4ARRAY(d, realtype/4)
             d = dVec(layerIndex)

             ! Update derivatives
             dgrowthb = dgrowthb + d*db
             db = dgrowth*db
             rnextb = rnextb + rb(layerindex+1, :)
             !rb(layerindex+1, :) = 0.0
             dpseudob = 0.0
             r0b = 0.0
             n0b = 0.0

             ! Retrieve cFactor and maxStretch of the current layer
             cFactor = cFactorVec(layerIndex)
             maxStretch = maxStretchVec(layerIndex)

             ! Compute the pseudo-marching distance
             dPseudo = d/cFactor

             ! Set do loop bound
             !CALL POPINTEGER4(ad_to)
             ad_to = cFactor

             DO indexsubiter=ad_to,1,-1

                ! Update index of the current subiteration
                minorIndex = minorIndex - 1

                ! Assign index to the m1 layer
                m1Index = max(1,minorIndex-2)

                nnextb = nnextb + n0b
                rnextb = rnextb + r0b
                r0b = 0.0

                r0b = rm1b

                s0b = sm1b

                n0b = 0.0 ! I Think this should be restarted on every loop
                sm1b = 0.0 ! I Think this should be restarted on every loop
                rm1b = 0.0 ! I Think this should be restarted on every loop

                !CALL POPCONTROL1B(branch)
                !IF (branch .EQ. 0) THEN

                IF (retainSpacing) THEN

                   !CALL POPREAL4ARRAY(rremeshed, realtype*3*numnodes/4)
                   rRemeshed = R_remeshed_in(minorIndex, :)

                   !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
                   rNext = R_final_in(minorIndex,:)

                   !CALL POPREAL4ARRAY(nnext, realtype*3*numnodes/4)
                   NNext = N_final_in(minorIndex, :,:)

                   rRemeshedb = 0.0

                   CALL PY_PROJECTION_B(rremeshed, rremeshedb, rnext, rnextb, &
                        &               nnext, nnextb, projID, numnodes)
                   projID = projID - 1

                   rnext2b = 0.0

                   ! Retrieve data
                   rNext2 = R_projected_in(minorIndex,:)
                   rRemeshed = R_remeshed_in(minorIndex,:)
                   normArcLength_dummy = normArcLength

                   DO intervalid=numguides+1,1,-1
                      node1 = breaknodes(intervalid)
                      node2 = breaknodes(intervalid+1)
                      numininterval = node2 - node1 + 1
                      normarclengthb_dummy = 0.0
                      CALL REDISTRIBUTE_NODES_BY_ARC_LENGTH_B(rnext2(3*(node1-1)+1:3*node2), &
                           &                                  rnext2b(3*(node1-1)+1:3*node2), &
                           &                                  normarclength_dummy(node1:node2), &
                           &                                  normarclengthb_dummy(node1:node2), &
                           &                                  numininterval, &
                           &                                  rremeshed(3*(node1-1)+1:3*node2), &
                           &                                  rremeshedb(3*(node1-1)+1:3*node2))

                      normarclengthb = normarclengthb + normarclengthb_dummy
                      !CALL POPREAL4ARRAY(normarclength(node2), realtype/4)
                      normarclengthb(node2) = 0.0
                      !CALL POPREAL4ARRAY(normarclength(node1), realtype/4)
                      normarclengthb(node1) = 0.0

                   END DO

                   ! The normals of the previous projection are discarded due to the
                   ! redistribution step. Therefore, we can set their seeds to zero
                   nnext2b = 0.0
                ELSE

                   nnext2b = nnextb

                   !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
                   rNext = R_final_in(minorIndex,:)

                   rnext2b = rnextb
                END IF

                !CALL POPREAL4ARRAY(rsmoothed, realtype*3*numnodes/4)
                rSmoothed = R_smoothed_in(minorIndex,:)

                !CALL POPREAL4ARRAY(rnext2, realtype*3*numnodes/4)
                rNext2 = R_projected_in(minorIndex,:)

                !CALL POPREAL4ARRAY(nnext, realtype*3*numnodes/4)
                NNext2 = N_projected_in(minorIndex,:,:)

                rsmoothedb_dummy = 0.0
                CALL PY_PROJECTION_B(rsmoothed, rsmoothedb_dummy, rnext2, rnext2b, &
                     &               nnext2, nnext2b, projID, numnodes)
                rsmoothedb = rsmoothedb_dummy !I believe this should not be accumulated
                projID = projID - 1

                eta = layerindex + 2

                !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
                rNext1 = R_initial_march_in(minorIndex,:)

                !CALL POPREAL4ARRAY(rsmoothed, realtype*3*numnodes/4)
                rSmoothed = R_smoothed_in(minorIndex,:)

                rNext1b_dummy = 0.0
                CALL SMOOTHING_MAIN_B(rnext1, rnext1b_dummy, eta, alphap0, &
                     &                numsmoothingpasses, numlayers, numnodes, &
                     &                rsmoothed, rsmoothedb)
                rNext1b = rNext1b_dummy ! I believe this should not be accumulated

                arg10 = layerindex - 1

                !CALL POPREAL4ARRAY(rnext1, realtype*3*numnodes/4)
                rNext1 = R_initial_march_in(minorIndex,:)

                !CALL POPREAL4ARRAY(n0, realtype*3*numnodes/4)
                N0 = N_final_in(minorIndex-1,:,:) ! These three were at the beggining of the do loop
                !CALL POPREAL4ARRAY(r0, realtype*3*numnodes/4)
                r0 = R_final_in(minorIndex-1,:)
                !CALL POPREAL4ARRAY(rm1, realtype*3*numnodes/4)
                rm1 = R_final_in(m1Index,:)

                !CALL POPREAL4ARRAY(s0, realtype*numnodes/4)
                S0 = S0_hist_in(minorIndex-1,:)
                Sm1 = Sm1_hist_in(minorIndex-1,:)

                r0b_dummy = 0.0
                n0b_dummy = 0.0
                s0b_dummy = 0.0
                rm1b_dummy = 0.0
                sm1b_dummy = 0.0
                CALL COMPUTEMATRICES_MAIN_B(r0, r0b_dummy, n0, n0b_dummy, &
                     &                      s0, s0b_dummy, rm1, rm1b_dummy, sm1, &
                     &                      sm1b_dummy, arg10, theta, sigmasplay, bc1, bc2, &
                     &                      numlayers, epse0, guideindices, &
                     &                      retainspacing, rnext1, rnext1b, numnodes, &
                     &                      numguides)
                r0b = r0b + r0b_dummy
                n0b = n0b + n0b_dummy
                s0b = s0b + s0b_dummy
                rm1b = rm1b + rm1b_dummy
                sm1b = sm1b + sm1b_dummy

                !CALL POPREAL4ARRAY(s0, realtype*numnodes/4)
                S0 = S0_hist_in(minorIndex-1,:)
                !CALL POPREAL4ARRAY(r0, realtype*3*numnodes/4)
                r0 = R_final_in(minorIndex-1,:)

                r0b_dummy = 0.0
                dpseudob_dummy = 0.0
                CALL AREAFACTOR_B(r0, r0b_dummy, dpseudo, dpseudob_dummy, nuarea, &
                     &            numareapasses, bc1, bc2, guideindices, &
                     &            numguides, numnodes, s0, s0b, maxstretch)
                r0b = r0b + r0b_dummy
                dpseudob = dpseudob + dpseudob_dummy

                !CALL POPREAL4ARRAY(sm1, realtype*numnodes/4)
                Sm1 = Sm1_hist_in(minorIndex-1,:)
                !CALL POPREAL4ARRAY(rm1, realtype*3*numnodes/4)
                rm1 = R_final_in(m1Index,:)

                rm1b_dummy = 0.0
                dpseudob_dummy = 0.0
                CALL AREAFACTOR_B(rm1, rm1b_dummy, dpseudo, dpseudob_dummy, nuarea, &
                     &            numareapasses, bc1, bc2, guideindices, &
                     &            numguides, numnodes, sm1, sm1b, maxstretch)
                rm1b = rm1b + rm1b_dummy
                dpseudob = dpseudob + dpseudob_dummy

                rnextb = 0.0
                nnextb = 0.0
             END DO
             !CALL POPREAL4ARRAY(dpseudo, realtype/4)
             db = db + dpseudob/cfactor
             !CALL POPINTEGER4ARRAY(cfactor, inttype/4)
             !CALL POPREAL4ARRAY(s0, realtype*numnodes/4)
             rnextb = rnextb + r0b
             nnextb = nnextb + n0b
          END DO
          rnextb = rnextb + rb(1, :)
          !rb(1, :) = 0.0

          !CALL POPCONTROL1B(branch)
          !branch = retainSpacing

          !IF (branch .EQ. 0) THEN
          IF (retainSpacing) THEN
             rstartb = 0.0
             DO intervalid=numguides+1,1,-1

                node1 = breaknodes(intervalid)
                node2 = breaknodes(intervalid+1)

                numininterval = node2 - node1 + 1

                !CALL POPREAL4ARRAY(normarclength(node1:node2), realtype*(node2-node1+1)/4)
                normArcLength_dummy = normArcLength

                CALL COMPUTE_ARC_LENGTH_B(rstart(3*(node1-1)+1:3*node2), &
                     &                    rstartb(3*(node1-1)+1:3*node2), &
                     &                    numininterval, &
                     &                    normarclength(node1:node2), &
                     &                    normarclengthb(node1:node2))
             END DO
          ELSE
             rstartb = 0.0
          END IF

          rnextb = rnextb + rm1b

          dstartb = db

          !CALL POPCONTROL1B(branch)
          !IF (branch .NE. 0) THEN
          IF (extension_given) THEN

             dstartb_dummy = 0.0
             dmaxb = 0.0
             CALL FINDRATIO_B(dmax, dmaxb, dstart, dstartb_dummy, numlayers, ratioguess, &
                  &           dgrowth, dgrowthb)
             dstartb = dstartb + dstartb_dummy
             radiusb = (marchparameter-1.)*dmaxb

             !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
             rNext = R_final_in(1,:)

             rnextb_dummy = 0.0
             CALL FINDRADIUS_B(rnext, rnextb_dummy, numnodes, radius, radiusb)
             rnextb = rnextb + rnextb_dummy
          END IF

          !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
          rNext = R_final_in(1,:)
          !CALL POPREAL4ARRAY(rnext, realtype*3*numnodes/4)
          NNext = N_final_in(1,:,:)

          rstartb_dummy = 0.0
          CALL PY_PROJECTION_B(rstart, rstartb_dummy, rnext, rnextb, nnext, nnextb, &
               &               projID, numnodes)
          rstartb = rstartb + rstartb_dummy
          projID = projID - 1

        END SUBROUTINE MARCH_B

        !=======================================================================

        subroutine releaseMemory()

          ! This subroutine just deallocates memory used by the marching code.
          ! Remember to call this function after you copy the outputs in Python.

          implicit none

          ! Deallocate output variables
          if (allocated(R_initial_march)) then
            deallocate(R_initial_march, R_smoothed, &
                       R_projected, R_remeshed, R_final, &
                       Sm1_hist, S0_hist, N_projected, N_final)
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

        call computeMatrices_main(r0, N0, S0, rm1, Sm1, layerIndex, theta,&
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
             &   bc2, guideindices, numguides, n, s, sd, maxstretch)

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

        call areafactor_d(r0, r0d, d, dd, nuarea, numareapasses, bc1, &
             &   bc2, guideindices, numguides, n, s, sd, maxstretch)

        end subroutine

        !=================================================================

        subroutine areafactor_test_b(r0, r0b, d, db, nuarea, numareapasses, bc1, &
             &   bc2, guideindices, numguides, n, s, sb, maxstretch)

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

        r0b = 0.0
        db = 0.0

        call areafactor_b(r0, r0b, d, db, nuarea, numareapasses, bc1, &
             &   bc2, guideindices, numguides, n, s, sb, maxstretch)

        end subroutine areafactor_test_b

        !=================================================================

        subroutine findradius_test(numNodes, r0, radius)

        use hypsurfmain, only: findRadius
        implicit none
        integer(kind=inttype), intent(in) :: numNodes
        real(kind=realtype), intent(in) :: r0(3*numNodes)
        real(kind=realtype), intent(out) :: radius

        call findRadius(r0, numNodes, radius)

        end subroutine
        
        !=================================================================

        subroutine findradius_test_d(numNodes, r0, r0d, radius, radiusd)

        use hypsurfmain_d, only: findRadius_d
        implicit none
        integer(kind=inttype), intent(in) :: numNodes
        real(kind=realtype), intent(in) :: r0(3*numNodes)
        real(kind=realtype), intent(in) :: r0d(3*numNodes)
        real(kind=realtype), intent(out) :: radius
        real(kind=realtype), intent(out) :: radiusd

        call findRadius_d(r0, r0d, numNodes, radius, radiusd)

        end subroutine

        !=================================================================

        subroutine findradius_test_b(numNodes, r0, r0b, radius, radiusb)

        use hypsurfmain_b, only: findRadius_b
        implicit none
        integer(kind=inttype), intent(in) :: numNodes
        real(kind=realtype), intent(in) :: r0(3*numNodes)
        real(kind=realtype), intent(out) :: r0b(3*numNodes)
        real(kind=realtype), intent(in) :: radius
        real(kind=realtype), intent(in) :: radiusb

        r0b = 0.0
        call findRadius_b(r0, r0b, numNodes, radius, radiusb)

        end subroutine

        !=================================================================

        subroutine dissipationcoefficients_test(layerindex, r0_xi, r0_eta, dsensor, &
             angle, numlayers, epse0, epse, epsi)

          implicit none

          integer(kind=inttype), intent(in) :: layerindex, numlayers
          real(kind=realtype), intent(in) :: dsensor, angle, epse0, r0_xi(3), &
               &   r0_eta(3)
          real(kind=realtype), intent(out) :: epse, epsi
          
          call dissipationcoefficients(layerindex, r0_xi, r0_eta, dsensor, &
               angle, numlayers, epse0, epse, epsi)
          
        end subroutine dissipationcoefficients_test
        
        !=================================================================
        
        subroutine dissipationcoefficients_d_test(layerindex, r0_xi, r0_xid, r0_eta&
             &   , r0_etad, dsensor, dsensord, angle, angled, numlayers, epse0, epse&
             &   , epsed, epsi, epsid)
          
          use hypsurfmain_d, only: dissipationcoefficients_d
          implicit none
          
          integer(kind=inttype), intent(in) :: layerindex, numlayers
          real(kind=realtype), intent(in) :: dsensor, angle, epse0, r0_xi(3), &
               &   r0_eta(3)
          real(kind=realtype), intent(in) :: dsensord, angled, r0_xid(3), &
               &   r0_etad(3)
          real(kind=realtype), intent(out) :: epse, epsi
          real(kind=realtype), intent(out) :: epsed, epsid
          
          call dissipationcoefficients_d(layerindex, r0_xi, r0_xid, r0_eta&
               &   , r0_etad, dsensor, dsensord, angle, angled, numlayers, epse0, epse&
               &   , epsed, epsi, epsid)
          
        end subroutine dissipationcoefficients_d_test
        
        !=================================================================

        subroutine dissipationcoefficients_b_test(layerindex, r0_xi, r0_xib, r0_eta&
             &   , r0_etab, dsensor, dsensorb, angle, angleb, numlayers, epse0, epse&
             &   , epseb, epsi, epsib)

          use hypsurfmain_b, only: dissipationcoefficients_b
          implicit none

          integer(kind=inttype), intent(in) :: layerindex, numlayers
          real(kind=realtype), intent(in) :: dsensor, angle, epse0, r0_xi(3), &
               &   r0_eta(3)
          real(kind=realtype), intent(out) :: dsensorb, angleb, r0_xib(3), r0_etab(3)
          real(kind=realtype), intent(in) :: epse, epsi
          real(kind=realtype), intent(in) :: epseb, epsib
          
          
          
          call dissipationcoefficients_b(layerindex, r0_xi, r0_xib, r0_eta&
               &   , r0_etab, dsensor, dsensorb, angle, angleb, numlayers, epse0, epse&
               &   , epseb, epsi, epsib)
          
        end subroutine dissipationcoefficients_b_test

        !=================================================================

        subroutine matrixBuilder_test(curr_index, bc1, bc2, r0, rm1, N0, S0, Sm1, &
             numLayers, epsE0, layerIndex, theta, numNodes, K, f)

          use hypsurfMain, only: matrixBuilder
          implicit none

          integer(kind=intType), intent(in) :: curr_index, numNodes, numLayers, layerIndex
          real(kind=realType), intent(in) :: r0(3*numNodes), rm1(3*numNodes), N0(3, numNodes)
          real(kind=realType), intent(in) :: Sm1(numNodes), epsE0, S0(numNodes), theta
          character*32, intent(in) :: bc1, bc2
          real(kind=realType), intent(inout) :: K(3*numNodes, 3*numNodes), f(3*numNodes)
        
          call matrixBuilder(curr_index, bc1, bc2, r0, rm1, N0, S0, Sm1, &
               numLayers, epsE0, layerIndex, theta, numNodes, K, f)

        end subroutine matrixBuilder_test

        !=================================================================

        subroutine matrixBuilder_d_test(curr_index, bc1, bc2, r0, r0d, rm1, rm1d, &
             &   n0, n0d, s0, s0d, sm1, sm1d, numlayers, epse0, layerindex, theta, &
             &   numnodes, k, kd, f, fd)

          use hypsurfMain_d, only: matrixBuilder_d
          implicit none

          integer(kind=inttype), intent(in) :: curr_index, numnodes, numlayers, layerindex
          real(kind=realtype), intent(in) :: r0(3*numnodes), rm1(3*numnodes), n0(3, numnodes)
          real(kind=realtype), intent(in) :: r0d(3*numnodes), rm1d(3*numnodes), n0d(3, numnodes)
          real(kind=realtype), intent(in) :: sm1(numnodes), epse0, s0(numnodes), theta
          real(kind=realtype), intent(in) :: sm1d(numnodes), s0d(numnodes)
          character(len=32), intent(in) :: bc1, bc2
          real(kind=realtype), intent(inout) :: k(3*numnodes, 3*numnodes), f(3*numnodes)
          real(kind=realtype), intent(inout) :: kd(3*numnodes, 3*numnodes), fd(3*numnodes)
        
          call matrixbuilder_d(curr_index, bc1, bc2, r0, r0d, rm1, rm1d, &
               &   n0, n0d, s0, s0d, sm1, sm1d, numlayers, epse0, layerindex, theta, &
               &   numnodes, k, kd, f, fd)

        end subroutine matrixBuilder_d_test

      end module
      
