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

      contains

        !=================================================================

        !=======================================================================

        subroutine march(py_projection, rStart, dStart, theta, sigmaSplay, bc1, bc2, plotQuality,&
        epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, R_initial_march, R_smoothed, R_final, fail, ratios)

        implicit none

        external py_projection
        real(kind=realType), intent(in) :: rStart(3*numNodes),dStart, theta, sigmaSplay
        integer(kind=intType), intent(in) :: numNodes, numLayers, numAreaPasses
        real(kind=realType), intent(in) :: epsE0, extension, nuArea
        character*32, intent(in) :: bc1, bc2
        real(kind=realType), intent(in) :: alphaP0, ratioGuess, cMax
        integer(kind=intType), intent(in) :: numSmoothingPasses, plotQuality

        real(kind=realType), intent(out) :: R_initial_march(numLayers, 3*numNodes)
        real(kind=realType), intent(out) :: R_smoothed(numLayers, 3*numNodes)
        real(kind=realType), intent(out) :: R_final(numLayers, 3*numNodes)
        integer(kind=intType), intent(out) :: fail
        real(kind=realType), intent(out) :: ratios(numLayers-1, numNodes-1)

        real(kind=realType) :: rNext(3*numNodes), NNext(3, numNodes)

        ! Dummy call to get f2py to work correctly
        call py_projection(rStart, rNext, NNext, numNodes)

        call march_main(py_projection, rStart, dStart, theta, sigmaSplay, bc1, bc2, plotQuality,&
        epsE0, alphaP0, extension, nuArea, ratioGuess, cMax, numSmoothingPasses, numAreaPasses,&
        numLayers, numNodes, R_initial_march, R_smoothed, R_final, fail, ratios)

        end subroutine march

        !=======================================================================

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

        !=================================================================

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
