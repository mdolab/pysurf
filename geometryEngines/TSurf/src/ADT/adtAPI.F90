!
!     ******************************************************************
!     *                                                                *
!     * File:          adtAPI.f90                                      *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-09-2006                                      *
!     * Last modified: 02-09-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtAPI
!
!     ******************************************************************
!     *                                                                *
!     * Module, which defines the API of the ADT routines. It is       *
!     * included in a module, such that an explicit interface is       *
!     * present.                                                       *
!     *                                                                *
!     ******************************************************************
!
      use adtBuild
      use adtSearch
      use adtUtils
      implicit none

      !=================================================================

      contains

        !===============================================================

        subroutine adtBuildSurfaceADT(nTria, nQuads,   nNodes,    &
                                      coor,  triaConn, quadsConn, &
                                      BBox,  useBBox,  comm,      &
                                      adtID)

!
!       ****************************************************************
!       *                                                              *
!       * This routine builds the 6 dimensional ADT, which stores the  *
!       * given surface grid. The memory intensive part of these       *
!       * arguments, the arrays with the coordinates and               *
!       * connectivities, are not copied. Instead pointers are set to  *
!       * these arrays. It is therefore the responsibility of the user *
!       * not to deallocate this memory before all the searches have   *
!       * been performed.                                              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nNodes:    Number of local nodes in the given grid.          *
!       * nTria:     Idem for the triangles.                           *
!       * nQuads:    Idem for the quadrilaterals.                      *
!       * BBox(3,2): The possible bounding box. Only elements within   *
!       *            this box will be stored in the ADT.               *
!       * useBBox:   Whether or not to use the bounding box.           *
!       * comm:      MPI-communicator for the global ADT.              *
!       * adtID:     The ID of the ADT.                                *
!       *                                                              *
!       * Subroutine intent(in), target arguments.                     *
!       * ----------------------------------------                     *
!       * coor(3,nNodes):      Nodal coordinates of the local grid.    *
!       * triaConn(3,nTria):   Local connectivity of the triangles.    *
!       * quadsConn(4,nQuads): Idem for the quadrilaterals.            *
!       *                                                              *
!       ****************************************************************
!
        implicit none

        !f2py intent(in) nTria, nQuads, nNodes, coor, triaConn, quadsConn
        !f2py intent(in) BBox, useBBox, comm, adtID
        !f2py depends(nNodes) coor
        !f2py depends(nTria) triaConn
        !f2py depends(nQuads) quadsConn

!
!       Subroutine arguments.
!
        integer, intent(in)          :: comm
        character(len=32), intent(in) :: adtID

        integer(kind=intType), intent(in) :: nTria
        integer(kind=intType), intent(in) :: nQuads
        integer(kind=intType), intent(in) :: nNodes

        logical, intent(in) :: useBBox

        integer(kind=intType), dimension(3,nTria), intent(in) :: triaConn
        integer(kind=intType), dimension(4,nQuads), intent(in) :: quadsConn

        real(kind=realType), dimension(3,2), intent(in) :: BBox

        real(kind=realType), dimension(3,nNodes), intent(in) :: coor

        !===============================================================

        ! Call the subroutine buildSurfaceADT to do the actual work.

        call buildSurfaceADT(nTria,    nQuads,    nNodes, coor,    &
                             triaConn, quadsConn, BBox,   useBBox, &
                             comm,     adtID)

        end subroutine adtBuildSurfaceADT

        !***************************************************************
        !***************************************************************

        subroutine adtBuildVolumeADT(nTetra,    nPyra,    nPrisms,    &
                                     nHexa,     nNodes,   coor,       &
                                     tetraConn, pyraConn, prismsConn, &
                                     hexaConn,  BBox,     useBBox,    &
                                     comm,      adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine builds the 6 dimensional ADT, which stores the  *
!       * given volume grid. The memory intensive part of these        *
!       * arguments, the arrays with the coordinates and               *
!       * connectivities, are not copied. Instead pointers are set to  *
!       * these arrays. It is therefore the responsibility of the user *
!       * not to deallocate this memory before all the searches have   *
!       * been performed.                                              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nNodes:    Number of local nodes in the given grid.          *
!       * nTetra:    Idem for the tetrahedra.                          *
!       * nPyra:     Idem for the pyramids.                            *
!       * nPrisms:   Idem for the prisms.                              *
!       * nHexa:     Idem for the hexahedra.                           *
!       * BBox(3,2): The possible bounding box. Only elements within   *
!       *            this box will be stored in the ADT.               *
!       * useBBox:   Whether or not to use the bounding box.           *
!       * comm:      MPI-communicator for the global ADT.              *
!       * adtID:     The ID of the ADT.                                *
!       *                                                              *
!       * Subroutine intent(in), target arguments.                     *
!       * ----------------------------------------                     *
!       * coor(3,nNodes):        Nodal coordinates of the local grid.  *
!       * tetraConn(4,nTetra):   Local connectivity of the tetrahedra. *
!       * pyraConn(5,nPyra):     Idem for the pyramids.                *
!       * prismsConn(6,nPrisms): Idem for the prisms.                  *
!       * hexaConn(8,nHexa):     Idem for the hexahedra.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer, intent(in)          :: comm
        character(len=*), intent(in) :: adtID

        integer(kind=intType), intent(in) :: nTetra
        integer(kind=intType), intent(in) :: nPyra
        integer(kind=intType), intent(in) :: nPrisms
        integer(kind=intType), intent(in) :: nHexa
        integer(kind=intType), intent(in) :: nNodes

        logical, intent(in) :: useBBox

        integer(kind=intType), dimension(:,:), intent(in) :: tetraConn
        integer(kind=intType), dimension(:,:), intent(in) :: pyraConn
        integer(kind=intType), dimension(:,:), intent(in) :: prismsConn
        integer(kind=intType), dimension(:,:), intent(in) :: hexaConn

        real(kind=realType), dimension(3,2), intent(in) :: BBox

        real(kind=realType), dimension(:,:), intent(in) :: coor

        !===============================================================

        ! Call the subroutine buildVolumeADT to do the actual work.

        call buildVolumeADT(nTetra,     nPyra,     nPrisms,   nHexa,    &
                            nNodes,     coor,      tetraConn, pyraConn, &
                            prismsConn, hexaConn,  BBox,      useBBox,  &
                            comm,       adtID)

        end subroutine adtBuildVolumeADT

        !***************************************************************
        !***************************************************************

        subroutine adtContainmentSearch(nCoor,     coor,        adtID,     &
                                        procID,    elementType, elementID, &
                                        uvw,       allxfs,      nInterpol, &
                                        arrDonor,  arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found the corresponding entry in procID is  *
!       * set to -1 to indicate failure.                               *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates of these points.                  *
!       * adtID:     The ADT to be searched.                           *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point.                       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                            elementType
        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: allxfs
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol

        !===============================================================

        ! Call the subroutine containmentSearch to do the actual work.

        call containmentSearch(nCoor,       coor,      adtID,    procID,    &
                               elementType, elementID, uvw,      allxfs,    &
                               nInterpol, arrDonor, arrInterpol)

        end subroutine adtContainmentSearch

        !***************************************************************
        !***************************************************************

        subroutine adtDeallocateADTs(adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine deallocates the memory for the given entry in   *
!       * the array ADTs and it tries to reallocate ADTs itself        *
!       * accordingly.                                                 *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * adtID: The entry in ADTs to be deallocated.                  *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        character(len=*), intent(in) :: adtID

        !===============================================================

        ! Call the subroutine deallocateADTs to do the actual work.

        call deallocateADTs(adtID)

        end subroutine adtDeallocateADTs

        !***************************************************************
        !***************************************************************

        subroutine adtFailSafeSearch(nCoor,    coor,        adtID,     &
                                     procID,   elementType, elementID, &
                                     uvw,      dist2,       allxfs,    &
                                     nInterpol,   arrDonor,  &
                                     arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found a minimum distance search is          *
!       * performed, such that always an interpolation can be          *
!       * performed. To indicate that the element does not contain the *
!       * point the element ID is negated.                             *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of coordinates for which the element must be   *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element, i.e. if  *
!       *              a minimum distance search had to be used.       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                              elementType

        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol

        real(kind=realType), dimension(:), intent(inout) :: dist2
        real(kind=realType), dimension(:,:), intent(inout) :: allxfs


        !===============================================================

        ! Call the subroutine failSafeSearch to do the actual work.

        call failSafeSearch(nCoor,       coor,      adtID,     procID,   &
                            elementType, elementID, uvw,       dist2,    &
                            allxfs,      nInterpol, arrDonor, &
                            arrInterpol)

        end subroutine adtFailSafeSearch

        !***************************************************************
        !***************************************************************

        subroutine adtMinDistanceSearch(nCoor,     nNodes,    coor,        &
                                        adtID,     procID,    elementType, &
                                        elementID, uvw,       dist2,       &
                                        allxfs,    nInterpol,   &
                                        arrDonor,  arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT which minimizes the distance to     *
!       * this point.                                                  *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of points for which the element must be        *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element.          *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights. If the tree       *
!       *              corresponds to a surface mesh the third entry   *
!       *              of this array will not be filled.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       * arrInterpol: Array with the interpolated data. Values will   *
!       *              only be replaced if we find a better candidate  *
!       *              at the current surface.                         *
!       *                                                              *
!       * allxfs: Array with projected points. Values will only be     *
!       *         replaced if we find a better candidate  at the       *
!       *         current surface.                                     *
!       ****************************************************************
!
        implicit none

        !f2py intent(in) nCoor, nNodes, coor, adtID, nInterpol, arrDonor
        !f2py intent(out) procID, elementType, elementID, uvw, arrInterpol
        !f2py intent(inout) dist2, allxfs, arrInterpol

!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nNodes,nInterpol
        character(len=32),     intent(in) :: adtID

        real(kind=realType), dimension(3,nCoor), intent(in) :: coor
        real(kind=realType), dimension(nInterpol,nNodes), intent(in) :: arrDonor

        integer, dimension(nCoor), intent(out) :: procID
        integer(kind=intType), dimension(nCoor), intent(out) :: elementID

        integer(kind=adtElementType), dimension(nCoor), intent(out) :: &
                                                            elementType

        real(kind=realType), dimension(3,nCoor), intent(out) :: uvw

        real(kind=realType), dimension(nCoor), intent(inout) :: dist2
        real(kind=realType), dimension(3,nCoor), intent(inout) :: allxfs
        real(kind=realType), dimension(nInterpol,nCoor), intent(inout) :: arrInterpol

        !===============================================================

        ! Call the subroutine minDistanceSearch to do the actual work.

        call minDistanceSearch(nCoor,       coor,      adtID,     procID,   &
                               elementType, elementID, uvw,       dist2,    &
                               allxfs,      nInterpol, arrDonor, &
                               arrInterpol)

        end subroutine adtMinDistanceSearch

        !***************************************************************
        !***************************************************************

        subroutine computeNodalNormals(nCoor, nTria, nQuads, coor, triaConn, &
                                       quadsConn, nodalNormals)
!
!       ****************************************************************
!       *                                                              *
!       *                                                              *
!       ****************************************************************
!
        implicit none

        !f2py intent(in) nCoor, nTria, nQuads, coor, triaConn, quadsConn
        !f2py intent(out) nodalNormals

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
        real(kind=realType), dimension(nCoor) :: connect_count
        real(kind=realType) :: normal1(3), normal2(3), normal3(3), normal4(3)
        integer(kind=intType) :: i, ind1, ind2, ind3, ind4
        real(kind=realType) :: x1(3), x2(3), x3(3), x4(3)
        real(kind=realType) :: x12(3), x23(3), x34(3), x41(3)


        !===============================================================

        ! Loop over triangle connectivities
        do i = 1,nTria
          ind1 = triaConn(1, i)
          ind2 = triaConn(2, i)
          ind3 = triaConn(3, i)

          x1 = coor(:, ind1)
          x2 = coor(:, ind2)
          x3 = coor(:, ind3)

          x12 = x1 - x2
          x23 = x2 - x3

          call cross_product(x12, x23, normal1)
          normal1 = normal1 / sqrt(dot_product(normal1, normal1))

          nodalNormals(:, ind1) = nodalNormals(:, ind1) + normal1
          nodalNormals(:, ind2) = nodalNormals(:, ind2) + normal1
          nodalNormals(:, ind3) = nodalNormals(:, ind3) + normal1

          connect_count(ind1) = connect_count(ind1) + 1
          connect_count(ind2) = connect_count(ind2) + 1
          connect_count(ind3) = connect_count(ind3) + 1
        end do

        ! Loop over quad connectivities
        do i = 1,nQuads
          ind1 = quadsConn(1, i)
          ind2 = quadsConn(2, i)
          ind3 = quadsConn(3, i)
          ind4 = quadsConn(4, i)

          x1 = coor(:, ind1)
          x2 = coor(:, ind2)
          x3 = coor(:, ind3)
          x4 = coor(:, ind4)

          x12 = x1 - x2
          x23 = x2 - x3
          x34 = x3 - x4
          x41 = x4 - x1

          call cross_product(x12, -x41, normal1)
          normal1 = normal1 / sqrt(dot_product(normal1, normal1))
          call cross_product(x23, -x12, normal2)
          normal2 = normal2 / sqrt(dot_product(normal2, normal2))
          call cross_product(x34, -x23, normal3)
          normal3 = normal3 / sqrt(dot_product(normal3, normal3))
          call cross_product(x41, -x34, normal4)
          normal4 = normal4 / sqrt(dot_product(normal4, normal4))

          nodalNormals(:, ind1) = nodalNormals(:, ind1) + normal1
          nodalNormals(:, ind2) = nodalNormals(:, ind2) + normal2
          nodalNormals(:, ind3) = nodalNormals(:, ind3) + normal3
          nodalNormals(:, ind4) = nodalNormals(:, ind4) + normal4

          connect_count(ind1) = connect_count(ind1) + 1
          connect_count(ind2) = connect_count(ind2) + 1
          connect_count(ind3) = connect_count(ind3) + 1
          connect_count(ind4) = connect_count(ind4) + 1
        end do

        do i = 1,3
            nodalNormals(i, :) = nodalNormals(i, :) / connect_count
        end do

        do i=1,nCoor
          normal1 = nodalNormals(:, i)
          nodalNormals(:, i) = normal1 / sqrt(dot_product(normal1, normal1))
        end do

        end subroutine computeNodalNormals

        subroutine cross_product(A, B, C)

          implicit none

          real(kind=realType), intent(in) :: A(3), B(3)
          real(kind=realType), intent(out) :: C(3)

          C(1) = A(2) * B(3) - A(3) * B(2)
          C(2) = A(3) * B(1) - A(1) * B(3)
          C(3) = A(1) * B(2) - A(2) * B(1)

        end subroutine cross_product


        !***************************************************************
        !***************************************************************

        subroutine adtGetNumberOfTrees(nADT)
!
!       ****************************************************************
!       * This subroutine just gets the maximum number of ADTs         *
!       * defined so far                                               *
!       ****************************************************************
!
! Ney Secco 08-2016

        implicit none

        !f2py intent(out) nADT

        ! OUTPUT VARIABLES
        integer(kind=intType), intent(out) :: nADT

        ! call subroutine defined in adtUtils.F90 that will do all the job
        call numberOfADTs(nADT)

      end subroutine adtGetNumberOfTrees

      end module adtAPI
