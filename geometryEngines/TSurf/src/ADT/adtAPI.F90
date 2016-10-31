!
!     ******************************************************************
!     *                                                                *
!     * File:          adtAPI.f90                                      *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-09-2006                                      *
!     * Last modified: 10-05-2016                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtAPI
!
!     ******************************************************************
!     *                                                                *
!     * Module which defines the API of the ADT routines. It is        *
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

        subroutine adtIntersectionSearch(nBBox,     inpBBox,        adtID,     &
                                         procID,    elementType, elementID, &
                                         BBoxPtr)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the detects elements whose bounding    *
!       * boxes intersect the user provided bounding boxes (inpBBox).  *
!       * This narrows down, for instance, candidates in intersection  *
!       * algorithms.                                                  *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nBBox:     Number of bounding boxes for which the            *
!       *            intersecting elements must be determined.         *
!       * inpBBox:   real(6,nBBox) containing the corners of each      *
!       *            bounding box (xmin, ymin, zmin, xmax, ymax, zmax).*
!       * adtID:     The ADT to be searched.                           *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the elements that intersect the given     *
!       *              bounding boxes are stored. You need to use      *
!       *              BBoxPtr to slice this array and get elements    *
!       *              that intersect a specific bounding box.         *
!       * elementType: The type of element which intersects the        *
!       *              bounding boxes.                                 *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which intersects the bounding boxes.            *
!       *              The sizes of procID, elementType, and elementID *
!       *              are all the same, and equal to the total number *
!       *              of elements that intersect any given BBox.      *
!       * BBoxPtr:     integer(nBBox+1). Pointers used to slice procID,*
!       *              elementType, and elementID.                     *
!       *              The indices from BBoxPtr(i) to BBoxPtr(i+1)-1   *
!       *              belong to inpBBox(i).                           *
!       *              If BBoxPtr(i) == BBoxPtr(i+1), then inpBBox(i)  *
!       *              has no intersections with the tree elements.    *
!       *                                                              *
!       ****************************************************************
!
!       ATTENTION!: This subroutine will only work if the trees in every
!       proc are the same. In other words, all elements should be defined
!       in all processors! In order to do this, you should use the same
!       coor, and connectivities when building the trees in each proc.
!
!       Ney Secco 2016-09
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nBBox

        real(kind=realType), dimension(6,nBBox), intent(in) :: inpBBox

        character(len=*), intent(in) :: adtID

        integer(kind=intType), dimension(:), allocatable, intent(out) :: procID

        integer(kind=adtElementType), dimension(:), allocatable, intent(out) :: elementType

        integer(kind=intType), dimension(:), allocatable, intent(out) :: elementID

        integer(kind=intType), dimension(nBBox+1), intent(out) :: BBoxPtr

        !===============================================================

        ! Call the subroutine intersectionSearch to do the actual work.

        call intersectionSearch(nBBox,       inpBBox,       adtID, &
                                procID,  elementType,   elementID, &
                                BBoxPtr)

      end subroutine adtIntersectionSearch

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
!       * nNodes: Number of donor nodes for interpolation information.
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

        subroutine adtMinDistanceSearch_d(nCoor,       nNodes,      coor,        coord, &
                                          adtID,       adtCoord,    procID,      elementType, &
                                          elementID,   uvw,         dist2,       &
                                          allxfs,      allxfsd,     nInterpol,   &
                                          arrDonor,    arrDonord, &
                                          arrInterpol, arrInterpold)
!
!       ****************************************************************
!       *                                                              *
!       * This routine computes derivatives of the minimum distance    *
!       * (projection) algorithm. Note that most inputs should be      *
!       * obtained with the execution of the original routine first.   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of points for which the element must be        *
!       *        determined.                                           *
!       * nNodes: Number of donor nodes for interpolation information. *
!       * coor: The coordinates of these points.                       *
!       * coord: Derivative seeds for the coordinates.                 *
!       * adtID: The ADT to be searched.                               *
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
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT.                                  *
!       * allxfs: Array with projected points.                         *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       * arrDonord: Derivative seeds of the donor data.               *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * allxfsd: Derivative seeds of projected points.               *
!       * arrInterpold: Derivative seeds of the interpolated data.     *
!       *                                                              *
!       ****************************************************************
!       Ney Secco 2016-10
!
        implicit none

        !f2py intent(in) nCoor, nNodes, coor, coord, adtID, nInterpol, arrDonor
        !f2py intent(in) procID, elementType, elementID, uvw, arrInterpol
        !f2py intent(in) dist2, allxfs, arrInterpol
        !f2py intent(out) allxfsd, arrInterpold

!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nNodes,nInterpol
        character(len=32),     intent(in) :: adtID

        real(kind=realType), dimension(3,nCoor), intent(in) :: coor, coord
        real(kind=realType), dimension(3,nNodes), intent(in) :: adtCoord
        real(kind=realType), dimension(nInterpol,nNodes), intent(in) :: arrDonor, arrDonord

        integer, dimension(nCoor), intent(in) :: procID
        integer(kind=intType), dimension(nCoor), intent(in) :: elementID

        integer(kind=adtElementType), dimension(nCoor), intent(in) :: &
                                                            elementType

        real(kind=realType), dimension(3,nCoor), intent(in) :: uvw

        real(kind=realType), dimension(nCoor), intent(in) :: dist2
        real(kind=realType), dimension(3,nCoor), intent(in) :: allxfs
        real(kind=realType), dimension(nInterpol,nCoor), intent(in) :: arrInterpol

        real(kind=realType), dimension(3,nCoor), intent(out) :: allxfsd
        real(kind=realType), dimension(nInterpol,nCoor), intent(out) :: arrInterpold

        !===============================================================

        ! Call the subroutine minDistanceSearch_d to do the actual work.

        call minDistanceSearch_d(nCoor,       coor,        coord,     adtID,    adtCoord, &
                                 procID,      elementType, elementID, uvw,      dist2,    &
                                 allxfs,      allxfsd,     nInterpol, arrDonor, arrDonord, &
                                 arrInterpol, arrInterpold)

        end subroutine adtMinDistanceSearch_d

        !***************************************************************
        !***************************************************************

        subroutine adtComputeNodalNormals(nCoor, nTria, nQuads, coor, triaConn, &
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
        use adtProjections
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

        ! Call function that will do the job.
        ! computeNodalNormals defined in adtProjections.F90
        call computeNodalNormals(nCoor, nTria, nQuads, coor, triaConn, &
                                 quadsConn, nodalNormals)

      end subroutine adtComputeNodalNormals

        !***************************************************************
        !***************************************************************

        subroutine adtComputeNodalNormals_d(nCoor, nTria, nQuads, coor, coord, triaConn, &
                                            quadsConn, nodalNormals, nodalNormalsd)
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
!       * coord: Derivative seed of the nodal coordinates.             *
!       * triaConn: Connectivity information for the triangle elments. *
!       * quadsConn: Connectivity information for the quad elments.    *
!       * nodalNormals: The interpolated normals at each coordinate    *
!       *               node (obtained by the original function).      *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * nodalNormalsd: Derivatives of the interpolated normals.      *
!       *                                                              *
!       ****************************************************************
!
        use adtProjections_d
        implicit none

        !f2py intent(in) nCoor, nTria, nQuads, coor, triaConn, quadsConn
        !f2py intent(out) nodalNormalsd, nodalNormals

!
!       Subroutine arguments.
!
        ! Input
        integer(kind=intType), intent(in) :: nCoor, nTria, nQuads
        real(kind=realType), dimension(3,nCoor), intent(in) :: coor
        real(kind=realType), dimension(3,nCoor), intent(in) :: coord
        integer(kind=intType), dimension(3,nTria), intent(in) :: triaConn
        integer(kind=intType), dimension(4,nQuads), intent(in) :: quadsConn
        real(kind=realType), dimension(3,nCoor), intent(out) :: nodalNormals

        ! Output
        real(kind=realType), dimension(3,nCoor), intent(out) :: nodalNormalsd

        ! Call function that will do the job.
        ! computeNodalNormals defined in adtProjections.F90
        call computeNodalNormals_d(nCoor, nTria, nQuads, coor, coord, triaConn, &
                                   quadsConn, nodalNormals, nodalNormalsd)

      end subroutine adtComputeNodalNormals_d

        !***************************************************************
        !***************************************************************
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

        ! call subroutine defined in adtUtils.F90 that will do the whole job
        call numberOfADTs(nADT)

      end subroutine adtGetNumberOfTrees

      end module adtAPI
