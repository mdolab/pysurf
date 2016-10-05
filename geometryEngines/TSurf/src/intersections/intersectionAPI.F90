module intersectionAPI

  use precision
  implicit none

  ! Here we define output variables of unknown shape as attributes, so we
  ! can read them from Python

  ! Intersection FE data
  ! coor: real(3,nNodesInt) -> Coordinates of the nodes that define the intersection.
  real(kind=realType), dimension(:,:), allocatable :: coor

  ! barsConn: integer(2,nBarsInt) -> Connectivities of the bar elements that define the intersection.
  integer(kind=intType), dimension(:,:), allocatable :: barsConn

  ! parentTria: integer(2,nBarsInt) -> Triangle pairs (from component A and component B) that
  !                                    generated each intersection bar element.
  integer(kind=intType), dimension(:,:), allocatable :: parentTria

contains

  subroutine computeIntersection(nNodesA, nTriaA, nQuadsA, &
    nNodesB, nTriaB, nQuadsB, &
    coorA, triaConnA, quadsConnA, &
    coorB, triaConnB, quadsConnB, &
    distTol, comm)

    ! This function computes the intersection curve between two
    ! triangulated surface components A and B.
    !
    ! INPUTS
    !
    ! nNodesA: integer -> number of nodes in component A
    !
    ! nTriaA: integer -> number of triangle elements in component A
    !
    ! nQuadsA: integer -> number of quad elements in component A
    !
    ! nNodesB: integer -> number of nodes in component B
    !
    ! nTriaB: integer -> number of triangle elements in component B
    !
    ! nQuadsB: integer -> number of quad elements in component B
    !
    ! coorA: real(3,nNodesA) -> Nodal coordinates used by component A
    !
    ! triaConnA: integer(3,nTriaA) -> Connectivities of triangle elements in component A
    !
    ! quadsConnA: integer(4,nQuadsA) -> Connectivities of quad elements in component A
    !
    ! coorB: real(3,nNodesB) -> Nodal coordinates used by component B
    !
    ! triaConnB: integer(3,nTriaB) -> Connectivities of triangle elements in component B
    !
    ! quadsConnB: integer(4,nQuadsB) -> Connectivities of quad elements in component B
    !
    ! distTol: real -> distance tolerance to merge nearby nodes in
    !          the intersection curve.
    !
    ! comm : integer -> MPI communicator.
    !
    ! OUTPUTS
    ! This subroutine has no explicit outputs. It updates the variables coor, barsConn, and parentTria
    ! which should be called from Python as attributes of the intersectionAPI module.
    !
    ! Ney Secco 2016-08

    use Intersection
    use Utilities ! This will bring condenseBarNodes_main
    use adtAPI ! This will bring adtBuildSurfaceADT and adtIntersectionSearch
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: nNodesA, nTriaA, nQuadsA
    integer(kind=intType), intent(in) :: nNodesB, nTriaB, nQuadsB
    real(kind=realType), dimension(3,nNodesA), intent(in) :: coorA
    integer(kind=intType), dimension(3,nTriaA), intent(in) :: triaConnA
    integer(kind=intType), dimension(4,nQuadsA), intent(in) :: quadsConnA
    real(kind=realType), dimension(3,nNodesB), intent(in) :: coorB
    integer(kind=intType), dimension(3,nTriaB), intent(in) :: triaConnB
    integer(kind=intType), dimension(4,nQuadsB), intent(in) :: quadsConnB
    real(kind=realType), intent(in) :: distTol
    integer(kind=intType), intent(in) :: comm
    !f2py intent(in) nNodesA, nTriaA, nQuadsA
    !f2py intent(in) nNodesB, nTriaB, nQuadsB
    !f2py intent(in) coorA, triaConnA, quadsConnA
    !f2py intent(in) coorB, triaConnB, quadsConnB
    !f2py intent(in) distTol, comm

    ! Working variables
    real(kind=realType), dimension(3,2) :: BBoxA, BBoxB, BBoxAB
    logical :: overlap
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_A, innerQuadsID_A
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_B, innerQuadsID_B
    integer(kind=intType), dimension(:,:), allocatable :: allTriaConnA, allTriaConnB
    integer(kind=intType) :: nInnerTriaA, nInnerQuadsA, nInnerTriaB, nInnerQuadsB
    integer(kind=intType) :: ii, jj
    integer(kind=intType) :: intersect
    integer(kind=intType) :: arraySize, nAllocations, nBarsConn
    integer(kind=intType), dimension(:,:), allocatable :: extBarsConn, extTempInt
    real(kind=realType), dimension(:,:), allocatable :: extCoor, extTempReal
    real(kind=realType), dimension(3) :: node1A, node2A, node3A
    real(kind=realType), dimension(3) :: node1B, node2B, node3B
    real(kind=realType), dimension(3) :: vecStart, vecEnd
    integer(kind=intType) :: nNodesInt, nBarsInt, nUniqueNodes
    integer(kind=intType), dimension(:,:), allocatable :: extParentTria

    character(len=32) :: adtID
    real(kind=realType), dimension(:,:), allocatable :: triaBBox
    real(kind=realType), dimension(6,0) :: quadsBBox
    integer(kind=intType), dimension(:), allocatable :: procID, elementID
    integer(kind=adtElementType), dimension(:), allocatable :: elementType
    integer(kind=intType), dimension(:), allocatable :: BBoxPtr
    integer(kind=intType) :: startIndex, endIndex, numBBoxTest, candidateID

    ! EXECUTION

    ! Compute bounding boxes for each component
    call computeBBox(coorA, BBoxA)
    call computeBBox(coorB, BBoxB)

    ! Compute bounding boxes intersection
    call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

    ! We can stop if there is no bounding box intersection
    if (.not. overlap) then
      print *,'Geometry objects do not overlap.'
      allocate(barsConn(0, 0), coor(0, 0))
      return
    else
      print *,'Geometry objects overlap.'
    end if

    ! Filter elements that are inside the intersected bounding box
    call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
    innerTriaID_A, innerQuadsID_A)
    call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
    innerTriaID_B, innerQuadsID_B)

    ! Get number of inner elements
    nInnerTriaA = size(innerTriaID_A)
    nInnerQuadsA = size(innerQuadsID_A)
    nInnerTriaB = size(innerTriaID_B)
    nInnerQuadsB = size(innerQuadsID_B)

    ! Print log
    print *,'Number of interior elements in A:'
    print *,nInnerTriaA + nInnerQuadsA,'of',nTriaA + nQuadsA
    print *,'Number of interior elements in B:'
    print *,nInnerTriaB + nInnerQuadsB,'of',nTriaB + nQuadsB

    ! Split quads into triangles
    call getAllTrias(triaConnA, quadsConnA, innerTriaID_A, innerQuadsID_A, allTriaConnA)
    call getAllTrias(triaConnB, quadsConnB, innerTriaID_B, innerQuadsID_B, allTriaConnB)

    ! Update number of interior triangles
    nInnerTriaA = size(allTriaConnA,2)
    nInnerTriaB = size(allTriaConnB,2)

    !========================================
    ! ADT filter
    !========================================
    ! This part of the code will generate an ADT using only elements inside the
    ! intersecting bounding box. This will help us eliminate unnecessary
    ! intersection checks.

    ! Create generic ID for the ADT
    adtID = 'intersectionFilter'

    ! Generate ADT with the component that has the largest number of elements.
    ! This ADT will contain only the triangles.
    ! adtBuildSurfaceADT defined in adtAPI.F90
    call adtBuildSurfaceADT(nInnerTriaB, 0, nNodesB, &
         coorB, allTriaConnB, quadsConnB, &
         BBoxAB, .false., comm, &
         adtID)

    ! Allocate arrays for bounding boxes
    allocate(triaBBox(6,nInnerTriaA), BBoxPtr(nInnerTriaA+1))

    ! Now we compute bounding boxes for each element of the component that
    ! does not have a tree. Remember that quadsBBox is a dummy variable, as
    ! we only deal with triangle elements.
    call computeBBoxPerElements(nNodesA, nInnerTriaA, 0, &
                                coorA, allTriaConnA, quadsConnA, &
                                triaBBox, quadsBBox)

    ! Now we use ADT to search for possible intersections.
    ! adtIntersectionSearch defined in adtAPI.F90
    call adtIntersectionSearch(nInnerTriaA, triaBBox, adtID, &
                               procID, elementType, elementID, &
                               BBoxPtr)

    ! Now we can deallocate the tree
    call adtDeallocateADTs(adtID)

    !========================================
    ! Pair-wise intersection tests
    !========================================

    ! Initialize extended connectivity arrays for the intersection curves.
    ! These arrays will be cropped to get the actual outputs.
    ! For now we will allocate arrays of size arraySize. We will increase them
    ! if necessary
    arraySize = 1000
    allocate(extCoor(3,arraySize), extBarsConn(2,arraySize))

    ! Allocate array that will store the IDs of the two triangles that generate
    ! each bar element. This is used by derivatives only.
    allocate(extParentTria(2,arraySize))

    ! Initialize values
    extCoor = 0.0
    extBarsConn = 0
    extParentTria = 0

    ! For now, we have just one allocation
    nAllocations = 1

    ! Initialize the number of intersection connectivities known so far.
    nBarsConn = 0

    ! Compute all triangle-triangle intersections.
    ! We will call every pair between triangles in A and B
    do ii = 1,nInnerTriaA

       ! Get initial and final position to slice tree outputs
       startIndex = BBoxPtr(ii)
       endIndex = BBoxPtr(ii+1)

       ! Check how many bounding boxes intersect the bounding box from
       ! the current element from component A. This is given by the difference
       ! of consecutive pointers in BBoxPtr.
       numBBoxTest = endIndex - startIndex

       if (numBBoxTest .le. 0) then ! We have no intersection, so jump to next element
          cycle
       end if

       ! Get nodes of the element in A
       node1A = coorA(:, allTriaConnA(1,ii))
       node2A = coorA(:, allTriaConnA(2,ii))
       node3A = coorA(:, allTriaConnA(3,ii))

       ! Now loop over the triangles of the other component that might
       ! intersect the current element.
       do jj = 1,numBBoxTest

          ! Get index of element in component B
          candidateID = elementID(startIndex + jj - 1)

          ! Get nodes of the element in B
          node1B = coorB(:, allTriaConnB(1,candidateID))
          node2B = coorB(:, allTriaConnB(2,candidateID))
          node3B = coorB(:, allTriaConnB(3,candidateID))

          ! Call triangle-triangle intersection routine
          call triTriIntersect(node1A, node2A, node3A, &
                               node1B, node2B, node3B, &
                               intersect, vecStart, vecEnd)
          
          ! Check if the triangles actually intersect
          if (intersect .eq. 1) then
             
             ! Increase the number of intersections elements known so far
             nBarsConn = nBarsConn + 1

             ! Check if we already extrapolated the memory allocated so far
             if ((nBarsConn .gt. size(extBarsConn,2)) .or. (2*nBarsConn .gt. size(extCoor,2))) then
                
                ! We need to allocate more memory
                nAllocations = nAllocations + 1
                
                ! REALLOCATING extCoor
                
                ! Create temporary array and initialize
                allocate(extTempReal(3,nAllocations*arraySize))
                extTempReal = 0.0
                
                ! Tranfer data from original array
                extTempReal(:,:(nAllocations-1)*arraySize) = extCoor
                
                ! Now move the new allocation back to extCoor.
                ! extTemp is deallocated in this process.
                call move_alloc(extTempReal, extCoor)
                
                ! REALLOCATING extBarsConn
                
                ! Create temporary array
                allocate(extTempInt(2,nAllocations*arraySize))
                extTempInt = 0
                
                ! Tranfer data from original array
                extTempInt(:,:(nAllocations-1)*arraySize) = extBarsConn
                
                ! Now move the new allocation back to extBarsConn.
                ! extTemp is deallocated in this process.
                call move_alloc(extTempInt, extBarsConn)
                
                ! REALLOCATING extParentTria
                
                ! Create temporary array
                allocate(extTempInt(2,nAllocations*arraySize))
                extTempInt = 0
                
                ! Tranfer data from original array
                extTempInt(:,:(nAllocations-1)*arraySize) = extParentTria
                
                ! Now move the new allocation back to extParentTria.
                ! extTemp is deallocated in this process.
                call move_alloc(extTempInt, extParentTria)
                
             end if
             
             ! Assign new nodes
             extCoor(:,2*nBarsConn-1) = vecStart
             extCoor(:,2*nBarsConn)   = vecEnd
             
             ! Assign new connectivity
             extBarsConn(:,nBarsConn) = [2*nBarsConn-1, 2*nBarsConn]
             
             ! Assign new parent triangles
             extParentTria(:,nBarsConn) = [ii, candidateID]
             
          end if
          
       end do
       
    end do
    
    ! Crop extended connectivity array
    allocate(barsConn(2,nBarsConn))
    barsConn(:,:) = extBarsConn(:,:nBarsConn)
    deallocate(extBarsConn)

    ! Crop extended parents array
    allocate(parentTria(2,nBarsConn))
    parentTria(:,:) = extParentTria(:,:nBarsConn)
    deallocate(extParentTria)

    ! Merge close nodes to get continuous FE data. The condensed coordinates (coor)
    ! will be returned to Python.
    nNodesInt = size(extCoor,2)
    nBarsInt = size(barsConn,2)
    call condenseBarNodes_main(nNodesInt, nBarsInt, distTol, &
    extCoor, barsConn, nUniqueNodes) ! Defined in utilities.F90

    ! Allocate array to hold only unique nodes
    allocate(coor(3,nUniqueNodes))
    coor(:,:) = extCoor(:,1:nUniqueNodes)
    deallocate(extCoor)

  end subroutine computeIntersection

  !=============================================

  subroutine computeIntersection_b(nNodesA, nTriaA, nQuadsA, &
    nNodesB, nTriaB, nQuadsB, &
    nNodesInt, nBarsInt, &
    coorA, coorAb, triaConnA, quadsConnA, &
    coorB, coorBb, triaConnB, quadsConnB, &
    coorInt, coorIntb, barsConnInt, &
    parentTriaInt, distTol)

    ! This function computes the backward derivatives of the intersection curve between two
    ! triangulated surface components A and B.
    ! This functions receives the seeds of the intersection points (coorIntb) and returns the
    ! derivatives of the component points (coorAb, and coorBb)
    !
    ! INPUTS
    !
    ! nNodesA: integer -> number of nodes in component A
    !
    ! nTriaA: integer -> number of triangle elements in component A
    !
    ! nQuadsA: integer -> number of quad elements in component A
    !
    ! nNodesB: integer -> number of nodes in component B
    !
    ! nTriaB: integer -> number of triangle elements in component B
    !
    ! nQuadsB: integer -> number of quad elements in component B
    !
    ! nNodesInt: integer -> number of nodes in the intersection curve
    !
    ! nBarsInt: integer -> number of bar elements in the intersection curve
    !
    ! coorA: real(3,nNodesA) -> Nodal coordinates used by component A
    !
    ! triaConnA: integer(3,nTriaA) -> Connectivities of triangle elements in component A
    !
    ! quadsConnA : integer(4,nQuadsA) -> Connectivities of quad elements in component A
    !
    ! coorB: real(3,nNodesB) -> Nodal coordinates used by component B
    !
    ! triaConnB: integer(3,nTriaB) -> Connectivities of triangle elements in component B
    !
    ! quadsConnB : integer(4,nQuadsB) -> Connectivities of quad elements in component B
    !
    ! coorInt: real(3,nNodesInt) -> Nodal coordinates used by the intersection curve
    !
    ! coorIntb: real(3,nNodesInt) -> Derivative seed of the nodal coordinates used by the intersection curve
    !
    ! barsConnInt: integer(2,nBarsInt) -> Connectivities of bar elements in the intersection curve
    !
    ! parentTriaInt: integer(2,nBarsInt) -> Triangle pairs (from component A and component B) that
    !                                       generated each intersection bar element.
    !
    ! distTol: real -> distance used to check if the start and end nodes of the intersection are flipped.
    !                  use the same distTol used for the original computeIntersection.
    !
    ! OUTPUTS
    !
    ! coorAb: real(3,nNodesA) -> Derivative seeds of the nodal coordinates used by component A
    !
    ! coorBb: real(3,nNodesB) -> Derivative seeds of the nodal coordinates used by component B
    !
    !
    ! Ney Secco 2016-09

    use Intersection
    use Intersection_b, only: triTriIntersect_b
    use Utilities ! This will bring condenseBarNodes_main
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: nNodesA, nTriaA, nQuadsA
    integer(kind=intType), intent(in) :: nNodesB, nTriaB, nQuadsB
    integer(kind=intType), intent(in) :: nNodesInt, nBarsInt
    real(kind=realType), dimension(3,nNodesA), intent(in) :: coorA
    integer(kind=intType), dimension(3,nTriaA), intent(in) :: triaConnA
    integer(kind=intType), dimension(4,nQuadsA), intent(in) :: quadsConnA
    real(kind=realType), dimension(3,nNodesB), intent(in) :: coorB
    integer(kind=intType), dimension(3,nTriaB), intent(in) :: triaConnB
    integer(kind=intType), dimension(4,nQuadsB), intent(in) :: quadsConnB
    real(kind=realType), dimension(3,nNodesInt), intent(in) :: coorInt ! Intersection nodes
    real(kind=realType), dimension(3,nNodesInt), intent(in) :: coorIntb ! Intersection nodes
    integer(kind=intType), dimension(2,nBarsInt), intent(in) :: barsConnInt ! Intersection bar connectivities
    integer(kind=intType), dimension(2,nBarsInt), intent(in) :: parentTriaInt ! Pair of triangles from A and B that generates each bar
    real(kind=realType), intent(in) :: distTol
    !f2py intent(in) nNodesA, nTriaA, nQuadsA
    !f2py intent(in) nNodesB, nTriaB, nQuadsB
    !f2py intent(in) nNodesInt, nBarsInt
    !f2py intent(in) coorA, triaConnA, quadsConnA
    !f2py intent(in) coorB, triaConnB, quadsConnB
    !f2py intent(in) coorInt, coorIntb, barsConnInt
    !f2py intent(in) parentTriaInt
    !f2py intent(in) distTol

    ! Output variables
    real(kind=realType), dimension(3,nNodesA), intent(out) :: coorAb
    real(kind=realType), dimension(3,nNodesB), intent(out) :: coorBb
    !f2py intent(out) coorAb, coorBb

    ! Working variables
    real(kind=realType), dimension(3,2) :: BBoxA, BBoxB, BBoxAB
    logical :: overlap
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_A, innerQuadsID_A
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_B, innerQuadsID_B
    integer(kind=intType), dimension(:,:), allocatable :: allTriaConnA, allTriaConnB
    integer(kind=intType) :: nInnerTriaA, nInnerQuadsA, nInnerTriaB, nInnerQuadsB
    integer(kind=intType) :: intersect
    real(kind=realType), dimension(3) :: node1A, node2A, node3A
    real(kind=realType), dimension(3) :: node1Ab, node2Ab, node3Ab
    real(kind=realType), dimension(3) :: node1B, node2B, node3B
    real(kind=realType), dimension(3) :: node1Bb, node2Bb, node3Bb
    real(kind=realType), dimension(3) :: vecStart, vecEnd
    real(kind=realType), dimension(3) :: vecStartb, vecEndb
    real(kind=realType), dimension(3) :: vecStartTest, vecEndTest
    integer(kind=intType) :: barID, parentA, parentB

    ! EXECUTION

    ! Compute bounding boxes for each component
    call computeBBox(coorA, BBoxA)
    call computeBBox(coorB, BBoxB)

    ! Compute bounding boxes intersection
    call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

    ! We can stop if there is no bounding box intersection
    if (.not. overlap) then
      print *,'Geometry objects do not overlap.'
      return
    else
      print *,'Geometry objects overlap.'
    end if

    ! Filter elements that are inside the intersected bounding box
    call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
    innerTriaID_A, innerQuadsID_A)
    call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
    innerTriaID_B, innerQuadsID_B)

    ! Get number of inner elements
    nInnerTriaA = size(innerTriaID_A)
    nInnerQuadsA = size(innerQuadsID_A)
    nInnerTriaB = size(innerTriaID_B)
    nInnerQuadsB = size(innerQuadsID_B)

    ! Print log
    print *,'Number of interior elements in A:'
    print *,nInnerTriaA + nInnerQuadsA,'of',nTriaA + nQuadsA
    print *,'Number of interior elements in B:'
    print *,nInnerTriaB + nInnerQuadsB,'of',nTriaB + nQuadsB

    ! Split quads into triangles
    call getAllTrias(triaConnA, quadsConnA, innerTriaID_A, innerQuadsID_A, allTriaConnA)
    call getAllTrias(triaConnB, quadsConnB, innerTriaID_B, innerQuadsID_B, allTriaConnB)

    ! Update number of interior triangles
    nInnerTriaA = size(allTriaConnA,2)
    nInnerTriaB = size(allTriaConnB,2)

    ! BACKWARD PASS

    ! Initialize derivatives
    coorAb = 0.0
    coorBb = 0.0

    ! Compute intersection point derivatives
    do barID = 1,nBarsInt

      ! Get parent triangles
      parentA = parentTriaInt(1,barID)
      parentB = parentTriaInt(2,barID)

      ! Get nodal coordinates of triangle in A
      node1A = coorA(:, allTriaConnA(1,parentA))
      node2A = coorA(:, allTriaConnA(2,parentA))
      node3A = coorA(:, allTriaConnA(3,parentA))

      ! Get nodal coordinates of triangle in B
      node1B = coorB(:, allTriaConnB(1,parentB))
      node2B = coorB(:, allTriaConnB(2,parentB))
      node3B = coorB(:, allTriaConnB(3,parentB))

      ! Get bar nodes
      vecStart = coorInt(:,barsConnInt(1,barID))
      vecEnd = coorInt(:,barsConnInt(2,barID))

      ! Get bar node seeds
      vecStartb = coorIntb(:,barsConnInt(1,barID))
      vecEndb = coorIntb(:,barsConnInt(2,barID))

      ! Run the forward function and check if vecStart and vecEnd are the same
      call triTriIntersect(node1A, node2A, node3A, &
           node1B, node2B, node3B, &
           intersect, vecStartTest, vecEndTest)

      ! If vecStart and vecEnd are flipped, we also need to flip the seeds
      ! They may flip due to the FEsort algorithm
      if ((norm2(vecEnd-vecStartTest) .le. 1e-7) .or. (norm2(vecStart-vecEndTest) .le. 1e-7)) then
         ! Here we flip seeds, and we use vecStartTest as buffer for the flip operation
         vecStartTest = vecStartb
         vecStartb = vecEndb
         vecEndb = vecStartTest
      end if

      ! We always detect intersections
      intersect = 1

      ! Initialize derivatives
      node1Ab = 0.0
      node2Ab = 0.0
      node3Ab = 0.0
      node1Bb = 0.0
      node2Bb = 0.0
      node3Bb = 0.0

      call triTriIntersect_b(node1A, node1Ab, node2A, node2Ab, node3A, node3Ab, &
      node1B, node1Bb, node2B, node2Bb, node3B, node3Bb, &
      intersect, vecStart, vecStartb, vecEnd, vecEndb)

      ! Assign derivatives to coorAb
      coorAb(:, allTriaConnA(1,parentA)) = coorAb(:, allTriaConnA(1,parentA)) + node1Ab
      coorAb(:, allTriaConnA(2,parentA)) = coorAb(:, allTriaConnA(2,parentA)) + node2Ab
      coorAb(:, allTriaConnA(3,parentA)) = coorAb(:, allTriaConnA(3,parentA)) + node3Ab

      ! Assign derivatives to coorBb
      coorBb(:, allTriaConnB(1,parentB)) = coorBb(:, allTriaConnB(1,parentB)) + node1Bb
      coorBb(:, allTriaConnB(2,parentB)) = coorBb(:, allTriaConnB(2,parentB)) + node2Bb
      coorBb(:, allTriaConnB(3,parentB)) = coorBb(:, allTriaConnB(3,parentB)) + node3Bb

    end do

  end subroutine computeIntersection_b

  !=============================================

  subroutine computeIntersection_d(nNodesA, nTriaA, nQuadsA, &
    nNodesB, nTriaB, nQuadsB, &
    nNodesInt, nBarsInt, &
    coorA, coorAd, triaConnA, quadsConnA, &
    coorB, coorBd, triaConnB, quadsConnB, &
    coorInt, coorIntd, barsConnInt, &
    parentTriaInt, distTol)

    ! This function computes the backward derivatives of the intersection curve between two
    ! triangulated surface components A and B.
    ! This functions receives the seeds of the component points (coorAd, and coorBd) and returns the
    ! derivatives of the intersection points (coorIntd)
    !
    ! INPUTS
    !
    ! nNodesA: integer -> number of nodes in component A
    !
    ! nTriaA: integer -> number of triangle elements in component A
    !
    ! nQuadsA: integer -> number of quad elements in component A
    !
    ! nNodesB: integer -> number of nodes in component B
    !
    ! nTriaB: integer -> number of triangle elements in component B
    !
    ! nQuadsB: integer -> number of quad elements in component B
    !
    ! nNodesInt: integer -> number of nodes in the intersection curve
    !
    ! nBarsInt: integer -> number of bar elements in the intersection curve
    !
    ! coorA: real(3,nNodesA) -> Nodal coordinates used by component A
    !
    ! coorAd: real(3,nNodesA) -> Derivative seeds of the nodal coordinates used by component A
    !
    ! triaConnA: integer(3,nTriaA) -> Connectivities of triangle elements in component A
    !
    ! quadsConnA : integer(4,nQuadsA) -> Connectivities of quad elements in component A
    !
    ! coorB: real(3,nNodesB) -> Nodal coordinates used by component B
    !
    ! coorBd: real(3,nNodesB) -> Derivative seeds of the nodal coordinates used by component B
    !
    ! triaConnB: integer(3,nTriaB) -> Connectivities of triangle elements in component B
    !
    ! quadsConnB : integer(4,nQuadsB) -> Connectivities of quad elements in component B
    !
    ! coorInt: real(3,nNodesInt) -> Nodal coordinates used by the intersection curve
    !
    ! barsConnInt: integer(2,nBarsInt) -> Connectivities of bar elements in the intersection curve
    !
    ! parentTriaInt: integer(2,nBarsInt) -> Triangle pairs (from component A and component B) that
    !                                       generated each intersection bar element.
    !
    ! distTol: real -> distance used to check if the start and end nodes of the intersection are flipped.
    !                  use the same distTol used for the original computeIntersection.
    !
    ! OUTPUTS
    !
    ! coorIntd: real(3,nNodesInt) -> Derivative seed of the nodal coordinates used by the intersection curve
    !
    !
    ! Ney Secco 2016-09

    use Intersection
    use Intersection_d, only: triTriIntersect_d
    use Utilities ! This will bring condenseBarNodes_main
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: nNodesA, nTriaA, nQuadsA
    integer(kind=intType), intent(in) :: nNodesB, nTriaB, nQuadsB
    integer(kind=intType), intent(in) :: nNodesInt, nBarsInt
    real(kind=realType), dimension(3,nNodesA), intent(in) :: coorA
    real(kind=realType), dimension(3,nNodesA), intent(in) :: coorAd
    integer(kind=intType), dimension(3,nTriaA), intent(in) :: triaConnA
    integer(kind=intType), dimension(4,nQuadsA), intent(in) :: quadsConnA
    real(kind=realType), dimension(3,nNodesB), intent(in) :: coorB
    real(kind=realType), dimension(3,nNodesB), intent(in) :: coorBd
    integer(kind=intType), dimension(3,nTriaB), intent(in) :: triaConnB
    integer(kind=intType), dimension(4,nQuadsB), intent(in) :: quadsConnB
    real(kind=realType), dimension(3,nNodesInt), intent(in) :: coorInt ! Intersection nodes
    integer(kind=intType), dimension(2,nBarsInt), intent(in) :: barsConnInt ! Intersection bar connectivities
    integer(kind=intType), dimension(2,nBarsInt), intent(in) :: parentTriaInt ! Pair of triangles from A and B that generates each bar
    real(kind=realType), intent(in) :: distTol
    !f2py intent(in) nNodesA, nTriaA, nQuadsA
    !f2py intent(in) nNodesB, nTriaB, nQuadsB
    !f2py intent(in) nNodesInt, nBarsInt
    !f2py intent(in) coorA, triaConnA, quadsConnA
    !f2py intent(in) coorB, triaConnB, quadsConnB
    !f2py intent(in) coorAd, coorBd
    !f2py intent(in) coorInt, barsConnInt
    !f2py intent(in) parentTriaInt
    !f2py intent(in) distTol

    ! Output variables
    real(kind=realType), dimension(3,nNodesInt), intent(out) :: coorIntd ! Intersection nodes
    !f2py intent(out) coorIntd

    ! Working variables
    real(kind=realType), dimension(3,2) :: BBoxA, BBoxB, BBoxAB
    logical :: overlap
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_A, innerQuadsID_A
    integer(kind=intType), dimension(:), allocatable :: innerTriaID_B, innerQuadsID_B
    integer(kind=intType), dimension(:,:), allocatable :: allTriaConnA, allTriaConnB
    integer(kind=intType) :: nInnerTriaA, nInnerQuadsA, nInnerTriaB, nInnerQuadsB
    integer(kind=intType) :: intersect
    real(kind=realType), dimension(3) :: node1A, node2A, node3A
    real(kind=realType), dimension(3) :: node1Ad, node2Ad, node3Ad
    real(kind=realType), dimension(3) :: node1B, node2B, node3B
    real(kind=realType), dimension(3) :: node1Bd, node2Bd, node3Bd
    real(kind=realType), dimension(3) :: vecStart, vecEnd
    real(kind=realType), dimension(3) :: vecStartd, vecEndd
    real(kind=realType), dimension(3) :: vecStartTest, vecEndTest
    integer(kind=intType) :: barID, parentA, parentB

    ! EXECUTION

    ! Compute bounding boxes for each component
    call computeBBox(coorA, BBoxA)
    call computeBBox(coorB, BBoxB)

    ! Compute bounding boxes intersection
    call computeBBoxIntersection(BBoxA, BBoxB, BBoxAB, overlap)

    ! We can stop if there is no bounding box intersection
    if (.not. overlap) then
      print *,'Geometry objects do not overlap.'
      return
    else
      print *,'Geometry objects overlap.'
    end if

    ! Filter elements that are inside the intersected bounding box
    call filterElements(coorA, triaConnA, quadsConnA, BBoxAB, &
    innerTriaID_A, innerQuadsID_A)
    call filterElements(coorB, triaConnB, quadsConnB, BBoxAB, &
    innerTriaID_B, innerQuadsID_B)

    ! Get number of inner elements
    nInnerTriaA = size(innerTriaID_A)
    nInnerQuadsA = size(innerQuadsID_A)
    nInnerTriaB = size(innerTriaID_B)
    nInnerQuadsB = size(innerQuadsID_B)

    ! Print log
    print *,'Number of interior elements in A:'
    print *,nInnerTriaA + nInnerQuadsA,'of',nTriaA + nQuadsA
    print *,'Number of interior elements in B:'
    print *,nInnerTriaB + nInnerQuadsB,'of',nTriaB + nQuadsB

    ! Split quads into triangles
    call getAllTrias(triaConnA, quadsConnA, innerTriaID_A, innerQuadsID_A, allTriaConnA)
    call getAllTrias(triaConnB, quadsConnB, innerTriaID_B, innerQuadsID_B, allTriaConnB)

    ! Update number of interior triangles
    nInnerTriaA = size(allTriaConnA,2)
    nInnerTriaB = size(allTriaConnB,2)

    ! BACKWARD PASS

    ! Initialize derivatives
    coorIntd = 0.0

    ! Compute intersection point derivatives
    do barID = 1,nBarsInt

      ! Get parent triangles
      parentA = parentTriaInt(1,barID)
      parentB = parentTriaInt(2,barID)

      ! Get nodal coordinates of triangle in A
      node1A = coorA(:, allTriaConnA(1,parentA))
      node2A = coorA(:, allTriaConnA(2,parentA))
      node3A = coorA(:, allTriaConnA(3,parentA))

      ! Get nodal coordinates of triangle in B
      node1B = coorB(:, allTriaConnB(1,parentB))
      node2B = coorB(:, allTriaConnB(2,parentB))
      node3B = coorB(:, allTriaConnB(3,parentB))

      ! Get bar nodes
      vecStart = coorInt(:,barsConnInt(1,barID))
      vecEnd = coorInt(:,barsConnInt(2,barID))

      ! We always detect intersections
      intersect = 1

      ! Assign derivatives to coorAb
      node1Ad = coorAd(:, allTriaConnA(1,parentA))
      node2Ad = coorAd(:, allTriaConnA(2,parentA))
      node3Ad = coorAd(:, allTriaConnA(3,parentA))

      ! Assign derivatives to coorBb
      node1Bd = coorBd(:, allTriaConnB(1,parentB))
      node2Bd = coorBd(:, allTriaConnB(2,parentB))
      node3Bd = coorBd(:, allTriaConnB(3,parentB))

      ! Initialize derivatives
      vecStartd = 0.0
      vecEndd = 0.0

      call triTriIntersect_d(node1A, node1Ad, node2A, node2Ad, node3A, node3Ad, &
      node1B, node1Bd, node2B, node2Bd, node3B, node3Bd, &
      intersect, vecStartTest, vecStartd, vecEndTest, vecEndd)

      ! If vecStart and vecEnd are flipped, we also need to flip the seeds.
      ! They may flip due to the FEsort algorithm
      if ((norm2(vecEnd-vecStartTest) .le. 1e-7) .or. (norm2(vecStart-vecEndTest) .le. 1e-7)) then
         ! Here we flip seeds, and we use vecStartTest as buffer for the flip operation
         vecStartTest = vecStartd
         vecStartd = vecEndd
         vecEndd = vecStartTest
      end if

      ! Get bar node seeds
      coorIntd(:,barsConnInt(1,barID)) = coorIntd(:,barsConnInt(1,barID)) + vecStartd
      coorIntd(:,barsConnInt(2,barID)) = coorIntd(:,barsConnInt(2,barID)) + vecEndd

    end do

  end subroutine computeIntersection_d

  !=============================================

  subroutine releaseMemory()

    ! This subroutine just deallocates memory used by the intersection code.
    ! Remember to call this function after you copy the outputs in Python.

    implicit none

    ! Deallocate output variables
    if (allocated(coor)) then
      deallocate(coor, barsConn, parentTria)
    end if

  end subroutine releaseMemory

  !=============================================

  subroutine testTri(V0, V1, V2, U0, U1, U2, intersect, vecStart, vecEnd)

    ! This function computes the intersection curve between two
    ! triangles.
    !
    ! John Jasa 2016-08

    use Intersection
    implicit none

    real(kind=realType), dimension(3), intent(in) :: V0, V1, V2, U0, U1, U2
    integer(kind=intType), intent(out) :: intersect
    real(kind=realType), dimension(3), intent(out) :: vecStart, vecEnd

    call triTriIntersect(V0, V1, V2, U0, U1, U2, intersect, vecStart, vecEnd)

  end subroutine testTri

end module intersectionAPI
