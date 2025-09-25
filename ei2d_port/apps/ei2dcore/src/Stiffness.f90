
!  StiffnessFE: Stiffness of a quadrilateral element of finite element method
!  StiffnessFD: Stiffness of an internal node of finite difference method
!
!  Construct a 5 x 5 local stiffness matrix of a single element. Then
!  reduce the 5 x 5 ElemStiff into 4 x 4 ReducedStiff using static 
!  condensation. See IP2DI report of Wannamaker (1992), a book in
!  Chinese by Luo and Zhang (1987), and Luiz Rijo's PhD thesis (1977).
! 
!  Though the complete conductivity array is passed, I use the conductivity
!  within a single element only. This subroutine may be modified to deal
!  with triangular element.
!-------------------------------------------------------------------------------
  subroutine StiffnessFE(NodeX, NodeY, Conductivity, ReducedStiff, abscix,     &
                         iElemX, iElemY, LocalNodes, CenterNodeX, CenterNodeY, &
                         ElemArea)

  use GlobalForw
  implicit none

  integer, intent(in) :: iElemX, iElemY
  real(Rkind), intent(in), dimension(1:gNumNodes) :: NodeX, NodeY
  real(Rkind), intent(in), dimension(1:gNumElem) :: Conductivity
  real(Rkind), intent(in) :: abscix
  integer, intent(out) :: LocalNodes(1:4)
  real(Rkind), intent(in), dimension(1:gNumElem)   :: CenterNodeX, CenterNodeY
  real(Rkind), intent(in), dimension(1:gNumElem*4) :: ElemArea
  real(Rkind), intent(out), dimension(1:4,1:4) :: ReducedStiff

  real(Rkind) :: tmp, tmp1, ElemStiff(1:5,1:5),b(1:3),c(1:3)
  integer     :: iElem, i, j, k, m, n1, n2

! Nodes are globally numbered from top to bottom and from left to right.
! But the LocalNodes of the current element iElem are numbered in a way
! illustrated in the figure below. Nodes 3 and 4 are numbered oppsite to
! the global numbering system. This will be reverted at the end of this 
! subroutine.
! 
!                                 1             4 
!                                   ------------
!                                  |            | 
!                                  |            |
!                                  |            |
!                                   ------------
!                                 2             3
!
  iElem = iElemY + gNumElemY * (iElemX - 1)
  LocalNodes(1) = (iElemX-1) * gNumNodeY + iElemY
  LocalNodes(2) = LocalNodes(1) + 1
  LocalNodes(3) = LocalNodes(2) + gNumNodeY
  LocalNodes(4) = LocalNodes(1) + gNumNodeY

!---------------------------------------------------------------------------
! Assemble local stiffness matrix of a quadrilateral element from four
! triangular elements. Note that ElemStiff is symmetric. The fifth node (k) 
! can be "removed" from global StiffMatrix by static condensation below
!---------------------------------------------------------------------------
  k = 5    
  ElemStiff = 0.0D0

  do i = 1, 4              ! Four triangular elements

     j = mod(i, 4) + 1     ! To provent j = k = 5
     n1 = LocalNodes(i)
     n2 = LocalNodes(j)    ! n3 is the center node

     b(1) = NodeY(n2) - CenterNodeY(iElem) 
     b(2) = CenterNodeY(iElem) - NodeY(n1) 
     b(3) = NodeY(n1) - NodeY(n2) 
     c(1) = CenterNodeX(iElem) - NodeX(n2) 
     c(2) = NodeX(n1) - CenterNodeX(iElem) 
     c(3) = NodeX(n2) - NodeX(n1) 

     ! tmp is a term introduced by Fourier transform
     ! Multiplication by element conductivity is put off to the end
     ! Note that conductivity is not applied in this do loop
     ! There could be a sign flip below due to the conflicts in Luo (1987)
     ! and Wannamaker (1992)
     tmp = ElemArea((iElem-1)*4+i)  
     tmp1 = 0.25D0 / tmp
     tmp  = tmp * abscix

     ElemStiff(i,i) = ElemStiff(i,i) + tmp1 * (b(1)*b(1)+c(1)*c(1)) + tmp*2.
     ElemStiff(i,j) = ElemStiff(i,j) + tmp1 * (b(2)*b(1)+c(2)*c(1)) + tmp
     ElemStiff(i,k) = ElemStiff(i,k) + tmp1 * (b(3)*b(1)+c(3)*c(1)) + tmp
     ElemStiff(j,j) = ElemStiff(j,j) + tmp1 * (b(2)*b(2)+c(2)*c(2)) + tmp*2.
     ElemStiff(j,k) = ElemStiff(j,k) + tmp1 * (b(3)*b(2)+c(3)*c(2)) + tmp
     ElemStiff(k,k) = ElemStiff(k,k) + tmp1 * (b(3)*b(3)+c(3)*c(3)) + tmp*2.

     ElemStiff(j,i) = ElemStiff(i,j)   ! Symmetric
     ElemStiff(k,i) = ElemStiff(i,k)
     ElemStiff(k,j) = ElemStiff(j,k)

  end do

! Switch Nodes 3 and 4, so global and local node numbering agree with each other
  LocalNodes(3) = LocalNodes(1) + gNumNodeY  
  LocalNodes(4) = LocalNodes(2) + gNumNodeY
!------------------------------------------------------------------------------------
!  Static Condensation. See Wannamaker (1992), Luo and Zhang (1987), 
!  Huebner, Thornton and Byrom (1995). Note that ReducedStiff is symmetric.
!  This transform reduces total number of nodes by almost 50%

  ReducedStiff = 0.0D0
  tmp = 1.0/ElemStiff(5,5)
  do i = 1, 4
     k = i
     if (i==3) k = 4        ! switch local nodes 3 and 4
     if (i==4) k = 3
     do j = i, 4            ! Symmetric stiffness
        m = j
        if (j==3) m = 4
        if (j==4) m = 3    
        ReducedStiff(i,j) = ElemStiff(k,m) - ElemStiff(k,5)*ElemStiff(5,m)*tmp  
     end do
  end do

  do i = 2, 4               ! Symmetric stiffness 
     do j = 1, i-1
        ReducedStiff(i,j) = ReducedStiff(j,i)
     end do
  end do

! ???: Assumption: Four triangular elements have the identical conductivity
! This may not be true for surveys with topography
  ReducedStiff = ReducedStiff * Conductivity(iElem)  

  return
  end subroutine StiffnessFE




!**************************************************************************
!==========================================================================
!
! subroutine StiffnessFD constructs stiffness of an internal node.
! Note that only Self, Top and Left coupling coefficients are calculated.
! Bottom and Right components of Self coupling will be added later in BCFD.f90
! according to symmetry of the stiffness matrix.
! Ref: Dey and Morrison (1979), Geophysical Prospecting, 27, 106-136
!
!  
!  Discretization by area               ------------
!  iNode is the center node with       |     |      |
!  four surrounding elements           |  n1 |  n3  |  dy1
!  n1, n2, n3, n4                      |     |      |
!                                      --------------
!                                      |     |      | 
!                                      |  n2 |  n4  |  dy2
!                                      |     |      |
!                                       ------------
!                                        dx1   dx2
!
!----------------------------------------------------------------------------------
  subroutine StiffnessFD(iNodeX, iNodeY, NodeX, NodeY, Cond, abscissa, Cl, Ct, Cp)

  use GlobalForw
  implicit none

  integer, intent(in) :: iNodeX, iNodeY
  real(Rkind), intent(in), dimension(1:gNumNodes) :: NodeX, NodeY
  real(Rkind), intent(in), dimension(1:gNumElem) :: Cond
  real(Rkind), intent(in) :: abscissa
  real(Rkind), intent(out) :: Cl, Ct, Cp

  integer :: iNode, n1, n2, n3, n4
  real(Rkind) :: dx1, dx2, dy1, dy2
!-----------------------------------------------------------------

  iNode = (iNodeX-1) * gNumNodeY + iNodeY

  ! four cells (elements) surrounding iNode
  n1 = (iNodeX-2) * gNumElemY + iNodeY - 1
  n2 = n1 + 1
  n3 = n1 + gNumElemY
  n4 = n3 + 1
    
  dx1 = abs(NodeX(iNode) - NodeX(iNode-gNumNodeY))
  dx2 = abs(NodeX(iNode+gNumNodeY) - NodeX(iNode))
  dy1 = abs(NodeY(iNode) - NodeY(iNode-1))
  dy2 = abs(NodeY(iNode+1) - NodeY(iNode))
  

  Cl = -(dy1 * Cond(n1) + dy2 * Cond(n2)) / dx1
  Ct = -(dx1 * Cond(n1) + dx2 * Cond(n3)) / dy1
  Cp = Cond(n1)*dx1*dy1 + Cond(n3)*dx2*dy1 + Cond(n4)*dx2*dy2 + Cond(n2)*dx1*dy2
  Cp = 0.5D0*abscissa*abscissa*Cp - Cl - Ct

  return
  end subroutine StiffnessFD
