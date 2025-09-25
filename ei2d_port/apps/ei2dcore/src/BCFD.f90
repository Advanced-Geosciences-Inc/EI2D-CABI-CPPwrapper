
!  This subroutine implements Dirichlet and mixed boundary conditions
!  for finite difference method. Note that stiffness is multiplied 
!  by a factor 2. This is made up after this routine.
!  Ref: Dey and Morrison (1979), Geophysical Prospecting, 27, 106-136
!  
!                                       ------------
!  iNode is the center node with       |     |      |
!  four surrounding elements           |  n1 |  n3  |  dy1
!  n1, n2, n3, n4                      |     |      |
!                                      --------------
!                                      |     |      | 
!                                      |  n2 |  n4  |  dy2
!                                      |     |      |
!                                       ------------
!                                        dx1   dx2

  subroutine BCFD(StiffMatrix,NodeX,NodeY,Cond,nBand,iAbscissa)

  use GlobalForw
  implicit none

! variables and arrays passed from Delphi GUI

  integer, intent(in) :: nBand
  real(Rkind), intent(inout), dimension(1:nBand*gNumNodes) :: StiffMatrix
  real(Rkind), intent(in), dimension(1:gNumNodes)  :: NodeX, NodeY
  real(Rkind), intent(in), dimension(1:gNumElem)   :: Cond
  integer, intent(in) :: iAbscissa

! variables and arrays declared in this subroutine

  integer :: i, k, iPos, iNode, n1, n2, n3, n4, ix, iy
  real(Rkind) :: Cl, Ct, Cr, Cb, Cp, abscissa2
  real(Rkind) :: dx1, dx2, dy1, dy2


  abscissa2 = gfAbscissa(iAbscissa) * gfAbscissa(iAbscissa) 
!----------------------------------------------------------------------
!  There are total 8 different cases of boundary condition: four corner
!  nodes and four side boudaries.
!
! (1) Top surface boundary condition is always of the Neumann Type.
!      This is mandatory B.C.

  do i = 2, gNumNodeX-1
     iNode = (i-1) * gNumNodeY + 1
     n2 = (i-2) * gNumElemY + 1
     n4 = n2 + gNumElemY
     dx1 = abs(NodeX(iNode) - NodeX(iNode-gNumNodeY))
     dx2 = abs(NodeX(iNode+gNumNodeY) - NodeX(iNode))
     dy2 = abs(NodeY(iNode+1) - NodeY(iNode))

     Cl = -dy2 * Cond(n2) / dx1
     Cr = -dy2 * Cond(n4) / dx2
     Cb = -(dx1*Cond(n2) + dx2*Cond(n4))/dy2
     Cp = 0.5D0*abscissa2*(Cond(n4)*dx2*dy2 + Cond(n2)*dx1*dy2) - Cl - Cr - Cb 

     ! Cholesky decomposition. Stores banded lower triangle.
     if (gnForwSolver == 0) then             
        iPos = (iNode-1)*nBand + 1
        StiffMatrix(iPos) = Cl
        StiffMatrix(iNode*nBand) = Cp

     ! Conjugate gradient. Stores upper triangle compactly 
     else if (gnForwSolver == 1) then      
        iPos = (iNode-1) * nBand + 1
        StiffMatrix(iPos) = Cp
        iPos = (iNode-gNumNodeY) * nBand
        StiffMatrix(iPos) = Cl
     end if
  end do

!----------------------------------------------------------------------
! (2) Top left corner node on the surface (node #1)

  dx2 = abs(NodeX(1+gNumNodeY) - NodeX(1))
  dy2 = abs(NodeY(2) - NodeY(1))
  Cb = -dx2 * Cond(1) / dy2
  Cr = -dy2 * Cond(1) / dx2
  Cp = 0.5D0*abscissa2*Cond(1)*dx2*dy2 - Cb - Cr
  if (gnBCType == 1) then                             ! mixed BC
     Cp = Cp + dy2*Cond(1)*gfEtaLeft(1, iAbscissa)
  end if

  ! Cholesky decomposition. Stores banded lower triangle.
  if (gnForwSolver == 0) then             
     StiffMatrix(nBand) = Cp

  ! Conjugate gradient. Stores upper triangle compactly 
  else if (gnForwSolver == 1) then      
     StiffMatrix(1) = Cp
  end if
!----------------------------------------------------------------------
! (3) Top right corner node on the surface

  iNode = (gNumNodeX-1) * gNumNodeY + 1
  n2 = (gNumElemX-1) * gNumElemY + 1

  dx1 = abs(NodeX(iNode) - NodeX(iNode-gNumNodeY))
  dy2 = abs(NodeY(iNode+1) - NodeY(iNode))
  Cl = -dy2 * Cond(n2) / dx1
  Cb = -dx1*Cond(n2)/dy2
  Cp = 0.50D0*abscissa2*Cond(n2)*dx1*dy2 - Cb - Cl

  if (gnBCType == 1) then    ! mixed BC
     Cp = Cp + dy2*Cond(n2)*gfEtaRight(1, iAbscissa)
  end if

  ! Cholesky decomposition. Stores banded lower triangle.
  if (gnForwSolver == 0) then             
     iPos = (iNode-1)*nBand + 1
     StiffMatrix(iPos) = Cl
     StiffMatrix(iNode*nBand) = Cp

  ! Conjugate gradient. Stores upper triangle compactly 
  else if (gnForwSolver == 1) then      
     iPos = (iNode-1) * nBand + 1
     StiffMatrix(iPos) = Cp
     iPos = (iNode-gNumNodeY) * nBand
     StiffMatrix(iPos) = Cl
  end if
!----------------------------------------------------------------------
! (4) Nodes located on the bottom edge

  do i = 2, gNumNodeX-1
     iNode = i * gNumNodeY
     n1 = (i-1) * gNumElemY
     n3 = i * gNumElemY
     dx1 = abs(NodeX(iNode) - NodeX(iNode-gNumNodeY))
     dx2 = abs(NodeX(iNode+gNumNodeY) - NodeX(iNode))
     dy1 = abs(NodeY(iNode) - NodeY(iNode-1))

     Cl = -dy1 * Cond(n1)/ dx1
     Cr = -dy1 * Cond(n3)/ dx2
     Ct = -(dx1 * Cond(n1) + dx2 * Cond(n3)) / dy1
     Cp = 0.50D0*abscissa2*(Cond(n1)*dx1*dy1 + Cond(n3)*dx2*dy1)
     Cp = Cp - Cl - Cr - Ct

     if (gnBCType == 1) then     ! mixed BC
        Cp = Cp + (dy1*Cond(n1) + dx2*Cond(n3))*gfEtaBottom(i,iAbscissa)
     end if

     ! Cholesky decomposition. Stores banded lower triangle.
     if (gnForwSolver == 0) then             
        iPos = (iNode-1)*nBand + 1
        StiffMatrix(iPos) = Cl
        StiffMatrix(iPos + gNumNodeY - 1) = Ct
        StiffMatrix(iNode*nBand) = Cp

     ! Conjugate gradient. Stores upper triangle compactly 
     else if (gnForwSolver == 1) then      
        iPos = (iNode-1) * nBand + 1
        StiffMatrix(iPos) = Cp
        iPos = (iNode-2) * nBand + 2
        StiffMatrix(iPos) = Ct
        iPos = (iNode-gNumNodeY) * nBand
        StiffMatrix(iPos) = Cl
     end if
  end do  
!----------------------------------------------------------------------
! (5) Bottom left corner node

  dx2 = abs(NodeX(2*gNumNodeY) - NodeX(gNumNodeY))
  dy1 = abs(NodeY(gNumNodeY-1) - NodeY(gNumNodeY))
  Cr = -dy1 * Cond(gNumElemY) / dx2
  Ct = -dx2*Cond(gNumElemY)/dy1
  Cp = 0.50D0*abscissa2*Cond(gNumElemY)*dx2*dy1 - Ct - Cr

  if (gnBCType == 1) then   ! mixed BC
     Cp = Cp + dx2*Cond(gNumElemY)*gfEtaBottom(1, iAbscissa) +  &
              dy1*Cond(gNumElemY)*gfEtaLeft(gNumNodeY, iAbscissa)
  end if

  ! Cholesky decomposition. Stores banded lower triangle.
  if (gnForwSolver == 0) then             
     iPos = (gNumNodeY-1)*nBand + gNumNodeY
     StiffMatrix(iPos)   = Ct
     StiffMatrix(iPos+1) = Cp

  ! Conjugate gradient. Stores upper triangle compactly 
  else if (gnForwSolver == 1) then      
     iPos = (gNumNodeY-1) * nBand + 1
     StiffMatrix(iPos) = Cp
     iPos = (iNode-2) * nBand + 2
     StiffMatrix(iPos) = Ct
  end if
!----------------------------------------------------------------------
! (6) Bottom right corner node - last node = gNumNode

  dx1 = abs(NodeX(gNumNodes) - NodeX(gNumNodes-gNumNodeY))
  dy1 = abs(NodeY(gNumNodes) - NodeY(gNumNodes-1))
  Cl = -dy1 * Cond(gNumElem) / dx1
  Ct = -dx1*Cond(gNumElem)/dy1
  Cp = 0.50D0*abscissa2*Cond(gNumElem)*dx1*dy1 - Ct - Cl
 
  if (gnBCType == 1) then   ! mixed BC      
     Cp = Cp + dx1*Cond(gNumElem)*gfEtaBottom(gNumNodeX, iAbscissa) +  &
            dy1*Cond(gNumElem)*gfEtaRight(gNumNodeY, iAbscissa)
  end if

  ! Cholesky decomposition. Stores banded lower triangle.
  if (gnForwSolver == 0) then             
     iPos = (gNumNodes-1)*nBand + 1
     StiffMatrix(iPos) = Cl
     StiffMatrix(iPos + gNumNodeY - 1) = Ct
     StiffMatrix(gNumNodes*nBand) = Cp

  ! Conjugate gradient. Stores upper triangle compactly 
  else if (gnForwSolver == 1) then      
     iPos = (gNumNodes-1) * nBand + 1
     StiffMatrix(iPos) = Cp
     iPos = (gNumNodes-2) * nBand + 2
     StiffMatrix(iPos) = Ct
     iPos = (gNumNodes-gNumNodeY) * nBand
     StiffMatrix(iPos) = Cl
  end if

!----------------------------------------------------------------------
! (7) Nodes on the left edge

  do iNode = 2, gNumNodeY-1
     dx2 = abs(NodeX(iNode+gNumNodeY) - NodeX(iNode))
     dy1 = abs(NodeY(iNode) - NodeY(iNode-1))
     dy2 = abs(NodeY(iNode+1) - NodeY(iNode))

     Cr = -(dy2*Cond(iNode) + dy1*Cond(iNode-1)) / dx2
     Ct = -dx2 * Cond(iNode-1) / dy1
     Cb = -dx2 * Cond(iNode) / dy2
     Cp = 0.50D0*abscissa2*(Cond(iNode-1)*dx2*dy1 + Cond(iNode)*dx2*dy2)
     Cp = Cp - Cb - Cr - Ct

     if (gnBCType == 1) then   ! mixed BC
        Cp = Cp + (dy2*Cond(iNode) + dy1*Cond(iNode-1))*gfEtaLeft(iNode,iAbscissa)
     end if

     ! Cholesky decomposition. Stores banded lower triangle.
     if (gnForwSolver == 0) then             
        iPos = (iNode-1)*nBand + gNumNodeY
        StiffMatrix(iPos) = Ct
        StiffMatrix(iNode*nBand) = Cp

     ! Conjugate gradient. Stores upper triangle compactly 
     else if (gnForwSolver == 1) then      
        iPos = (iNode-1) * nBand + 1
        StiffMatrix(iPos) = Cp
        iPos = (iNode-2) * nBand + 2
        StiffMatrix(iPos) = Ct
     end if
  end do  
!----------------------------------------------------------------------
! (8) Nodes on the right edge

  do i = 2, gNumNodeY-1
     iNode = (gNumNodeX-1)*gNumNodeY + i
     n1 = (gNumElemX-1) * gNumElemY + i - 1
     n2 = n1 + 1
     dx1 = abs(NodeX(iNode-gNumNodeY) - NodeX(iNode))
     dy1 = abs(NodeY(iNode) - NodeY(iNode-1))
     dy2 = abs(NodeY(iNode+1) - NodeY(iNode))

     Cl = -(dy2*Cond(n2) + dy1*Cond(n1)) / dx1
     Ct = -dx1 * Cond(n1) / dy1
     Cb = -dx1 * Cond(n2) / dy2
     Cp = 0.50D0*abscissa2*(Cond(n1)*dx1*dy1 + Cond(n2)*dx1*dy2)
     Cp = Cp - Cb - Cl - Ct

     if (gnBCType == 1) then  ! mixed BC
        Cp = Cp + (dy2*Cond(n2) + dy1*Cond(n1))*gfEtaRight(i,iAbscissa)
     end if

     ! Cholesky decomposition. Stores banded lower triangle.
     if (gnForwSolver == 0) then             
        iPos = (iNode-1)*nBand + 1
        StiffMatrix(iPos) = Cl
        StiffMatrix(iPos + gNumNodeY - 1) = Ct
        StiffMatrix(iNode*nBand) = Cp

     ! Conjugate gradient. Stores upper triangle compactly 
     else if (gnForwSolver == 1) then      
        iPos = (iNode-1) * nBand + 1
        StiffMatrix(iPos) = Cp
        iPos = (iNode-2) * nBand + 2
        StiffMatrix(iPos) = Ct
        iPos = (iNode-gNumNodeY) * nBand
        StiffMatrix(iPos) = Cl
     end if
  end do  
!---------------------------------------------------------------------------
!  Add symmetric coupling coefficients to diagonal terms of interior nodes

  do ix = 2, gNumNodeX-1            
     do iy = 2, gNumNodeY-1

        iNode = (ix-1) * gNumNodeY + iy

        ! Cholesky decomposition. Stores banded lower triangle.
        if (gnForwSolver == 0) then             
           iPos = iNode*nBand
           n1 = (iNode+gNumNodeY-1)*nBand + 1
           n2 = (iNode+1) * nBand - 1

        ! Conjugate gradient. Stores upper triangle compactly 
        else if (gnForwSolver == 1) then      
           iPos = (iNode-1) * nBand + 1
           n1 = iNode * nBand
           n2 = iPos + 1
        end if

        StiffMatrix(iPos) = StiffMatrix(iPos) - StiffMatrix(n1) - StiffMatrix(n2)

     end do
  end do
!----------------------------------------------------------------------
  if (gnBCType == 0) then                        ! Dirichlet 
!----------------------------------------------------------------------
!  Implement Dirichlet boundary conditions. First, both side boundaries
!  The method we used can be found in pages 46-53 of Huebner, Thornton 
!  and Byrom (1995). We set diagonal terms of stiffness matrix at 
!  boundary nodes to a large number 1.0d+15. This will force the potential 
!  to be around 0. For single precision processing, 1.0d+10 may be used
!  to replace 1.0d+15.
!----------------------------------------------------------------------
     k = gNumNodes - gNumNodeY    
     do i = 1, gNumNodeY               
        ! Left boundary
        if (gnForwSolver == 0) then             !  CD
           iPos  = i * nBand
        else                                    !  CG
           iPos = (i-1)*nBand + 1
        end if
        StiffMatrix(iPos)  = StiffMatrix(iPos)  * 1.0d+15
        ! Right Boundary
        if (gnForwSolver == 0) then
           iPos = (i+k) * nBand
        else 
           iPos = (i+k-1)*nBand + 1
        end if
        StiffMatrix(iPos) = StiffMatrix(iPos) * 1.0d+15
     end do 

     !  Now, set bottom Dirichlet boundary conditions
     do i = 2, gNumNodeX-1        
        if (gnForwSolver == 0) then              ! CD
           iPos = i * gNumNodeY * nBand
        else                                     ! CG
           iPos = (i*gNumNodeY-1)*nBand + 1
        end if
        StiffMatrix(iPos) = StiffMatrix(iPos) * 1.0d+15
     end do

  end if
!----------------------------------------------------------------------
  return
  end subroutine BCFD
