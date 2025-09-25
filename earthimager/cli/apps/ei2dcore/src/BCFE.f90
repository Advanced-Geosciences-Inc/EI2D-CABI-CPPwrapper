
!  This subroutine implements Dirichlet or mixed boundary condition
!  for finite element method

  subroutine BCFE(StiffMatrix,NodeX,NodeY,Cond,nBand,iAbscissa)

  use GlobalForw
  implicit none

! variables and arrays passed from Delphi GUI

  integer, intent(in) :: nBand  
  real(Rkind), intent(inout), dimension(1:nBand*gNumNodes) :: StiffMatrix
  real(Rkind), intent(in), dimension(1:gNumNodes)  :: NodeX, NodeY
  real(Rkind), intent(in), dimension(1:gNumElem)   :: Cond
  integer, intent(in) :: iAbscissa

! variables and arrays declared in this subroutine

  integer :: i, j, k, n, iPos
  real(Rkind) :: tmp
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
        else 
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
        else 
           iPos = (i*gNumNodeY-1)*nBand + 1
        end if
        StiffMatrix(iPos) = StiffMatrix(iPos) * 1.0d+15
     end do
!----------------------------------------------------------------------
  else if (gnBCType == 1) then           ! mixed boundary condition 
!----------------------------------------------------------------------
!  The only reference of implementation below is Luo and Zhang (1987).
!  I scaled the element contribution by a factor 1/2 due to its 
!  discrepancy of stiffness from Wannamaker's algorithm
!----------------------------------------------------------------------
     do i = 1, gNumNodeY - 1                    ! Left boundary            
        if (gnForwSolver == 0) then             ! CD
           iPos  = i * nBand
           j = iPos + nBand
           k = j - 1
        else                                    ! CG
           iPos = (i-1)*nBand + 1
           j = iPos + nBand
           k = iPos + 1
        end if
        tmp = abs(NodeY(i) - NodeY(i+1)) * Cond(i) / 3.0D0
        StiffMatrix(iPos) = StiffMatrix(iPos) + gfEtaLeft(i,iAbscissa) * tmp
        StiffMatrix(j)    = StiffMatrix(j) + gfEtaLeft(i+1,iAbscissa)  * tmp
        StiffMatrix(k)    = StiffMatrix(k) + gfEtaLeft(i+1,iAbscissa)  * tmp * 0.50D0
     end do

     n = gNumNodes - gNumNodeY    
     do i = 1, gNumNodeY - 1                    ! Right boundary            
        if (gnForwSolver == 0) then             ! CD
           iPos = (i+n) * nBand
           j = iPos + nBand
           k = j - 1
        else                                    ! CG
           iPos = (i+n-1)*nBand + 1
           j = iPos + nBand
           k = iPos + 1
        end if
        tmp = abs(NodeY(n+i) - NodeY(n+i+1)) * Cond(gNumElemY*(gNumElemX-1)+i) / 3.0D0
        StiffMatrix(iPos) = StiffMatrix(iPos) + gfEtaRight(i,iAbscissa) * tmp
        StiffMatrix(j)    = StiffMatrix(j) + gfEtaRight(i+1,iAbscissa)  * tmp
        StiffMatrix(k)    = StiffMatrix(k) + gfEtaRight(i+1,iAbscissa)  * tmp * 0.50D0
     end do

     do i = 1, gNumNodeX-1                       ! bottom boundary 
        if (gnForwSolver == 0) then              ! CD
           iPos = i * gNumNodeY * nBand
           j = (i+1) * gNumNodeY * nBand
           k = j - gNumNodeY
        else                                     ! CG
           iPos = (i*gNumNodeY-1)*nBand + 1
           j = ((i+1)*gNumNodeY-1)*nBand + 1
           k = iPos + 3
        end if
        tmp = abs(NodeX(i*gNumNodeY) - NodeX((i+1)*gNumNodeY)) * Cond(i*gNumElemY) / 3.0D0
        StiffMatrix(iPos) = StiffMatrix(iPos) + gfEtaBottom(i,iAbscissa) * tmp
        StiffMatrix(j)    = StiffMatrix(j) + gfEtaBottom(i+1,iAbscissa)  * tmp
        StiffMatrix(k)    = StiffMatrix(k) + gfEtaBottom(i+1,iAbscissa)  * tmp * 0.50D0
     end do
  end if                               ! end of boundary conditions

!-----------------------------------------------------------------------------------------
  return
  end subroutine BCFE
