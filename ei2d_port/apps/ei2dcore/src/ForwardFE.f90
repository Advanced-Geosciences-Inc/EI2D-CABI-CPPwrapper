
!  This is a driver routine for forward modeling using finite element method.
!  We follow an algorithm described in Wannamaker (1992), Luo and Zhang (1987), 
!  and Luiz Rijo's PhD thesis (1977).

!  gnForwSolver = 0, Cholesky decomposition;  = 1, conjugate gradient.

!  There is a factor 2 discrepancy between Luo and Wannamaker's algorithms.
!  The first factor 2 difference is in the stiffness matrix and the second
!  one is in the inverse FT.
!---------------------------------------------------------------------------

  subroutine ForwardFE(CallBackForw, NodeX, NodeY, Cond, VIcalc, Jacobian,  & 
                       ElecNodeID, StingCMD, ParamX1, ParamX2, ParamY1,     &
                       ParamY2,InfElec, CenterNodeX, CenterNodeY, ElemArea, &
                       GetJacobian)

!DEC$ ATTRIBUTES C, DLLEXPORT :: ForwardFE
!DEC$ ATTRIBUTES ALIAS: "forwardfe" :: ForwardFE

  use GlobalForw
  implicit none

  interface
    subroutine CallBackForw(CallBackStatus)
       !DEC$ ATTRIBUTES C :: CallBackForw
       integer, intent(inout), dimension(0:1) :: CallBackStatus
    end subroutine
  end interface

  external CallBackForw


! variables and arrays passed from Delphi GUI

  integer :: GetOffsetIndex
  logical :: IsInfElectrode
  real(Rkind), intent(in),  dimension(1:gNumNodes) :: NodeX, NodeY
  real(Rkind), intent(in),  dimension(1:gNumElem)  :: Cond
  real(Rkind), intent(out), dimension(1:gNumData)  :: VIcalc
  real(Rkind), intent(out), dimension(1:gNumData*gNumParam) :: Jacobian
  integer, intent(in), dimension(1:gNumElectrodes) :: ElecNodeID
  integer, intent(in), dimension(1:gNumData*4)     :: StingCMD
  integer, intent(in), dimension(1:gNumParamX) :: ParamX1,ParamX2
  integer, intent(in), dimension(1:gNumParamY) :: ParamY1,ParamY2
  real(Rkind), intent(in), dimension(1:gNumElem)   :: CenterNodeX, CenterNodeY
  real(Rkind), intent(in), dimension(1:gNumElem*4) :: ElemArea
  integer, intent(in), dimension(1:gNumInfElectrodes) :: InfElec
  integer, intent(in) :: GetJacobian

! variables and arrays declared in this subroutine

  integer, dimension(0:1) :: CallBackStatus
  real(Rkind), parameter :: PI = 3.1415926
  integer :: i, j, k, n, kk, iAbscissa, iNode, iData 
  integer :: iElectrode, iPos, iElemX, iElemY
  integer :: PosAM, PosAN, PosBM, PosBN, nBand, StiffLen
  real(Rkind) :: abscissa, abscix, FTweight
  integer, dimension(1:4) :: LocalNodes
  real(Rkind), dimension(1:4,1:4) :: ElemStiff
  real(Rkind), allocatable, dimension(:) :: StiffMatrix, ElemStiffAll
  integer, allocatable, dimension(:) :: IndexOffset
  real(Rkind) :: vAM, vAN, vBM, vBN

  real(Rkind), allocatable, dimension(:)  :: VoverIft, DiagStiff, S
!-----------------------------------------------------------------------
! allocatable arrays are allocated on heap, but automatic arrays are
! on stack which often causes stack overflow.
  allocate(VoverIft(1:gNumNodes), stat=i)
  allocate(DiagStiff(1:gNumNodes), stat=i)
  allocate(S(1:gNumNodes), stat=i)
  DiagStiff = 0D0

! Stiffness matrix consumes different amount of memory depending on
! forward modeling method and forward solver. Finite difference method
! with a conjugate gradient solver uses the least memory.

  if (gnForwSolver == 0) then         !  Cholesky decomposition (CD)
     nBand = gnBand
  else if (gnForwSolver == 1) then    !  Conjugate gradient (CG) method
     nBand = gNumDiagonal
     allocate(IndexOffset(gNumDiagonal), stat=n)
     if (n /= 0) return                  ! ???: err message desired
     IndexOffset = (/0, 1, gNumNodeY-1, gNumNodeY, gNumNodeY+1/)
  end if

  StiffLen = nBand * gNumNodes
  allocate(StiffMatrix(StiffLen), stat=n)
  if (n /= 0) return                  ! ???: err message desired

  if (allocated(gfVoverI)) deallocate(gfVoverI)
  allocate(gfVoverI(gNumNodes*gNumElectrodes), stat=n)
  if (n /= 0) return                  ! ???: err message desired 
  
  if (GetJacobian > 0) then
     allocate(ElemStiffAll(gNumElem*16), stat=n)
     if (n /= 0) return                  ! ???: err message desired
	 ElemStiffAll = 0D0
  end if 
   
!-----------------------------------------------------------------------
!  Start forward modeling using finite element method. There are two
!  major loops in forward modeling. The first one loops through the 
!  number of abscissas, and the second one loops through the number of
!  electrodes. 

  VIcalc = 0D0
  gfVoverI = 0D0
  abscissa = 0D0
  if (GetJacobian > 0) Jacobian = 0D0

  OuterLoop : do iAbscissa = 1, gNumAbscissa
     abscissa = gfAbscissa(iAbscissa)
     abscix   = abscissa * abscissa / 12.0D0      ! used in stiffness
     FTweight = gfFTWeight(iAbscissa)      
     StiffMatrix = 0D0                
!------------------------------------------------------------------
! Construct stiffness matrix element by element (quadrilateral).
! Stiffness matrix is stored in a one-dimensional array which
! results in a fast execution. The StiffMatrix is divided into
! gNumNodes chunks. Each chunk consists of gnBand entries.
! Though only 5 of gnBand entries are non-zero, it is preferred 
! that StiffMatrix has a length of gNumNodes * gnBand since later
! Cholesky decomposition will fill the whole array.
! For conjugate gradient method, we need store only 3 non-zero 
! diagonals. 
!---------------------------------------------------------------------
     do iElemX = 1, gNumElemX            
        do iElemY = 1, gNumElemY

           call StiffnessFE(NodeX, NodeY, Cond, ElemStiff, abscix, iElemX,  &
                        iElemY, LocalNodes, CenterNodeX, CenterNodeY, ElemArea)

           ! We may save some RAM here ! due to symmetry of ElemStiff
           if (GetJacobian > 0) then
              iPos = 16 * ((iElemX-1) * gNumElemY + iElemY - 1) - 4   
              do j = 1, 4                               
                 do i = 1, 4
                    k = iPos + 4 * j + i
                    ElemStiffAll(k) = ElemStiff(i, j)
                 end do
              end do
           end if

           if (gnForwSolver == 0) then           ! Cholesky decomposition
              do i = 1,4                         ! Stores banded lower triangle 
                 do j = 1,i
                    iPos = (nBand-1)*LocalNodes(i) + LocalNodes(j)
                    StiffMatrix(iPos) = StiffMatrix(iPos) + ElemStiff(i,j)
                 end do
              end do

           else if (gnForwSolver == 1) then      ! Conjugate gradient
              do i = 1,4                         ! stores upper triangle compactly 
                 do j = i,4
                    k = LocalNodes(j) - LocalNodes(i)
                    iPos = nBand*(LocalNodes(i)-1) + GetOffsetIndex(k,IndexOffset)
                    StiffMatrix(iPos) = StiffMatrix(iPos) + ElemStiff(i,j)
                 end do
              end do
           end if

        end do
     end do
!------------------------------------------------------------------------
!  Implement boundary condition

     call BCFE(StiffMatrix,NodeX,NodeY,Cond,nBand,iAbscissa)
!-------------------------------------------------------------------------
!  Cholesky decomposition (CD) of symmetric stiffness matrix. Stiffness 
!  stays constant for different electrode locations. So only one
!  decomposition is needed for all electrodes. This is the advantage
!  of CD over conjugate gradient (CG) method
     if (gnForwSolver == 0) then        ! Cholesky decomposition
        call Decompose(StiffMatrix)

!  Scaling of stiffness matrix for SSOR preconditioning
!  Ref.:  Spitzer (1995) and Weller et al. (1996).
!  The scaling is equivalent to diagonal preconditioning
!  Turn diagonals of stiffness matrix into unity. This is achieved
!  by multiplying i-th row and column with 1/sqrt(stiffness(i,i))

     else if (gnForwSolver == 1) then 
        do i = 1, gNumNodes
           j = nBand*(i-1) + 1    ! diagonal of StiffMatix
           DiagStiff(i) = 1.0D0 / sqrt(StiffMatrix(j))
           StiffMatrix(j) = 1.0D0
        end do
        do i = 2, nBand             
           do j = 1, gNumNodes - IndexOffset(i)
              iPos = nBand*(j-1)+i
              StiffMatrix(iPos) = StiffMatrix(iPos) *   &
                         DiagStiff(j) * DiagStiff(j+IndexOffset(i))
           end do
        end do
     end if
!----------------------------------------------------------------
!  For each electrode, we will calculate an electric potential 
!  distribution produced by a point source over the entire grid.
! 
     gfVoverI = 0D0
     VoverIft = 0D0

     do iElectrode = 1, gNumElectrodes
        CallBackStatus(0) = -1 * (gNumElectrodes * iAbscissa + iElectrode)
        call CallBackForw(CallBackStatus)
        if (CallBackStatus(1) == 1) exit OuterLoop
!-------------------------------------------------------------
!  Solve for transformed V/I by by forward and then backward 
!  substitution based on the above Cholesky decomposition
!  Store all forward solutions in tranform domain for Jacobian 
!  calculation

        if (ElecNodeID(iElectrode) < 0) then
           k = (iElectrode-1) * gNumNodes
           do iNode = 1, gNumNodes
              iPos = k + iNode
              gfVoverI(iPos) = 0D0
           end do

        else

           if (gnForwSolver == 0) then
              VoverIft = 0D0
              VoverIft(ElecNodeID(iElectrode)) = 1.0D0
              call Solve(StiffMatrix,VoverIft)

           else if (gnForwSolver == 1) then
              S = 0D0
              S(ElecNodeID(iElectrode)) = DiagStiff(ElecNodeID(iElectrode)) 
              VoverIft = VoverIft / DiagStiff
              call ForwPCG(StiffMatrix,VoverIft,S,StiffLen)
              VoverIft = VoverIft * DiagStiff
           end if

!  Now store all calculated potential fields in a huge array for
!  Jacobian calculation

           k = (iElectrode-1) * gNumNodes
           do iNode = 1, gNumNodes
              iPos = k + iNode
              gfVoverI(iPos) = VoverIft(iNode)
           end do

        end if

     end do
!---------------------------------------------------------------------
!  Calculate inverse Fourier transform for calculated V/I. This is the 
!  final results of the forward modeling at current abscissa. Abscissas 
!  and weights are set in SetAbscissa.f90. We use Gauss-LaGuerre
!  quadrature (Wannamaker, 1992; LaBrecque et al., Geophysics, 1996)
!---------------------------------------------------------------------
     do iData = 1, gNumData
        kk = (iData-1) * 4
        i = StingCMD(kk + 1) * gNumNodes         !  Electrode A
        j = StingCMD(kk + 2) * gNumNodes         !  Electrode B
        k = StingCMD(kk + 3) + 1                 !  Electrode M
        n = StingCMD(kk + 4) + 1                 !  Electrode N

        if (ElecNodeID(k) > 0) then
            PosAM = i + ElecNodeID(k)
            PosBM = j + ElecNodeID(k)
            vAM = gfVoverI(PosAM)
            vBM = gfVoverI(PosBM)
        else
            vAM = 0D0
            vBM = 0D0
        end if

        if (ElecNodeID(n) > 0) then
            PosBN = j + ElecNodeID(n)
            PosAN = i + ElecNodeID(n)
            vAN = gfVoverI(PosAN)
            vBN = gfVoverI(PosBN)
        else
            vAN = 0D0
            vBN = 0D0
        end if

        VIcalc(iData) = VIcalc(iData) + FTweight * (vAM + vBN - vAN - vBM) 

     end do
!-----------------------------------------------------------------------------
!  Calculate Jacobian matrix for resistivity inversion. Jacobian will be also 
!  used in IP forward modeling and inversion.
 
     if(GetJacobian > 0) then
	    call JacobFE2(Jacobian, StingCMD, ParamX1, ParamX2, & 
                              ParamY1, ParamY2, FTweight, ElemStiffAll)
        ElemStiffAll = 0D0
     end if
!-----------------------------------------------------------------------------
  end do OuterLoop

  if (allocated(StiffMatrix)) deallocate(StiffMatrix)
  if (allocated(gfVoverI)) deallocate(gfVoverI)
  if (allocated(ElemStiffAll)) deallocate(ElemStiffAll)
  if (allocated(IndexOffset)) deallocate(IndexOffset)

  if (allocated(VoverIft)) deallocate(VoverIft)
  if (allocated(DiagStiff)) deallocate(DiagStiff)
  if (allocated(S)) deallocate(S)
!-----------------------------------------------------------------------------
  return
  end subroutine ForwardFE
