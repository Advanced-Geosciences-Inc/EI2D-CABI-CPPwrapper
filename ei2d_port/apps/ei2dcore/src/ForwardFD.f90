
!  This is a driver routine for forward modeling using finite difference method.
!  We follow an algorithm described in Dey and Morrison (Geophys. Prosp. 1979). 
!  gnForwSolver = 0, Cholesky decomposition;  = 1, conjugate gradient.
!  Note that the only difference from Dey and Morrison is that the stiffness is
!  multiplied by a factor 2 and is made up later. I used an algorithm of
!  Discretization By Area 
!-------------------------------------------------------------------------------

  subroutine ForwardFD(CallBackForw, NodeX, NodeY, Cond, VIcalc, Jacobian, ElecNodeID, &
                       StingCMD, ParamX1, ParamX2, ParamY1, ParamY2, InfElec, &
                       GetJacobian)

!DEC$ ATTRIBUTES C, DLLEXPORT :: ForwardFD
!DEC$ ATTRIBUTES ALIAS: "forwardfd" :: ForwardFD

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
  integer, intent(in), dimension(1:gNumInfElectrodes) :: InfElec
  integer, intent(in) :: GetJacobian

! variables and arrays declared in this subroutine

  integer, dimension(0:1) :: CallBackStatus
  integer :: i, j, k, n, kk, iAbscissa, iNode, iData 
  integer :: iElectrode, iPos, iNodeX, iNodeY
  integer :: PosAM, PosAN, PosBM, PosBN, nBand, StiffLen
  real(Rkind) :: abscissa, FTweight
  real(Rkind), allocatable, dimension(:) :: StiffMatrix
  integer, allocatable, dimension(:) :: IndexOffset
  real(Rkind) :: Cl, Ct, Cp   ! Left, Top and Self coupling coefficients
  real(Rkind) :: vAM, vAN, vBM, vBN

  real(Rkind), allocatable, dimension(:)  :: dx, dy, VoverIft, DiagStiff, S
!----------------------------------------------------------------------------
! allocatable arrays are allocated on heap, but automatic arrays are
! on stack which often causes stack overflow.
  allocate(VoverIft(1:gNumNodes), stat=i)
  allocate(DiagStiff(1:gNumNodes), stat=i)
  allocate(S(1:gNumNodes), stat=i)
  allocate(dx(1:gNumElem), stat=i)
  allocate(dy(1:gNumElem), stat=i)

  dx = 0D0
  dy = 0D0
  DiagStiff = 0D0
  S = 0D0

! Stiffness matrix consumes different amount of memory depending on
! forward modeling method and forward solver. Finite difference method
! with a conjugate gradient solver uses the least memory.

  if (gnForwSolver == 0) then         !  Cholesky decomposition (CD)
     nBand = gnBand
  else if (gnForwSolver == 1) then    !  Conjugate gradient (CG) method
     nBand = gNumDiagonal
     allocate(IndexOffset(gNumDiagonal), stat=n)
     if (n /= 0) return                  ! ???: err message desired
     IndexOffset = (/0, 1, gNumNodeY/)
  end if

  StiffLen = nBand * gNumNodes
  allocate(StiffMatrix(StiffLen), stat=n)
  if (n /= 0) return                  ! ???: err message desired

  if (allocated(gfVoverI)) deallocate(gfVoverI)
  allocate(gfVoverI(gNumNodes*gNumElectrodes), stat=n)
  if (n /= 0) return                  ! ???: err message desired 
  gfVoverI = 0D0
!-----------------------------------------------------------------------
!  Start forward modeling using finite difference method. There are two
!  major loops in forward modeling. The first one loops through the 
!  number of abscissas, and the second one loops through the number of
!  electrodes. 

  VIcalc = 0D0
  if (GetJacobian > 0) Jacobian = 0D0


  do i = 1, gNumElemX
     do j = 1, gNumElemY
        k = j + (i-1) * gNumElemY
        n = j + (i-1) * gNumNodeY
        dx(k) = abs(NodeX(n+gNumNodeY) - NodeX(n))
        dy(k) = abs(NodeY(n+1)         - NodeY(n))
     end do 
  end do

  OuterLoop : do iAbscissa = 1,gNumAbscissa    ! the outermost loop 
     abscissa = gfAbscissa(iAbscissa)
     FTweight = gfFTWeight(iAbscissa)
     StiffMatrix = 0D0               
!---------------------------------------------------------------------
! Construct stiffness matrix node by node. Stiffness matrix is stored 
! in a one-dimensional array which results in a fast execution. The 
! StiffMatrix is divided into gNumNodes chunks. Each chunk consists of 
! gnBand entries. Though only 3 of gnBand entries are non-zero, it is 
! preferred that StiffMatrix has a length of gNumNodes * gnBand since 
! later Cholesky decomposition will fill the whole array. For conjugate 
! gradient method, we need store only 3 non-zero diagonals. 
! Note that nodal indices exclude boundary nodes which are handled in BCFD.
!---------------------------------------------------------------------
     do iNodeX = 2, gNumNodeX-1            
        do iNodeY = 2, gNumNodeY-1

           call StiffnessFD(iNodeX, iNodeY, NodeX, NodeY, Cond,   &
                            abscissa, Cl, Ct, Cp)
           iNode = (iNodeX-1) * gNumNodeY + iNodeY

           ! Cholesky decomposition. Stores banded lower triangle.
           if (gnForwSolver == 0) then             
               iPos = (iNode-1) * nBand + 1
               StiffMatrix(iPos)   = Cl
               iPos = iNode*nBand
               StiffMatrix(iPos-1) = Ct
               StiffMatrix(iPos)   = Cp

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
     end do
!-----------------------------------------------------------------------
!  Implement boundary condition

     call BCFD(StiffMatrix,NodeX,NodeY,Cond,nBand,iAbscissa)
     StiffMatrix = 0.50D0 * StiffMatrix 
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------
!  For each electrode, we will calculate an electric potential 
!  distribution produced by a point source over the entire grid.
! 
     VoverIft = 0D0

     do iElectrode = 1, gNumElectrodes
        CallBackStatus(0) = gNumElectrodes * iAbscissa + iElectrode
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

     if(GetJacobian > 0) call JacobFD(Jacobian, dx, dy, Cond, StingCMD, &  
         ParamX1, ParamX2, ParamY1, ParamY2, abscissa, FTweight)

!	 if (iAbscissa < gNumAbscissa) gfVoverIOld = gfVoverI
  end do OuterLoop     ! iAbscissa

  if (allocated(StiffMatrix)) deallocate(StiffMatrix)
  if (allocated(gfVoverI)) deallocate(gfVoverI)
  if (allocated(IndexOffset)) deallocate(IndexOffset)

  if (allocated(VoverIft)) deallocate(VoverIft)
  if (allocated(DiagStiff)) deallocate(DiagStiff)
  if (allocated(S)) deallocate(S)
  if (allocated(dx)) deallocate(dx)
  if (allocated(dy)) deallocate(dy)

  return 
  end subroutine ForwardFD
