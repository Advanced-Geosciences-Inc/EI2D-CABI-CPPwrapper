!
!  Preconditioned conjugate gradient relaxation method
!  Solving a linearized SPD system from smooth model Occam's inversion.
!
!  Nonlinear inversion consists of two nested loops. The
!  outer loop updates model bit by bit and evaluates the convergence.
!  The inner loop is an iterative relaxation loop which searches for the 
!  actual model update. The alternative is to do a Cholesky decomposition
!  which is an exact method without any loop but pretty slow.
! 
!  Ref.: Golub and van Loan, 1989, 1996. Matrix Computation
!        Barrett et al., 1994. Temperates for the solution of ...
!
!  Note that DataWeight is the reciprocal of the standard deviation of data, 
!  and Jacobian is stored row by row in a 1-D array
!
!  Jacobian times a matrix consists of only a single line of code. So I 
!  decide not to place it in a subroutine which has run-time overhead.
!  Jt * Vec has a similar scenario.

!  For nonlinear IP inversion, 
!     ObsData = -(V/I)res / [(V/I)ip]**2,  a scale factor of IP Jacobian
!     CalcData = gfIPMeas - gfIPCalc

!  For linear IP inversion w/ positivity constraint, 
!     ObsData = IPJacobFactor = -1.0 / gfVICalc
!     CalcData = gfIPMeas - gfIPCalc

!--------------------------------------------------------------------------------
  subroutine InvPCGOC(ObsData, CalcData, DataWeight, ModelParam, PriorModel,    &
                      ModelUpdate, Jacobian, RoughX0, RoughY0, RoughX, RoughY,  &
                      Lagrange, DampFactor, ResIPFlag, IterNum)


!DEC$ ATTRIBUTES C, DLLEXPORT :: InvPCGOC
!DEC$ ATTRIBUTES ALIAS: "invpcgoc" :: InvPCGOC
!DEC$ ATTRIBUTES C, DLLIMPORT :: RoughV0
!DEC$ ATTRIBUTES ALIAS: "roughv0" :: RoughV0

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelUpdate
  real(Rkind), intent(in), dimension(1:gNumData)   :: ObsData, CalcData, DataWeight
  real(Rkind), intent(in), dimension(1:gNumParam)  :: ModelParam, PriorModel
  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumParam) :: RoughX0, RoughY0
  real(Rkind), intent(out), dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(in) :: Lagrange, DampFactor
  integer, intent(in) :: IterNum
  integer, intent(in) :: ResIPFlag   ! 0=Res, 1=NonlinearIP, 2=LinearIPPos

! local variables

  integer :: i, k, iCG
  logical bPos
  real(Rkind) :: rVecNorm, rVecNorm0, alpha, beta
  real(Rkind) :: rtrOld, rtrNew
  real(Rkind), allocatable, dimension(:) :: DataVec, DataWeight2
  real(Rkind), allocatable, dimension(:) :: rVec, pVec, zVec, ApVec, DampVec, PosVec
!-----------------------------------------------------------------------------
  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  allocate(gfDiagPrecond(1:gNumParam), stat=i)

! allocatable arrays are allocated on heap, but automatic arrays are
! on stack which often causes stack overflow.
  allocate(DataVec(1:gNumData), stat=i)
  allocate(DataWeight2(1:gNumData), stat=i)
  allocate(rVec(1:gNumParam), stat=i)
  allocate(pVec(1:gNumParam), stat=i)
  allocate(zVec(1:gNumParam), stat=i)
  allocate(ApVec(1:gNumParam), stat=i)
  allocate(DampVec(1:gNumParam), stat=i)
  allocate(PosVec(1:gNumParam), stat=i)

  bPos = (gnIPPosMethod == 1) .AND. (ResIPFlag == 1)
  PosVec = 1.0D0
  if (bPos) call GetPositivityVec(ModelParam, PosVec)
  DataVec = 0D0
  DataWeight2 = 0D0
  rVec = 0D0
  pVec = 0D0
  zVec = 0D0
  ApVec = 0D0
  DampVec = 0D0
  gfDiagPrecond  = 0D0

  k = (gNumData - 1) * gNumParam
  if (ResIPFlag == 0) then
     DataWeight2 = DataWeight ** 2
     DataVec = DataWeight             ! for gfDiagPrecond 
  else
     DataWeight2 = (DataWeight * ObsData) ** 2
     DataVec = DataWeight * ObsData   ! for gfDiagPrecond 
  end if

! Obtain diagonal preconditioner, Roughness part will be calculated in 4-3 below
  forall (i = 1:gNumParam) gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*DataVec)**2) 
  if (bPos) gfDiagPrecond = (PosVec ** 2) * gfDiagPrecond 

! Condition roughness matrix with Jacobian
  if (ResIPFlag == 1) then
      DampVec = gfDiagPrecond ** (1.5*gfModResoFactor)
  else
      DampVec = gfDiagPrecond ** gfModResoFactor
  end if 
  RoughX  = RoughX0 * DampVec
  RoughY  = RoughY0 * DampVec
  DampVec = DampFactor * DampVec

! Calculate the righthand side (RHS). It takes four steps.
! 4-1. Get the weighted data residual.
  if (ResIPFlag == 0) then
     DataVec = (ObsData - CalcData) * DataWeight2
  else
     DataVec = CalcData * (DataWeight**2)
     DataVec = DataVec * ObsData
  end if

! 4-2. Transponsed Jacobian times a vector (DataVec)
  forall (i=1:gNumParam) rVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)
  if (bPos) rVec = rVec * PosVec

! 4-3. Roughness times a vector. Use zVec as a tmp storage of input and output
!      Actual augument 0 flags calculating roughness part of gfDiagPrecond.
  zVec = PriorModel - ModelParam
! pVec = RtR * zVec, zVec = diag{RtR},  
  call RoughV2(zVec, RoughX, RoughY, pVec)
! Update diagonal preconditioner
  gfDiagPrecond = gfDiagPrecond + Lagrange * zVec + DampVec

! 4-4. The righthand side
  rVec = rVec + Lagrange * pVec 

! Initial residual
  rVecNorm0 = sqrt(dot_product(rVec, rVec))
  if (rVecNorm0 < gcTinyReal) rVecNorm0 = gcTinyReal

  ModelUpdate = 0D0
  rtrNew      = 0D0

! conjugate gradient relaxation (inner loop)

  do iCG = 1, gnMaxNumIterInvCG

!    Preconditioning
     zVec = rVec / gfDiagPrecond

     rtrOld = rtrNew
     rtrNew = dot_product(rVec, zVec)

     if (iCG == 1) then
        pVec = zVec
     else
        beta = rtrNew / rtrOld
        pVec = zVec + beta * pVec
     end if

!    Calculate the kernel (Jt Wd J + a*R + lamda*I) * pVec, i.e., Ap in CG algorithm
!    Jacobian times a vector
     if (bPos) then
        zVec = PosVec * pVec
        forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), zVec)
	   else
        forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), pVec)
     end if
     DataVec = DataVec * DataWeight2

!    Transponsed Jacobian times a vector
     ApVec = 0D0
     forall (i=1:gNumParam) ApVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)
     if (bPos) ApVec = ApVec * PosVec

     ApVec = ApVec + DampVec * pVec

!    zVec = RtR * pVec. Use zVec as a tmp storage, pVec remain unchanged
     call RoughV0(pVec, RoughX, RoughY, zVec)

!    Add roughness term
     ApVec = ApVec + Lagrange * zVec

     alpha = rtrNew / dot_product(pVec, ApVec)

!    Update ModelUpdate and search direction (residual) vector
     ModelUpdate = ModelUpdate + alpha * pVec
     rVec        = rVec - alpha * ApVec

!    Evaluate the stop residual
     rVecNorm    = dot_product(rVec, rVec)
     rVecNorm    = sqrt(rVecNorm) / rVecNorm0
     if (rVecNorm < gcStopResidual) exit

  end do

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  if (allocated(DataVec)) deallocate(DataVec)
  if (allocated(DataWeight2)) deallocate(DataWeight2)
  if (allocated(rVec)) deallocate(rVec)
  if (allocated(pVec)) deallocate(pVec)
  if (allocated(zVec)) deallocate(zVec)
  if (allocated(ApVec)) deallocate(ApVec)
  if (allocated(DampVec)) deallocate(DampVec)
  if (allocated(PosVec)) deallocate(PosVec)

  return
  end subroutine InvPCGOC


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Linear IP inversion: solve for x = m - m0, then m = x + m0
! [Jt * Wd * J + Lagrange * R + DampFac * I] (m - m0) = Jt * Wd * (d - J * m0)
! ------------------------------------------  ------    ---------------------
!            A                                   x    =         b

  subroutine InvPCGOCIP(IPMeas, IPCalc, VICalc, DataWeight, ModelParam, PriorModel,  &
                 Jacobian, RoughX0, RoughY0, RoughX, RoughY, Lagrange, DampFactor)
 
!DEC$ ATTRIBUTES C, DLLEXPORT :: InvPCGOCIP
!DEC$ ATTRIBUTES ALIAS: "invpcgocip" :: InvPCGOCIP
!DEC$ ATTRIBUTES C, DLLIMPORT :: RoughV0
!DEC$ ATTRIBUTES ALIAS: "roughv0" :: RoughV0

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelParam
  real(Rkind), intent(in),  dimension(1:gNumData)  :: IPMeas, IPCalc, VICalc, DataWeight
  real(Rkind), intent(in),  dimension(1:gNumParam) :: PriorModel
  real(Rkind), intent(in),  dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in),  dimension(1:gNumParam) :: RoughX0, RoughY0
  real(Rkind), intent(out), dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(in) :: Lagrange, DampFactor

! local variables

  integer, parameter :: NumIterLinearIP = 200
  real(Rkind), parameter :: StopResidual = 1.0D-5

  integer :: i, k, iCG
  real(Rkind) :: rVecNorm, rVecNorm0, alpha, beta
  real(Rkind) :: rtrOld, rtrNew
  real(Rkind), allocatable, dimension(:) :: DataVec, NewDataWeight, NewDataWeight2
  real(Rkind), allocatable, dimension(:) :: rVec, pVec, zVec, ApVec, DampVec

!-----------------------------------------------------------------------------  
! allocatable arrays are allocated on heap, but automatic arrays are
! on stack which often causes stack overflow.
  allocate(DataVec(1:gNumData), stat=i)
  allocate(NewDataWeight(1:gNumData), stat=i)
  allocate(NewDataWeight2(1:gNumData), stat=i)
  allocate(rVec(1:gNumParam), stat=i)
  allocate(pVec(1:gNumParam), stat=i)
  allocate(zVec(1:gNumParam), stat=i)
  allocate(ApVec(1:gNumParam), stat=i)
  allocate(DampVec(1:gNumParam), stat=i)

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  allocate(gfDiagPrecond(1:gNumParam), stat=i)

  DataVec = 0D0
  NewDataWeight = 0D0
  NewDataWeight2 = 0D0
  rVec = 0D0
  pVec = 0D0
  zVec = 0D0
  ApVec = 0D0
  gfDiagPrecond  = 0D0

  k = (gNumData - 1) * gNumParam
  NewDataWeight  = -1.0D0 *  DataWeight / (VICalc + gcTinyReal)
  NewDataWeight2 = NewDataWeight ** 2

! Calculate the righthand side (RHS) r0 = b - A*x0. See definition of A, x & b above
! x0 = m - m0, the CG starts from x0 = 0, so A * x0 = 0; r0 = b
! Jacobian scaling has a sign flip, i.e., J_ip = J_res / (-gfVICalc)
! b = Jt Wd [d - J m0] with Jacobian scaling, IPCalc = J m0
  DataVec = (IPCalc - IPMeas) * NewDataWeight2 * VICalc

! Transponsed Jacobian times a vector (DataVec), get r0 = b
  forall (i=1:gNumParam) rVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

! Obtain diagonal preconditioner diag(Jt Wd J)
  forall (i = 1:gNumParam)
     gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*NewDataWeight)**2) 
  end forall
  
! Roughness reweighting with  diag(Jt Wd J), 
  DampVec = gfDiagPrecond ** (1.25*gfModResoFactor) 
  RoughX  = RoughX0 * DampVec
  RoughY  = RoughY0 * DampVec
  DampVec = DampFactor * DampVec
! Obtain diag(Rt R). Use zVec as a tmp storage of input and output
  zVec = PriorModel  
  call RoughV2(zVec, RoughX, RoughY, pVec)
  gfDiagPrecond = gfDiagPrecond + Lagrange * zVec  + DampVec

  ModelParam  = 0
! Initial residual
  rVecNorm0 = sqrt(dot_product(rVec, rVec))
  if (rVecNorm0 < gcTinyReal) rVecNorm0 = gcTinyReal
  rtrNew = 0D0

! conjugate gradient relaxation (inner loop)
  do iCG = 1, NumIterLinearIP

!    Preconditioning
     zVec = rVec / gfDiagPrecond

     rtrOld = rtrNew
     rtrNew = dot_product(rVec, zVec)

     if (iCG == 1) then
        pVec = zVec
     else
        beta = rtrNew / rtrOld
        pVec = zVec + beta * pVec
     end if

!    Calculate the kernel (Jt Wd J + a R) * pVec, i.e., Ap in CG algorithm
!    Jacobian times a vector
     forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), pVec)
     DataVec = DataVec * NewDataWeight2

!    Transponsed Jacobian times a vector
     ApVec = 0D0
     forall (i=1:gNumParam) ApVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

     ApVec = ApVec + DampVec *  pVec
     
     ! zVec = RtR * pVec, pVec remains unchanged.
     call RoughV0(pVec, RoughX, RoughY, zVec)

!    Add roughness term
     ApVec = ApVec + Lagrange * zVec

     alpha = rtrNew / dot_product(pVec, ApVec)

!    Update ModelUpdate and search direction (residual) vector
     ModelParam = ModelParam + alpha * pVec
     rVec       = rVec - alpha * ApVec

!    Evaluate the stop residual
     rVecNorm    = dot_product(rVec, rVec)
     rVecNorm    = sqrt(rVecNorm) / rVecNorm0
     if (rVecNorm < StopResidual) exit
  end do
  ModelParam = ModelParam + PriorModel
  
  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  if (allocated(DataVec)) deallocate(DataVec)
  if (allocated(NewDataWeight)) deallocate(NewDataWeight)
  if (allocated(NewDataWeight2)) deallocate(NewDataWeight2)
  if (allocated(rVec)) deallocate(rVec)
  if (allocated(pVec)) deallocate(pVec)
  if (allocated(zVec)) deallocate(zVec)
  if (allocated(ApVec)) deallocate(ApVec)
  if (allocated(DampVec)) deallocate(DampVec)

  return
  end subroutine InvPCGOCIP
  