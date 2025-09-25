!
!  Preconditioned conjugate gradient relaxation method
!  Solving a linearized SPD system from L1 norm robust inversion.
!
!  Note that nonlinear inversion consists of two nested loops. The
!  outer loop updates model bit by bit and evaluates the convergence.
!  The inner loop is an iterative relaxation loop which searches for the 
!  actual model update. The alternative is to do a Cholesky decomposition
!  which is an exact method without any loop but pretty slow.
! 
!  Ref.: Golub and van Loan, 1989, 1996. Matrix Computation
!        Barrett et al., 1994. Temperates for the solution of ...
!
!  DataWeight is the reciprocal of the standard deviation of data, 
!  and Jacobian is stored row by row in a 1-D array
!
!  Jacobian times a matrix consists of only a single line of code. So I 
!  decide not to place it in a subroutine which has run-time overhead.
!  Jt * Vec has a similar scenario.

!  For nonlinear IP inversion, 
!     ObsData = -(V/I)res / [(V/I)ip]**2,  a scale factor of IP Jacobian
!     CalcData = ObsData - ClacData
!
!-------------------------------------------------------------------------------

  subroutine InvPCGRobust(ObsData, CalcData, DataWeight, ModelParam, PriorModel,   &
                          ModelUpdate, Jacobian, RoughX0, RoughY0, RoughX, RoughY,  &
                          Lagrange, ResIPFlag, IterNum)



!DEC$ ATTRIBUTES C, DLLEXPORT :: InvPCGRobust
!DEC$ ATTRIBUTES ALIAS: "invpcgrobust" :: InvPCGRobust
!DEC$ ATTRIBUTES C, DLLIMPORT :: RoughV0
!DEC$ ATTRIBUTES ALIAS: "roughv0" :: RoughV0

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelUpdate
  real(Rkind), intent(in), dimension(1:gNumData)  :: ObsData, CalcData, DataWeight
  real(Rkind), intent(in), dimension(1:gNumParam) :: ModelParam, PriorModel
  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumParam) :: RoughX0, RoughY0
  real(Rkind), intent(out), dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(in) :: Lagrange
  integer, intent(in) :: ResIPFlag, IterNum

! local variables

  integer :: i, k, iCG, nospace
  real(Rkind) :: rVecNorm, rVecNorm0, alpha, beta, gamma
  real(Rkind) :: rtrOld, rtrNew, EpsilonD, EpsilonM

  real(Rkind), allocatable, dimension(:) :: DataVec, ReweightData, DataVec2
  real(Rkind), allocatable, dimension(:) :: rVec, pVec, zVec, ApVec, DampVec
  ! R matrices in Farquharson (SEGJ 2003) R_i, i = s, x, y
  real(Rkind), allocatable, dimension(:) :: ReweightMx, ReweightMy, ReweightMs
!---------------------------------------------------------------------------------------
! allocatable arrays are allocated on heap, but automatic arrays are
! on stack which often causes stack overflow.
  allocate(DataVec(1:gNumData), stat=i)
  allocate(ReweightData(1:gNumData), stat=i)
  allocate(DataVec2(1:gNumData), stat=i)
  allocate(rVec(1:gNumParam), stat=i)
  allocate(pVec(1:gNumParam), stat=i)
  allocate(zVec(1:gNumParam), stat=i)
  allocate(ApVec(1:gNumParam), stat=i)
  allocate(ReweightMx(1:gNumParam), stat=i)
  allocate(ReweightMy(1:gNumParam), stat=i)
  allocate(ReweightMs(1:gNumParam), stat=i)
  allocate(DampVec(1:gNumParam), stat=i)

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  allocate(gfDiagPrecond(1:gNumParam), stat=nospace)    ! global

  DataVec = 0D0
  DataVec2 = 0D0
  ReweightData = 0D0
  rVec = 0D0
  pVec = 0D0
  zVec = 0D0
  ApVec = 0D0
  ReweightMx = 0D0
  ReweightMy = 0D0
  ReweightMs = 0D0
  gfDiagPrecond  = 0D0

  k = (gNumData - 1) * gNumParam

! Obtain ReweightData - Rd matrix in Farquharson and Oldenburg (1998)
  if (ResIPFlag == 0) then 
     DataVec = (ObsData - CalcData) * DataWeight
  else
     DataVec = CalcData * DataWeight
  end if

! R matrix in Farquharson (1998, 2003)
  if (IterNum == 1) then   
     ReweightData = 1.0D0
  else if (gnRobustMethod == 0) then
     ! Lp - norm 
     do i = 1, gNumData
        gamma = max(gfEpsilonD, abs(DataVec(i)))
        ReweightData(i) = 1.0D0 / gamma
     end do
  else if (gnRobustMethod == 1) then
     ! Huber M-Estimator 
     do i = 1, gNumData
        gamma = max(gfEpsilonD, abs(DataVec(i)))
        ReweightData(i) = 2.0D0 * gfEpsilonD / gamma 
     end do
  else if (gnRobustMethod == 2) then
     ! Ekblom 
     EpsilonD = gfEpsilonD * gfEpsilonD
     ReweightData = 1.0D0 / sqrt(DataVec * DataVec + EpsilonD)
  else if (gnRobustMethod == 3) then
     ! Focusing 
     EpsilonD = gfEpsilonD * gfEpsilonD
     ReweightData = 2.0D0 * EpsilonD /(DataVec * DataVec + EpsilonD)
  end if

  if (ResIPFlag == 0) then
     DataVec2 = DataWeight * sqrt(ReweightData)
  else
     DataVec2 = ObsData * DataWeight * sqrt(ReweightData)
  end if

! Obtain diagonal preconditioner
! Roughness part will be calculated in 4-3 below
  forall (i = 1:gNumParam)
     gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*DataVec2)**2) 
  end forall

!--------------------------------------------------------------------
! Calculate the righthand side (RHS). It takes four steps.
! 4-1. Get the weighted data residual. Note that DataVec is data misfit
!      weighted by DataWeight (standard deviation). See DataVec above.
  if (ResIPFlag == 0) then
     DataVec = DataVec * DataWeight * ReweightData
  else
     DataVec = DataVec * DataWeight * ReweightData * ObsData 
  end if

! 4-2. Transponsed Jacobian times a vector (DataVec)
  forall (i=1:gNumParam) rVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

! 4-3. Roughness times a vector. Use zVec as a tmp storage
!      Note that zVec = PriorModel - ModelParam. See above
!      I also calculate Model Conditioner (ReweightModel) and 
!      Roughness parf of diagonal preconditioner (gfDiagPrecond)

! Condition roughness with Jacobian and R matrix in Farquharson. This step may be redundant
  DampVec = gfDiagPrecond ** gfModResoFactor
  RoughX = RoughX0 * DampVec
  RoughY = RoughY0 * DampVec

  zVec = PriorModel - ModelParam
  if (dot_product(zVec, zVec) < gcTinyReal) then
     nospace = 0   ! nospace is a tmp storage
  else
     nospace = 1
     call RoughL1X(zVec, RoughX, ReweightMx)
     call RoughL1Y(zVec, RoughY, ReweightMy)
  end if

! R matrix in Farquharson (1998, 2003)
  if ((IterNum==1) .or. (nospace==0))then
     ReweightMx = 1.0D0
     ReweightMy = 1.0D0
     ReweightMs = 1.0D0     
  else if (gnRobustMethod == 0) then
     ! Lp - norm 
     do i = 1, gNumParam
        gamma = max(gfEpsilonM, abs(ReweightMx(i)))
        ReweightMx(i) = 1.0D0 / gamma
        gamma = max(gfEpsilonM, abs(ReweightMy(i)))
        ReweightMy(i) = 1.0D0 / gamma
        gamma = max(gfEpsilonM, abs(zVec(i)))
        ReweightMs(i) = 1.0D0 / gamma
     end do
  else if (gnRobustMethod == 1) then
     ! Huber M-Estimator 
     do i = 1, gNumParam
        gamma = max(gfEpsilonM, abs(ReweightMx(i)))
        ReweightMx(i) = 2.0D0 * gfEpsilonM / gamma 
        gamma = max(gfEpsilonM, abs(ReweightMy(i)))
        ReweightMy(i) = 2.0D0 * gfEpsilonM / gamma 
        gamma = max(gfEpsilonM, abs(zVec(i)))
        ReweightMs(i) = 2.0D0 * gfEpsilonM / gamma 
     end do
  else if (gnRobustMethod == 2) then
     ! Ekblom 
     EpsilonM = gfEpsilonM * gfEpsilonM
     ReweightMx = 1.0D0 / sqrt(ReweightMx * ReweightMx + EpsilonM)
     ReweightMy = 1.0D0 / sqrt(ReweightMy * ReweightMy + EpsilonM)
     ReweightMs = 1.0D0 / sqrt(zVec * zVec + EpsilonM)
  else if (gnRobustMethod == 3) then
     ! Focusing 
     EpsilonM = gfEpsilonM * gfEpsilonM
     ReweightMx = 2.0D0 * EpsilonM /(ReweightMx * ReweightMx + EpsilonM)
     ReweightMy = 2.0D0 * EpsilonM /(ReweightMy * ReweightMy + EpsilonM)
     ReweightMs = 2.0D0 * EpsilonM /(zVec * zVec + EpsilonM)
  end if

! Condition roughness with R matrix in Farquharson. This step may be redundant
  RoughX = RoughX * sqrt(ReweightMx)
  RoughY = RoughY * sqrt(ReweightMy)
  call RoughV3(zVec, RoughX, RoughY, pVec)

! zVec below is the diag(RtR)
  gfDiagPrecond = gfDiagPrecond + Lagrange * (zVec + ReweightMs)

! 4-4. The righthand side. pVec = RtR * (PriorModel - ModelParam)
  rVec = rVec + Lagrange * (pVec + ReweightMs * (PriorModel - ModelParam))

! Initial residual
  rVecNorm0 = sqrt(dot_product(rVec, rVec))
  if (rVecNorm0 < gcTinyReal) rVecNorm0 = gcTinyReal

  ModelUpdate = 0D0
  rtrNew      = 0D0

  if (ResIPFlag == 0) then
     ReweightData = ReweightData * (DataWeight**2) 
  else
     ReweightData = ReweightData * ((DataWeight * ObsData)**2)
  end if

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

!    Calculate the kernel (Jt Rd Wd J + a Rm R) * pVec, i.e., Ap in CG algorithm
!    Jacobian times a vector
     forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), pVec)
     DataVec = DataVec * ReweightData

!    Transponsed Jacobian times a vector
     ApVec = 0D0
     forall (i=1:gNumParam) ApVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

!    zVec = RtR * pVec, pVec remains unchanged.
     call RoughV0(pVec, RoughX, RoughY, zVec)

!    Add roughness term
     ApVec = ApVec + Lagrange * (zVec + ReweightMs * pVec)  

     alpha = rtrNew / dot_product(pVec, ApVec)

!    Update ModelUpdate and search direction (residual) vector
     ModelUpdate = ModelUpdate + alpha * pVec
     rVec        = rVec - alpha * ApVec

!    Evaluate the stop residual
     rVecNorm = dot_product(rVec, rVec)
     rVecNorm = sqrt(rVecNorm) / rVecNorm0
     if (rVecNorm < gcStopResidual) exit

  end do

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  if (allocated(DataVec)) deallocate(DataVec)
  if (allocated(ReweightData)) deallocate(ReweightData)
  if (allocated(DataVec2)) deallocate(DataVec2)
  if (allocated(rVec)) deallocate(rVec)
  if (allocated(pVec)) deallocate(pVec)
  if (allocated(zVec)) deallocate(zVec)
  if (allocated(ApVec)) deallocate(ApVec)
  if (allocated(ReweightMx)) deallocate(ReweightMx)
  if (allocated(ReweightMy)) deallocate(ReweightMy)
  if (allocated(ReweightMs)) deallocate(ReweightMs)
  if (allocated(DampVec)) deallocate(DampVec)

  return
  end subroutine InvPCGRobust
