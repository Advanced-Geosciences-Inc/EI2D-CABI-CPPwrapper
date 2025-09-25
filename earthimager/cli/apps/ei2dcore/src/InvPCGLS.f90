!
!  Preconditioned conjugate gradient relaxation method
!  Solving a linearized SPD system from damped least squares inversion.
!  Note that nonlinear inversion consists of two nested loops. The
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
!
!  For nonlinear IP inversion, 
!     ObsData = -(V/I)res / [(V/I)ip]**2,  a scale factor of IP Jacobian
!     CalcData = ObsData - ClacData
!
!-----------------------------------------------------------------------------

  subroutine InvPCGLS(ObsData, CalcData, DataWeight, ModelParam, PriorModel,    &
                    ModelUpdate, Jacobian, DampingFactor, ResIPFlag, IterNum)

!DEC$ ATTRIBUTES C, DLLEXPORT :: InvPCGLS
!DEC$ ATTRIBUTES ALIAS: "invpcgls" :: InvPCGLS

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelUpdate
  real(Rkind), intent(in), dimension(1:gNumData)  :: ObsData, CalcData, DataWeight
  real(Rkind), intent(in), dimension(1:gNumParam) :: ModelParam, PriorModel
  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in) :: DampingFactor
  integer, intent(in) :: ResIPFlag, IterNum

! local variables

  integer :: i, k, iCG
  real(Rkind) :: rVecNorm, rVecNorm0, alpha, beta
  real(Rkind) :: rtrOld, rtrNew
  real(Rkind), allocatable, dimension(:) :: DataVec, DataWeight2
  real(Rkind), allocatable, dimension(:) :: rVec, pVec, zVec, ApVec, ModelDamping
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
  allocate(ModelDamping(1:gNumParam), stat=i)

  DataVec = 0D0
  DataWeight2 = 0D0
  rVec = 0D0
  pVec = 0D0
  zVec = 0D0
  ApVec = 0D0
  ModelDamping = 0D0
  gfDiagPrecond  = 0D0

  if (ResIPFlag == 0) then
     DataWeight2 = DataWeight ** 2
     DataVec = DataWeight            ! for gfDiagPrecond 
  else
     DataWeight2 = (DataWeight * ObsData) ** 2
     DataVec = DataWeight * ObsData  ! for gfDiagPrecond 
  end if

  k = (gNumData - 1) * gNumParam

! Obtain diagonal preconditioner
! ???: Like roughness factor, damping factor can be an array and may vary
! from parameter to parameter depending on the diag(JtJ)

  forall (i = 1:gNumParam)
     gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*DataVec)**2) 
  end forall

  ModelDamping = DampingFactor * (gfDiagPrecond ** gfModResoFactor)
!  gfDiagPrecond  = gfDiagPrecond + DampingFactor   
  gfDiagPrecond  = gfDiagPrecond + ModelDamping   

! Calculate the righthand side (RHS). First, get the weighted data residual.
  if (ResIPFlag == 0) then
     DataVec = (ObsData - CalcData) * DataWeight2
  else
     DataVec = CalcData * DataWeight**2
     DataVec = DataVec * ObsData
  end if

! Transponsed Jacobian times a vector (DataVec)

  forall (i=1:gNumParam) rVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

! The righthand side. 
  pVec = PriorModel - ModelParam
  zVec = ModelDamping * pVec
  rVec = rVec + zVec 

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

!    Calculate the kernel (Jt Wd J + a I) times pVec, i.e., Ap in CG algorithm
!    First, Jacobian times pVec
     forall (i=1:gNumData) DataVec(i) =   &
             dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), pVec)
     DataVec = DataVec * DataWeight2

!    Transponsed Jacobian times a vector
     ApVec = 0D0
     forall (i=1:gNumParam) ApVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

!    Add damping term
     ApVec = ApVec + ModelDamping * pVec
!     ApVec = ApVec + DampingFactor * pVec

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
  if (allocated(ModelDamping)) deallocate(ModelDamping)

  return
  end subroutine InvPCGLS
