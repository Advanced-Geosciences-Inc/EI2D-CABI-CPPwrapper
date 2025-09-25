!
!  Preconditioned conjugate gradient relaxation method
!  Solving a linearized SPD system of focusing inversion.
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
!  Note that DataWeight is the reciprocal of the standard deviation of data, 
!  and Jacobian is stored row by row in a 1-D array
!
!  Jacobian times a matrix consists of only a single line of code. So I 
!  decide not to place it in a subroutine which has run-time overhead.
!  Jt * Vec has a similar scenario.
!---------------------------------------------------------------------------
  subroutine InvPCGFocus(ObsData, CalcData, DataWeight, LogCond, PriorModel,  &
                      ModelUpdate, Jacobian, ModelWeight0, ModelWeight,       &
                      Lagrange, ResIPFlag)


!DEC$ ATTRIBUTES C, DLLEXPORT :: InvPCGFocus
!DEC$ ATTRIBUTES ALIAS: "invpcgfocus" :: InvPCGFocus

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelUpdate
  real(Rkind), intent(in), dimension(1:gNumData)  :: ObsData, CalcData, DataWeight
  real(Rkind), intent(in), dimension(1:gNumParam) :: LogCond, PriorModel, ModelWeight0
  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelWeight
  real(Rkind), intent(in) :: Lagrange
  integer, intent(in) :: ResIPFlag

! local variables

  integer :: i, k, iCG
  real(Rkind) :: rVecNorm, rVecNorm0, alpha, beta
  real(Rkind) :: rtrOld, rtrNew
  real(Rkind), dimension(1:gNumData)  :: DataVec
  real(Rkind), dimension(1:gNumParam) :: rVec, pVec, zVec, ApVec
!--------------------------------------------------------------------
  k = (gNumData - 1) * gNumParam

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  allocate(gfDiagPrecond(1:gNumParam), stat=i)
  if (i /= 0) return                     ! ???: err message desired
  gfDiagPrecond = 0D0


! Obtain diagonal preconditioner
  forall (i = 1:gNumParam)
     gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*DataWeight)**2) 
  end forall

! Condition modelweight matrix with Jacobian
! Use zVec as a tmp storage
  zVec = PriorModel - LogCond 
  ModelWeight = ModelWeight0 * (gfDiagPrecond**gfModResoFactor) / (zVec*zVec + gcTinyReal)
  gfDiagPrecond = gfDiagPrecond + Lagrange * ModelWeight

! Calculate the righthand side (RHS). It takes four steps.
! 3-1. Get the weighted data residual.
  DataVec = ObsData - CalcData
  DataVec = DataVec * DataWeight * DataWeight

! 3-2. Transponsed Jacobian times a vector (DataVec)
  forall (i=1:gNumParam) rVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

! 3-3. The righthand side
  rVec = rVec + Lagrange * ModelWeight * zVec 

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

!    Calculate the kernel (Jt Wd J + a R) * pVec, i.e., Ap in CG algorithm
!    Jacobian times a vector
     forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), pVec)
     DataVec = DataVec * DataWeight * DataWeight

!    Transponsed Jacobian times a vector
     ApVec = 0D0
     forall (i=1:gNumParam) ApVec(i) = dot_product(Jacobian(i:k+i:gNumParam), DataVec)

!    Add stabilizer
     ApVec = ApVec + Lagrange * ModelWeight * pVec

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

  return
  end subroutine InvPCGFocus



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Conjugate Gradient Method to Solve a Rectangular System. 
!  Solve equations (32) and (33) in Oleg Portniaguine (1999, PhD 
!  dissertation) using his algorithm (22). 
!  Linerized equation: A x = b. Note that A has a dimension of
!  gNumParam Rows by (gNumParam+gNumData) columns. x has a dimension of
!  (gNumParam+gNumData) rows and 1 columm. The right hand side is a 
!  data vector. I split the parameter vector x into a param-domain vector 
!  and a data-domain vector
!
!  rVec   data-domain residual vector 
!  gVec   gradient vector that has two parts: gVec1 in param domain and 
!         gVec2 in data domain
!  pVec   Search direction vector that also has two parts: pVec1 in param
!         domain and pVec2 in data domain
!  ModelUpdate and NoiseUpdate are the vectors saught in the sub. These
!         two vectors are in fact one vector in the algorithm
!
!  
!-------------------------------------------------------------------------

  subroutine InvCGFocus2(ObsData, CalcData, DataWeight, LogCond,       &
                      PriorModel, ModelUpdate, Jacobian, ModelWeight0, &
                      ModelWeight, SqrtAlpha, ResIPFlag)


!DEC$ ATTRIBUTES C, DLLEXPORT :: InvCGFocus2
!DEC$ ATTRIBUTES ALIAS: "invcgfocus2" :: InvCGFocus2
!DEC$ ATTRIBUTES REFERENCE :: SqrtAlpha

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(out), dimension(1:gNumParam) :: ModelUpdate, ModelWeight
  real(Rkind), intent(in),  dimension(1:gNumData)  :: ObsData, CalcData, DataWeight
  real(Rkind), intent(in),  dimension(1:gNumParam) :: LogCond, PriorModel, ModelWeight0
  real(Rkind), intent(in),  dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(inout) :: SqrtAlpha
  integer, intent(in) :: ResIPFlag

! local variables

  integer :: i, k, iCG
  real(Rkind) :: NoiseNorm, rtrOld, rtrNew, gtgOld, gtgNew
  real(Rkind) :: tmp, ApNorm, a, b, c, step, Noise
  real(Rkind), dimension(1:gNumData)  :: rVec, pVec2, gVec2, NoiseUpdate
  real(Rkind), dimension(1:gNumParam) :: pVec1, gVec1

!-------------------------------------------------------------------------

  k = (gNumData - 1) * gNumParam
 
  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)
  allocate(gfDiagPrecond(1:gNumParam), stat=i)
  if (i /= 0) return                     ! ???: err message desired
  gfDiagPrecond = 0D0

! Determine model weight and residual vector rVec before starting loops

! Obtain diagonal of JtJ to construct model weight
  forall (i = 1:gNumParam)
     gfDiagPrecond(i) = sum((Jacobian(i:k+i:gNumParam)*DataWeight)**2) 
  end forall
! Construct modelweight matrix from diag(JtJ). gVec1 is tmp variables
!???:  gVec1 = abs(PriorModel - LogCond) 
!???:  tmp   = maxval(gVec1) + gcTinyReal
!???:  ModelWeight = gVec1 / (ModelWeight0 * (gfDiagPrecond**gfModResoFactor) * tmp) 
  ModelWeight = 1D0

! Jacobian times a vector: J(m-m0). gVec1 is a tmp variables
!???:  gVec1 = LogCond - PriorModel
!???:  forall (i=1:gNumData) rVec(i) =   &
!???:         dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), gVec1)

! Right hand side -b ==> rVec
!???:  rVec = (CalcData - ObsData - rVec) * DataWeight
  rVec = (CalcData - ObsData) * DataWeight

! Initialization

  Noise       = 0.03D0    ! 3% noise
  ModelUpdate = 0D0
  NoiseUpdate = 0D0
  SqrtAlpha   = 0D0
  gtgOld      = 1D0
  rtrOld      = 1.0D+37
  pVec1       = 0D0
  pVec2       = 0D0

  do iCG = 1, gnMaxNumIterInvCG

     rtrNew = dot_product(rVec, rVec)
     NoiseNorm = SqrtAlpha * SqrtAlpha * dot_product(NoiseUpdate, NoiseUpdate)
     if (rtrNew < 0.1 * NoiseNorm) exit

     ! Gradient vector gVec = A' * r. Note A has two parts.
     ! The first part of A is Wd * J * inverse(Wm). A'=(inverse(Wm))'*J'*Wd'
     ! gVec2 is a tmp var
     gVec2 = rVec * DataWeight

     ! Transponsed Jacobian times a data-domain vector, J' * gVec2
     gVec1 = 0D0
     forall (i=1:gNumParam) gVec1(i) = dot_product(Jacobian(i:k+i:gNumParam), gVec2)

     gVec1 = gVec1 * ModelWeight
     gVec2 = SqrtAlpha * rVec

     gtgNew = dot_product(gVec1, gVec1) + dot_product(gVec2, gVec2)

     if (gtgNew < gcTinyReal) exit
     if (rtrNew >= rtrOld) exit

     tmp = gtgNew / gtgOld    
     pVec1 = gVec1 + pVec1 * tmp
     pVec2 = gVec2 + pVec2 * tmp

     gtgOld = gtgNew
     rtrOld = rtrNew

     ! The first part of Ap. Here I reuse gVec1 and gVec2 as fVec
     gVec1 = pVec1 * ModelWeight
     ! Jacobian times gVec1
     forall (i=1:gNumData) gVec2(i) =   &
           dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), gVec1)
     gVec2 = gVec2 * DataWeight

     ! The second part of Ap
     gVec2 = gVec2 + SqrtAlpha * pVec2
     ApNorm = dot_product(gVec2, gVec2)
     
     if (iCG == 1) then             !  Find SqrtAlpha at the first iteration
        a = rtrNew * (1D0- Noise)
        b = gtgNew * (1D0- 2D0* Noise)
        c = -Noise * ApNorm
        SqrtAlpha = sqrt( (-b + sqrt(b*b - 4D0*a*c)) / (2D0*a) )
        pVec2 = SqrtAlpha * rVec
        gtgOld = dot_product(pVec1, pVec1) + dot_product(pVec2, pVec2)
        gVec2=gVec2 + SqrtAlpha * pVec2
        ApNorm = dot_product(gVec2, gVec2)
     end if

     step = dot_product(gVec2, rVec) / ApNorm
     ModelUpdate = ModelUpdate - step * pVec1
     NoiseUpdate = NoiseUpdate - step * pVec2
     rVec = rVec - step * gVec2

  end do

  ModelUpdate = ModelUpdate * ModelWeight - LogCond + PriorModel
  SqrtAlpha = SqrtAlpha * SqrtAlpha
  ModelWeight = 1D0 / (ModelWeight + gcTinyReal)

  if (allocated(gfDiagPrecond)) deallocate(gfDiagPrecond)

  return
  end subroutine InvCGFocus2
