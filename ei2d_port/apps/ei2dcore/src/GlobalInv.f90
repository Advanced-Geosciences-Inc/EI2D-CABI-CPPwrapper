!  GlobalInv.f90 
!
  module GlobalInv
  implicit none
  save

  integer, parameter :: Rkind = 8
  real(Rkind), parameter :: gcTinyReal = 1.0D-19
  real(Rkind), parameter :: gcStopResidual = 1.0D-3    ! nonlinear CG loop
! upper limit of IP positivity constraint, the lower limit is set to 0
  real(Rkind), parameter :: gcMaxIP = 1D0              

  integer :: gNumData              ! Number of data (meausrements)
  integer :: gNumElemX             ! Number of elements in X (horizontal) direction
  integer :: gNumElemY             ! Number of elements in Y (vertical) direction
  integer :: gNumParamX            ! Number of X Parameters in the inverse mesh
  integer :: gNumParamY            ! Number of Y Parameters in the inverse mesh
  integer :: gNumParam             ! Number of Parameters in the inverse mesh
  integer :: gnMaxNumIterInvCG     ! Max number of iteration of inverse CG solver

  ! Inversion Method. 0=Forward Only, 1=Least square,
  ! 2=Occam's, 3=Robust, 4=Focusing, 5=Time-elapse, 6=Mixed norm
  integer :: gnInvMethod           
  integer :: gnIPInvMethod           
  integer :: gnIPPosMethod         ! IP Positivity Constraint method 0=No, 1=Yes
  integer :: gnRobustMethod        ! 0 = Lp, 1 = M-Est, 2 = Ekblom, 3 = Focusing

  real(Rkind) :: gfModResoFactor   ! This factor is intended to produce a uniform model resolution
  real(Rkind) :: gfEpsilonD        ! Square of a small number in Ekblom robust algorithm, D = Data
  real(Rkind) :: gfEpsilonM        ! Square of a small number in Ekblom robust algorithm, D = Model

! global arrays
  real(Rkind), allocatable, dimension(:) :: gfDiagPrecond  ! Diagonal preconditioner of Inv PCG

  end module GlobalInv

