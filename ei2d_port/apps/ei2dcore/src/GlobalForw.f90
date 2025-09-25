!  GlobalForw.f90 
!
  module GlobalForw
  implicit none
  save

  integer, parameter :: Rkind = 8
  integer, parameter :: Real4 = 4

  integer :: gNumData              ! Number of data (meausrements)
  integer :: gNumParamX            ! Number of X Parameters in the inverse mesh
  integer :: gNumParamY            ! Number of Y Parameters in the inverse mesh
  integer :: gNumParam             ! Number of Parameters in the inverse mesh
  integer :: gNumElectrodes        ! Number of electrodes used
  integer :: gNumInfElectrodes     ! Number of infinity electrodes
  integer :: gNumElemX             ! Number of elements in X direction
  integer :: gNumElemY             ! Number of elements in Y direction
  integer :: gNumElem              ! Total number of elements in the mesh
  integer :: gNumNodeX             ! Number of nodes in X direction
  integer :: gNumNodeY             ! Number of nodes in Y direction
  integer :: gNumNodes             ! Total number of nodes in the mesh

  integer :: gnForwModMeth         ! Forward Modeling method  0=FD  1=FE
  integer :: gnForwSolver          ! Forward matrix equation solver  0=CD  1=CG
  integer :: gnInvMethod           ! Inversion Method. 0=Forward Only, 1=Least square,
                                     ! 2=Occam's, 3=Robust, 4=Focusing, 5=Stochastic
                                     ! 6=Time-elapse, 7=Mixed norm
  integer :: gnBand                ! half bandwidth of stiffness matrix
  integer :: gNumAbscissa          ! number of abscissas of numerical integration for inverse FT
  integer :: gNumDiagonal          ! Number of unique diagonals in stiffness matrix
  integer :: gnBCType              ! 0 = Dirichlet, 1 = mixed
  integer :: gnForwCGIter          ! max number of iteration of forward CG (200)
  real(Rkind) :: gfForwCGResid     ! forward CG stop residual (1d-6)

! global arrays
  real(Rkind), allocatable, dimension(:) :: &
               gfAbscissa,  &      ! Numerical integration points
               gfFTWeight,  &      ! Numerical integration weights
               gfStiff             ! 1-D Stiffness Matrix of compact and efficient storage

! Eta factor for mixed B.C.
  real(Rkind), allocatable, dimension(:, :) :: gfEtaLeft, gfEtaRight, gfEtaBottom  

! Array to store complete forward solution at a single abscissa. It is 
! updated from abscissa to abscissa. For efficient access, I use a 1D array
!  real(RKind), allocatable, dimension(:) :: gfVoverIOld 
  real(RKind), allocatable, dimension(:) :: gfVoverI 

  end module GlobalForw
