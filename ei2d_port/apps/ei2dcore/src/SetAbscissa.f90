!
!  Set abscissas and weights for inverse Fourier transform. We use
!  Gauss-LaGuerre Quadrature (Wannamaker, 1992; LaBrecque et al.,
!  Geophysics, v. 61, p538, 1996)
!  ???: Important: There are two types of LaGuerre weights in Abramowitz
!  and Stegun (1970): W and W * exp(X). We listed the latter here.
!  Carnahan, B., Luther, H.A., and Wilkes, J.O., 1969, 
!      Applied Numerical Methods. p. 104-105

  subroutine SetAbscissa(ForwAccuracy, MinTxRxSep)

  use GlobalForw
  implicit none

  integer, intent(in)     :: ForwAccuracy
  real(Rkind), intent(in) :: MinTxRxSep
  real(Rkind), parameter  :: PI = 3.14159265359
  integer :: i, j, nospace, NumAbscG, NumAbscL
  real(Rkind) :: GaussLimit, tmp, InvPI

!------------------------------------------------------------------
! First, figure out number of abscissas to use

  if (ForwAccuracy == 0) then
     gNumAbscissa = 8
  else if (ForwAccuracy == 2) then
     gNumAbscissa = 16
  else
     gNumAbscissa = 12
  end if

! Now initialize abscissas and weights of inverse Fourier transform
! gNumAbscissa = gNumAbscG (Gauss points) + gNumAbscL (LaGuerre points)
! Weight Factors of Gauss-Laguerre Quadrature below is the product of 
! based-version factor times exp(z).	

  if (allocated(gfAbscissa)) deallocate(gfAbscissa)
  allocate(gfAbscissa(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  gfAbscissa = 0.0D0

  if (allocated(gfFTWeight)) deallocate(gfFTWeight)
  allocate(gfFTWeight(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  gfFTWeight = 0.0D0

  select case (gNumAbscissa)
  case (6)
     NumAbscG   = 4
     NumAbscL   = 2
     GaussLimit = 2.7D0
     gfAbscissa = (/0D0, 0D0, 0.339981043584856D0, 0.861136311594053D0,              &
                    0.585786437627D0, 3.414213562373D0/)
     gfFTWeight = (/0D0, 0D0, 0.652145154862546D0, 0.347854845137454D0,              &
                     1.53332603312D0, 4.45095733505D0/)
  case (8)
     NumAbscG   = 6
     NumAbscL   = 2
     GaussLimit = 2.8D0
     gfAbscissa = (/0D0, 0D0, 0D0, 0.238619186083197D0, 0.661209386466265D0,          &
                    0.932469514203152D0, 0.585786437627D0, 3.414213562373D0/)
     gfFTWeight = (/0D0, 0D0, 0D0, 0.467913934572691D0, 0.360761573048139D0,          &
                   0.171324492379170D0, 1.53332603312D0, 4.45095733505D0/)
  case (10)
     NumAbscG   = 6
     NumAbscL   = 4
     GaussLimit = 2.6D0
     gfAbscissa = (/0D0, 0D0, 0D0, 0.238619186083197D0, 0.661209386466265D0,          &
                    0.932469514203152D0, 0.322547689619D0,1.745761101158D0,           &
                    4.536620296921D0,9.395070912301D0/)
     gfFTWeight = (/0D0, 0D0, 0D0, 0.467913934572691D0, 0.360761573048139D0,          &
                    0.171324492379170D0, 0.832739123838D0, 2.04810243845D0,           &
                    3.63114630582D0,6.48714508441D0/) 
  case (12)
     NumAbscG   = 8
     NumAbscL   = 4
     GaussLimit = 2.7D0
     gfAbscissa = (/0D0, 0D0, 0D0, 0D0, 0.183434642495650D0, 0.525532409916329D0,      &
                    0.796666477413627D0, 0.960289856497536D0, 0.322547689619D0,        &
                    1.745761101158D0, 4.536620296921D0, 9.395070912301D0/)
     gfFTWeight = (/0D0, 0D0, 0D0, 0D0, 0.362683783378362D0, 0.313706645877887D0,      &
                    0.222381034453374D0, 0.101228536290376D0, 0.832739123838D0,        &
                    2.04810243845D0, 3.63114630582D0,6.48714508441D0/)
  case (14)
     NumAbscG   = 8
     NumAbscL   = 6
     GaussLimit = 2.5D0
     gfAbscissa = (/0D0, 0D0, 0D0, 0D0, 0.183434642495650D0, 0.525532409916329D0,      &
                    0.796666477413627D0, 0.960289856497536D0, 0.222846604179D0,        &
                    1.188932101673D0, 2.992736326059D0, 5.775143569105D0,              & 
                    9.837467418383D0, 15.982873980602D0/)  
     gfFTWeight = (/0D0, 0D0, 0D0, 0D0, 0.362683783378362D0, 0.313706645877887D0,      &
                    0.222381034453374D0, 0.101228536290376D0, 0.573535507423D0,        &
                    1.36925259071D0, 2.26068459338D0, 3.35052458236D0,                 &
                    4.88682680021D0, 7.84901594560D0/)
  case (16)
     NumAbscG   = 10
     NumAbscL   = 6
     GaussLimit = 3.0D0
     gfAbscissa = (/0D0, 0D0, 0D0, 0D0, 0D0, 0.148874338981631D0, 0.433395394129247D0,  &
                    0.679409568299024D0, 0.865063366688985D0, 0.973906528517172D0,      &
                    0.222846604179D0, 1.188932101673D0, 2.992736326059D0,               &
                    5.775143569105D0, 9.837467418383D0, 15.982873980602D0/)
     gfFTWeight = (/0D0, 0D0, 0D0, 0D0, 0D0, 0.295524224714753D0, 0.269266719309996D0,  &
                    0.219086362515982D0, 0.149451349150581D0, 0.066671344308688D0,      &
                    0.573535507423D0, 1.36925259071D0, 2.26068459338D0,                 &
                    3.35052458236D0, 4.88682680021D0, 7.84901594560D0/) 
  end select

!--------------------------------------------------------------------------------
! Perform transformation X = KX**2 of the integration variable
! for Gauss-Legendre transform to treat log singularity.

! Original abscissas are symmetric in [-1,1], populate another half of 
! abscissa and weight arrays
  do i = 1, NumAbscG/2      
     gfAbscissa(i) = -gfAbscissa(NumAbscG-i+1) 
     gfFTWeight(i) =  gfFTWeight(NumAbscG-i+1)
  end do

! Carnahan, B., Luther, H.A., and Wilkes, J.O., 1969, 
! Applied Numerical Methods. p. 104-105
! Transformed abscissa: (Zi + 1) / 2
  do i = 1, NumAbscG
     gfAbscissa(i) = 0.5D0 * (gfAbscissa(i) + 1.0D0)
!     gfFTWeight(i) = 0.5D0 * gfFTWeight(i)
  end do

  tmp  = GaussLimit / MinTxRxSep    !  b * b
  InvPI = 1.0D0 / PI
  do i = 1, NumAbscG
     gfFTWeight(i) = tmp * InvPI * gfFTWeight(i) * gfAbscissa(i)
     gfAbscissa(i) = tmp * gfAbscissa(i)**2
  end do

!  Abscissas for Gauss-LaGuerre transform are shifted from (0, infinity)
!  to (GaussLimit, Infinity). The code below by Wannamaker, commented out,
!  seems incorrect.
!  tmp = 1.0 / MinTxRxSep
!  do i = 1, NumAbscL
!     j = i+NumAbscG
!     gfAbscissa(j) = (gfAbscissa(j) + GaussLimit) * tmp
!     gfFTWeight(j) = 2.0 * gfFTWeight(j) * tmp * InvPI 
!     gfFTWeight(j) = gfFTWeight(j) * tmp * InvPI 
!  end do

! Carnahan, B., Luther, H.A., and Wilkes, J.O., 1969, 
! Applied Numerical Methods. p. 115

  tmp = sqrt( GaussLimit / MinTxRxSep)
  do i = 1, NumAbscL
     j = i+NumAbscG
     gfAbscissa(j) = gfAbscissa(j) + tmp
     gfFTWeight(j) = gfFTWeight(j) * exp(-tmp) * InvPI 
  end do

  return
  end subroutine SetAbscissa





!================================================================
!  Based on Luo and Zhang, 1987

  subroutine SetAbscissa1(ForwAccuracy, MinTxRxSep, MaxTxRxSep)

  use GlobalForw
  implicit none

  integer, intent(in)     :: ForwAccuracy
  real(Rkind), intent(in) :: MinTxRxSep, MaxTxRxSep
  real(Rkind), parameter  :: K = 2.3D0
  integer :: i, N, nospace
  real(Rkind) :: MinAb, MaxAb

!------------------------------------------------------------------
! First, figure out number of abscissas to use

  MinAb = 0.3D0 / (MaxTxRxSep + 1.0D-8)
  MaxAb = 3.0D0 / (MinTxRxSep + 1.0D-8)
  N = int(1.5+log(MaxAb/MinAb) / log(K)) 

  if (ForwAccuracy == 0) then
     gNumAbscissa = N-3
  else if (ForwAccuracy == 2) then
     gNumAbscissa = N+3
  else
     gNumAbscissa = N
  end if

  if (gNumAbscissa < 5) gNumAbscissa = 5 
  if (gNumAbscissa > 32) gNumAbscissa = 32 

! Now initialize abscissas of inverse Fourier transform

  if (allocated(gfAbscissa)) deallocate(gfAbscissa)
  allocate(gfAbscissa(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  gfAbscissa = 0.0D0

  do i = 1, gNumAbscissa      
     gfAbscissa(i) = MinAb * (K ** (i-1))  
  end do

  return
  end subroutine SetAbscissa1

