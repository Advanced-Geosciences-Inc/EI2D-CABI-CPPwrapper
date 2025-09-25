
! FindAbscissa sets number of absicssas and the actual abscissas
! using the algorithm in Luo and Zhang (1987) and Queralt et al (1991)
! This algorithm is based on logarithmic approximation for small 
! wavenumber and piecewise exponential approximation for larger ones
!----------------------------------------------------------------------

  subroutine FindAbscissa(ForwAccuracy, MinTxRxSep, MaxTxRxSep)

  use GlobalForw, only : Rkind,gNumAbscissa,gfAbscissa,gfFTWeight
  implicit none

  integer, intent(in)     :: ForwAccuracy
  real(Rkind), intent(in) :: MinTxRxSep, MaxTxRxSep
  real(Rkind), parameter  :: Kfactor = 2.5
  
  real(Rkind) :: Kmin, Kmax
  integer :: i, NoSpace

  Kmin = 0.1D0 / MaxTxRxSep
  Kmax = 3.0D0 / MinTxRxSep

  gNumAbscissa = nint(log(Kmax/Kmin) / log(Kfactor) + 1)
  if (gNumAbscissa < 5)  gNumAbscissa = 5
  if (gNumAbscissa > 16) gNumAbscissa = 16

  if (allocated(gfAbscissa)) deallocate(gfAbscissa)
  allocate(gfAbscissa(1:gNumAbscissa), stat=nospace)
  if (NoSpace /= 0) return                     ! ???: err message desired
  gfAbscissa = 0.0D0

  do i = 1, gNumAbscissa
     gfAbscissa(i) = Kmin * (Kfactor ** (i-1))
  end do

  return
  end subroutine FindAbscissa
!-------------------------------------------------------------------------
! Program for selecting k and g in Fourier inverse transformation
! FindAbscissa2 below uses Xu Shizhe's optimization algorithm but
! it is not stable so far (12/21/2001, Yang) 

  subroutine FindAbscissa2(ForwAccuracy, MinTxRxSep, MaxTxRxSep)

  use GlobalForw, only : Rkind,gNumAbscissa,gfAbscissa,gfFTWeight
  implicit none

  integer, intent(in)     :: ForwAccuracy
  real(Rkind), intent(in) :: MinTxRxSep, MaxTxRxSep
  real(Rkind), parameter  :: PI = 3.14159265359D0
  integer, parameter :: NR = 41

  real(Rkind), dimension(1:NR) :: u, u1, c, d, r
  real(Rkind), allocatable :: p(:), rhs(:),  a(:,:), a1(:,:), uk(:,:), b(:,:)

  real(Rkind) :: x, BesselK0, fi, fiOld, MinR, MaxR, dr, Damping
  integer :: i, j, NumIter, iter, nospace, ReverseCount 

  if (ForwAccuracy == 0) then
     gNumAbscissa = 5
  else if (ForwAccuracy == 1) then
     gNumAbscissa = 7
  else if (ForwAccuracy == 2) then
     gNumAbscissa = 9
  else
     gNumAbscissa = 7
  end if
  gNumAbscissa = 8

  if (allocated(gfAbscissa)) deallocate(gfAbscissa)
  allocate(gfAbscissa(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  gfAbscissa = 0.0D0

  if (allocated(gfFTWeight)) deallocate(gfFTWeight)
  allocate(gfFTWeight(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  gfFTWeight = 0.0D0

  allocate(rhs(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  rhs = 0.0D0

  allocate(p(1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  p = 0.0D0

  allocate(a(1:NR, 1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  a = 0.0D0

  allocate(a1(1:NR, 1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  a1 = 0.0D0

  allocate(uk(1:NR, 1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  uk = 0.0D0

  allocate(b(1:gNumAbscissa, 1:gNumAbscissa), stat=nospace)
  if (nospace /= 0) return                     ! ???: err message desired
  b = 0.0D0

  select case (gNumAbscissa)
  case (5)
     gfAbscissa = (/0.004758, 0.0407011, 0.1408855, 0.393225, 1.08838/)             
  case (7)
     gfAbscissa = (/0.000694, 0.0007186, 0.0037328, 0.01765, 0.0822884, 0.40969, 2.179599/)          
  case (8)
     gfAbscissa =(/0.00001, 0.00004, 0.0002, 0.001, 0.006, 0.049, 0.3, 2.0/)
  case (9)
     gfAbscissa = (/0.00001, 0.00008, 0.0005, 0.003, 0.02, 0.1, 0.4, 2.0, 12.0/)         
  end select

  ! Find NR's r from Min and Max TX-RX separation.
  if (MinTxRxSep < 0.01) then
     MinR = 0.01D0
  else
     MinR = 0.8D0 * MinTxRxSep
  end if

  if (MaxTxRxSep > 1000.0) then
     MaxR = 1000.0D0
  else
     MaxR = 1.2D0 * MaxTxRxSep
  end if

  MinR = 0.1D0
  MaxR = 10000.0D0
  dr = log(MaxR/MinR) / (NR - 1)
  Damping = log(MinR)              ! tmp variable
  do i = 1, NR
     r(i) = exp(Damping + dr * (i-1))
  end do

  Damping = 1.0D-03
  ReverseCount = 0
  NumIter = 60
  do iter = 0, NumIter

     ! assemble the array a(m,n) and a1(m,n)
     do i = 1, NR
        c(i) = 1.0D0
        do j = 1, gNumAbscissa
           x = r(i) * gfAbscissa(j)
           a(i,j)  = r(i) * BesselK0(x)
           a1(i,j) = r(i) * BesselK0(1.1D0*x)
        end do
     end do

! to get g(n) under condition that k is given
! using the least square method
     ! call acbg(a,c,b,rhs,NR,gNumAbscissa)
     b = matmul(transpose(a), a)
     rhs = matmul(transpose(a), c)
     do i = 1, gNumAbscissa
        b(i,i) = b(i,i) + damping
     end do

! solve the linear equation system
     p = 0.
     call CholeskyDC(b, gNumAbscissa, p)
     call ChleskySolve(b, gNumAbscissa, p, rhs, gfFTWeight)

    ! to get u from a and gfFTWeight
     u = matmul(a, gfFTWeight) 

     ! calculate the error
     fiOld = fi
     fi=0.0D0
     do i=1,NR
        fi=fi+(1D0-u(i))**2
     end do
     fi = sqrt(fi/NR)
     if(iter == NumIter) exit
     if ((iter > 1) .and. (fi > fiOld)) ReverseCount = ReverseCount + 1
     if (ReverseCount > 1) exit

     ! construct array uk(m,n)
     do j=1,gNumAbscissa

        do i=1,NR
           d(i)=a(i,j)
           a(i,j)=a1(i,j)
        end do

        ! call acbg(a,c,b,rhs,NR,gNumAbscissa)
        b = matmul(transpose(a), a)
        rhs = matmul(transpose(a), c)
        do i = 1, gNumAbscissa
           b(i,i) = b(i,i) + damping
        end do
        p = 0.
        call CholeskyDC(b, gNumAbscissa, p)
        call ChleskySolve(b, gNumAbscissa, p, rhs, gfFTWeight)
        u1 = matmul(a, gfFTWeight) 

        do i=1,NR
            a(i,j)=d(i)
        end do

        do i=1,NR
           uk(i,j) = (u1(i)-u(i)) / (0.10D0*gfAbscissa(j))
        end do

     end do

     ! to get delta k
     do i = 1,NR
        c(i)=1-u(i)
     end do
     ! call acbg(uk,c,b,rhs,NR,gNumAbscissa)
     b = matmul(transpose(uk), a)
     rhs = matmul(transpose(uk), c)
     do i = 1, gNumAbscissa
        b(i,i) = b(i,i) + damping
     end do
     p = 0.0D0
     call CholeskyDC(b, gNumAbscissa, p)
     call ChleskySolve(b, gNumAbscissa, p, rhs, gfFTWeight)   

     ! delta k have been obtained, then add delta k to k
     do i = 1, gNumAbscissa
        gfAbscissa(i) = gfAbscissa(i) + 0.20D0*gfFTWeight(i)
        if(gfAbscissa(i) <= 0.) gfAbscissa(i) = gfAbscissa(i) - 0.20D0*gfFTWeight(i)
     end do
     Damping = Damping * 0.7D0

  end do

!  gfAbscissa = (/0.000694, 0.0007186, 0.0037328, 0.01765, 0.0822884, 0.40969, 2.179599, 10.89892/)         
!  gfFTWeight = (/0.0001491, 0.0008157, 0.0037312, 0.0173954, 0.0807139, 0.4407338, 2.242198, 11.42917/)         

  gfFTWeight = 2.0D0 * gfFTWeight / PI
!  gfAbscissa = gfAbscissa / MinTxRxSep


  if (allocated(rhs)) deallocate(rhs)
  if (allocated(p)) deallocate(p)
  if (allocated(a)) deallocate(a)
  if (allocated(a1)) deallocate(a1)
  if (allocated(uk)) deallocate(uk)
  if (allocated(b)) deallocate(b)

  return
  end subroutine FindAbscissa2

!---------------------------------------------------------
  subroutine CholeskyDC(a, n, p)
  use GlobalForw, only : Rkind
  implicit none

  integer, intent(in) :: n
  real(Rkind), intent(inout) :: a(n, n), p(n)
  integer :: i, j, k
  real(Rkind) :: sum

  do i = 1, n
     do j = i, n
        sum = a(i, j)
        do k = i-1, 1, -1
           sum = sum - a(i,k) * a(j,k)
        end do
        if (i == j) then
           ! if (sum <= 0.) pause 'Cholesky decomposition in FindAbsc failed'
           p(i) = sqrt(sum)
        else
           a(j,i) = sum / p(i)
        end if
     end do
  end do

  return
  end subroutine CholeskyDC

!---------------------------------------------------------
  subroutine ChleskySolve(a, n, p, b, x)
  use GlobalForw, only : Rkind
  implicit none

  integer, intent(in) :: n
  real(Rkind), intent(in) :: a(1:n, 1:n), p(n), b(n)
  real(Rkind), intent(out) :: x(n)
  real(Rkind) :: sum
  integer :: i, k

  do i = 1, n         ! Forward substitution
     sum = b(i)
     do k = i-1, 1, -1
        sum = sum - a(i,k) * x(k)
     end do
     x(i) = sum / p(i)
  end do

  do i = n, 1, -1     ! Backward substitution
     sum = x(i)
     do k = i+1, n
        sum = sum - a(k,i) * x(k)
     end do
     x(i) = sum / p(i)
  end do

  return
  end subroutine ChleskySolve

!----------- END ----------------------------------



