!----------------------------------------------------------------------
! Update Jacobian using Broyden Quasi-Newton method
! Loke and Barker, 1996, Geophysical Prospecting, p. 135
! Tsourlos, Panagiotis, PhD dissertation, 1995. p. 249
!-----------------------------------------------------------------------

  subroutine UpdateJacobian(p, DataDiff, Jacobian) 

!DEC$ ATTRIBUTES C, DLLEXPORT :: UpdateJacobian
!DEC$ ATTRIBUTES ALIAS: "updatejacobian" :: UpdateJacobian

  use GlobalInv
  implicit none

! variables and arrays passed from Delphi GUI

  real(Rkind), intent(inout), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumParam) :: p       ! Model update
  real(Rkind), intent(in), dimension(1:gNumData)  :: DataDiff

! Local Variable
  integer :: i, j, k
  real(Rkind) :: pNorm
  real(Rkind), dimension(1:gNumData) :: DataVec

  pNorm = dot_product(p, p) 

  ! Jacobian times modelupdate vector 
  forall (i=1:gNumData) DataVec(i) =   &
            dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), p)

  ! u vector (in data space) in Loke and Barker (1996) p. 135
  pNorm = 1.0D0 / pNorm 
  do i = 1, gNumData
     DataVec(i) = (DataDiff(i) - DataVec(i)) * pNorm
  end do 

  ! Assmeble new Jacobian
  do i = 1, gNumData
     do j = 1, gNumParam
        k = (i-1) * gNumParam + j 
        Jacobian(k) = Jacobian(k) + DataVec(i) * p(j) 
     end do
  end do

  return
  end subroutine UpdateJacobian
