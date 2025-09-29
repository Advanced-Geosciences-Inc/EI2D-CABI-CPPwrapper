
!  This file contains utility functions and subroutines for inversion
!--------------------------------------------------------------------------
! subroutine to write debug info to a file named debug.ech
! ???: This routine should be deleted later

  subroutine CheckPointInv(k1, k2)
  implicit none

  integer, intent(in) :: k1, k2
  integer :: iech

  iech = 18
  open(unit=iech, file='C:\cvs\Projects\Volente\debug.ech', status="unknown")
  write(iech, *) 'Stop at ', k1, k2
  close(iech)

  end subroutine CheckPoint

!------------------------------------------------------------------------------
! Jacobian times a vector

  subroutine JacobVec(InVec, Jacobian, OutVec)

!DEC$ ATTRIBUTES C, DLLEXPORT :: JacobVec
!DEC$ ATTRIBUTES ALIAS: "jacobvec" :: JacobVec

  use GlobalInv
  implicit none

  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumParam) :: InVec
  real(Rkind), intent(out), dimension(1:gNumData) :: OutVec

  integer :: i

  forall (i=1:gNumData) OutVec(i) =   &
          dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), InVec)

  return
  end subroutine JacobVec
!------------------------------------------------------------------------------
! Transposed Jacobian times a vector

  subroutine TJacobVec(InVec, Jacobian, OutVec)

!DEC$ ATTRIBUTES C, DLLEXPORT :: TJacobVec
!DEC$ ATTRIBUTES ALIAS: "tjacobvec" :: TJacobVec

  use GlobalInv
  implicit none

  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumData) :: InVec
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

  integer :: i, k

  k = (gNumData - 1) * gNumParam
  forall (i=1:gNumParam) OutVec(i) = dot_product(Jacobian(i:k+i:gNumParam), InVec)

  return
  end subroutine TJacobVec
!------------------------------------------------------------------------------
! Max row sum of Jt Wd J

  subroutine GetInitialLagrange(InVec, Jacobian, x)

  use GlobalInv
  implicit none

  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumData) :: InVec
  real(Rkind), intent(out) :: x

  real(Rkind), allocatable, dimension(:) :: tmp1, tmp2
  integer :: i, k

  allocate(tmp1(1:gNumParam), stat=i)
  tmp1 = 1.0D0
  allocate(tmp2(1:gNumData), stat=i)


  forall (i=1:gNumData) tmp2(i) =   &
          dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), tmp1)
  tmp2 = tmp2 * InVec

  k = (gNumData - 1) * gNumParam
  forall (i=1:gNumParam) tmp1(i) = dot_product(Jacobian(i:k+i:gNumParam), tmp2)
  x = MAXVAL(tmp1)

  if (allocated(tmp1)) deallocate(tmp1)
  if (allocated(tmp2)) deallocate(tmp2)

  return
  end subroutine GetInitialLagrange
!------------------------------------------------------------------------------
! Oldenburg and Li (1994), Inversion of IP data.
! Const parameters must agree with those in IPElem2IPParamPositivity in Delphi
  subroutine GetPositivityVec(InVec, PosVec)

  use GlobalInv
  implicit none

  real(Rkind), intent(in),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(out), dimension(1:gNumParam) :: PosVec

  if (gnIPPosMethod == 1) then   ! m = ln(eta)
     PosVec = exp(InVec)
  else
     PosVec = 1.0D0
  end if

  return
  end subroutine GetPositivityVec
!------------------------------------------------------------------------------
! Linear IP forward solution
  subroutine LinearIPForwSoln(VICalc, IPParam, Jacobian, IPCalc)

!DEC$ ATTRIBUTES C, DLLEXPORT :: LinearIPForwSoln
!DEC$ ATTRIBUTES ALIAS: "linearipforwsoln" :: LinearIPForwSoln

  use GlobalInv
  implicit none

  real(Rkind), intent(in), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in), dimension(1:gNumParam) :: IPParam
  real(Rkind), intent(in), dimension(1:gNumData)  :: VICalc
  real(Rkind), intent(out), dimension(1:gNumData) :: IPCalc

  integer :: i

  forall (i=1:gNumData) IPCalc(i) =   &
          dot_product(Jacobian((i-1)*gNumParam+1:i*gNumParam), IPParam)
  IPCalc = -IPCalc / (VICalc + gcTinyReal);

  return
  end subroutine LinearIPForwSoln
!------------------------------------------------------------------------------
