
!  These subroutines calculate several combinations of 
!  Roughness times a vector
!  Note that R = {... -1 ... +1 ... }
! 
!  RoughV0:    OutVec = RtR * InVec
!  RoughV1:    OutVec = R * InVec
!  RoughV2:    OutVec = RtR * InVec, InVec = diag{RtR}   
!  RoughV3:    OutVec = RtR * InVec, OutVec2 = R * InVec, InVec = diag{RtR}
!  RoughV3 is the same as RoughV2   6/27/07

!  For clarity, there are some redundancy in these routines. I do not want 
!  to live with too many IFs and flags
!------------------------------------------------------------------------------

! RoughVec0: OutVec = RtR * InVec
  subroutine RoughV0(InVec, RoughX, RoughY, OutVec)

!DEC$ ATTRIBUTES C, DLLEXPORT :: RoughV0
!DEC$ ATTRIBUTES ALIAS: "roughv0" :: RoughV0

  use GlobalInv
  implicit none

  real(Rkind), intent(in),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

  integer :: i, p1, ipx, ipy
  real(Rkind) :: tmp


! Start with the X-direction roughness. 
! RoughX * RoughX may be redundant. However, the power factor over 
! data weights may be removed if we use R * R instead of RoughX alone.

  OutVec = 0D0
  do ipx = 1, gNumParamX-1
     do ipy = 1, gNumParamY
        i = (ipy - 1) * gNumParamX + ipx
        p1 = i + 1                            ! next param
        tmp = RoughX(i) * RoughX(i) * (InVec(i) - InVec(p1))
        OutVec(i)  = OutVec(i)  + tmp
        OutVec(p1) = OutVec(p1) - tmp
     end do
  end do


! Y-direction Roughness * Vector. This case is more complicated 
! since each layer has different number of parameters.

  do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY-1
        i = (ipy - 1) * gNumParamX + ipx
        p1 = i + gNumParamX                   ! next param
        tmp = RoughY(i) * RoughY(i) * (InVec(i) - InVec(p1))
        OutVec(i)  = OutVec(i)  + tmp
        OutVec(p1) = OutVec(p1) - tmp
     end do
  end do

  return
  end subroutine RoughV0


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  RoughV1:    OutVec = R * InVec, for Robust Inversion only

  subroutine RoughV1(InVec, RoughX, RoughY, OutVec)

!DEC$ ATTRIBUTES C, DLLEXPORT :: RoughV1
!DEC$ ATTRIBUTES ALIAS: "roughv1" :: RoughV1

  use GlobalInv
  implicit none

  real(Rkind), intent(in),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

  integer :: i, ipx, ipy


! X direction
  OutVec = 0D0
  do ipx = 1, gNumParamX-1
     do ipy = 1, gNumParamY
        i = (ipy - 1) * gNumParamX + ipx
        OutVec(i) = RoughX(i) * (InVec(i+1) - InVec(i))
	 end do
  end do


! Y direction
  do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY-1
        i = (ipy - 1) * gNumParamX + ipx
        OutVec(i) = RoughY(i) * (InVec(i+gNumParamX) - InVec(i))
	 end do
  end do
 
  return
  end subroutine RoughV1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  RoughV2:    OutVec = RtR * InVec, InVec = diag{RtR}

  subroutine RoughV2(InVec, RoughX, RoughY, OutVec)

  use GlobalInv
  implicit none

  real(Rkind), intent(inout),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec
! Local Variables

  integer :: i, p1, ipx, ipy
  real(Rkind) :: R2, tmp
  real(Rkind), allocatable, dimension(:) :: DiagRtR


  allocate(DiagRtR(1:gNumParam), stat=i)
  DiagRtR = 0D0

! Start with the X-direction roughness. 
! RoughX * RoughX may be redundant. However, the power factor over 
! data weights may be removed if we use R * R instead of RoughX alone.

  DiagRtR = 0D0
  OutVec = 0D0

  do ipx = 1, gNumParamX-1
     do ipy = 1, gNumParamY
        i = (ipy - 1) * gNumParamX + ipx
        R2 = RoughX(i) * RoughX(i)
        p1 = i + 1                              ! next param
        tmp = R2 * (InVec(i) - InVec(p1))
        OutVec(i)  = OutVec(i)  + tmp
        OutVec(p1) = OutVec(p1) - tmp
        DiagRtR(i)  = DiagRtR(i) + R2 
        DiagRtR(p1) = DiagRtR(p1) + R2
     end do
  end do

! Y-direction Roughness * Vector. This case is more complicated 
! since each layer has different number of parameters.

  do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY-1
        i = (ipy - 1) * gNumParamX + ipx
        R2 = RoughY(i) * RoughY(i)
        p1 = i + gNumParamX                  ! next param
        tmp = R2 * (InVec(i) - InVec(p1))
        OutVec(i)  = OutVec(i)  + tmp
        OutVec(p1) = OutVec(p1) - tmp
        DiagRtR(i)  = DiagRtR(i) + R2 
        DiagRtR(p1) = DiagRtR(p1) + R2
     end do
  end do

  InVec = DiagRtR
  if (allocated(DiagRtR)) deallocate(DiagRtR)

  return
  end subroutine RoughV2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  RoughV3:    OutVec = RtR * InVec, InVec = diag{RtR}
  subroutine RoughV3(InVec, RoughX, RoughY, OutVec)

  use GlobalInv
  implicit none

  real(Rkind), intent(inout),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: RoughX, RoughY
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

! Local Variables

  integer :: i, p1, ipx, ipy
  real(Rkind) :: R2, dif, tmp
  real(Rkind), allocatable, dimension(:) :: DiagRtR


  
! InVec will carry the result of diagRtR back, but we should not 
! mesh up InVec from the beginning
  allocate(DiagRtR(1:gNumParam), stat=i)
  DiagRtR = 0D0

! Start with the X-direction roughness. 
! RoughX * RoughX may be redundant. However, the power factor over 
! data weights may be removed if we use R * R instead of RoughX alone.

  DiagRtR = 0D0
  OutVec  = 0D0

  do ipx = 1, gNumParamX-1
     do ipy = 1, gNumParamY
        i = (ipy - 1) * gNumParamX + ipx
        R2 = RoughX(i) * RoughX(i)
        p1 = i + 1                       ! next param
		dif = InVec(p1) - InVec(i)             
        tmp = R2 * dif
        OutVec(i)  = OutVec(i)  - tmp
        OutVec(p1) = OutVec(p1) + tmp

        DiagRtR(i)  = DiagRtR(i) + R2 
        DiagRtR(p1) = DiagRtR(p1) + R2
     end do
  end do


! Y-direction Roughness * Vector. This case is more complicated 
! since each layer has different number of parameters.

  do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY-1
        i = (ipy - 1) * gNumParamX + ipx
        R2 = RoughY(i) * RoughY(i)
        p1 = i + gNumParamX                       ! next param
		dif = InVec(p1) - InVec(i)             
        tmp = R2 * dif
        OutVec(i)  = OutVec(i)  - tmp
        OutVec(p1) = OutVec(p1) + tmp

        DiagRtR(i)  = DiagRtR(i) + R2 
        DiagRtR(p1) = DiagRtR(p1) + R2
     end do
  end do

  InVec = DiagRtR

  if (allocated(DiagRtR)) deallocate(DiagRtR)

  return
  end subroutine RoughV3


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  RoughL1X:    OutVec = R * InVec, for Robust Inversion only
  subroutine RoughL1X(InVec, R, OutVec)
  use GlobalInv
  implicit none

  real(Rkind), intent(in),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: R
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

  integer :: i, ipx, ipy

! X direction
  OutVec = 0D0
  do ipx = 1, gNumParamX-1
     do ipy = 1, gNumParamY
        i = (ipy - 1) * gNumParamX + ipx
        OutVec(i) = R(i) * (InVec(i+1) - InVec(i))
	 end do
  end do
 
  return
  end subroutine RoughL1X

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  RoughL1Y:    OutVec = R * InVec, for Robust Inversion only
  subroutine RoughL1Y(InVec, R, OutVec)
  use GlobalInv
  implicit none

  real(Rkind), intent(in),  dimension(1:gNumParam) :: InVec
  real(Rkind), intent(in),  dimension(1:gNumParam) :: R
  real(Rkind), intent(out), dimension(1:gNumParam) :: OutVec

  integer :: i, ipx, ipy

! Y direction
  OutVec = 0D0
  do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY-1
        i = (ipy - 1) * gNumParamX + ipx
        OutVec(i) = R(i) * (InVec(i+gNumParamX) - InVec(i))
	 end do
  end do
 
  return
  end subroutine RoughL1Y

!=========  end of this file ==================

