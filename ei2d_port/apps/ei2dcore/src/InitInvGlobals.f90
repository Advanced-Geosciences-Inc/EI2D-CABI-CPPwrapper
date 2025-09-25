!  InitInvGlobals.f90 
!
! Expose subroutine InitInvGlobals to users of this DLL
!
  subroutine InitInvGlobals(NumData, NumElemX, NumElemY, NumParamX, NumParamY, InvMethod,   &
              IPInvMethod, MaxNumIterInvCG, IPPosMeth, ModResoFactor, EpsilonD, EpsilonM)

!DEC$ ATTRIBUTES C, DLLEXPORT :: InitInvGlobals
!DEC$ ATTRIBUTES ALIAS: "initinvglobals" :: InitInvGlobals

  use GlobalInv
  implicit none

  integer, intent(in) :: NumData, NumElemX, NumElemY, NumParamX, NumParamY, MaxNumIterInvCG
  integer, intent(in) :: InvMethod, IPInvMethod, IPPosMeth
  real(Rkind), intent(in) :: ModResoFactor, EpsilonD, EpsilonM

! global variables: array dimensions
  gNumData           = NumData 
  gNumParamX         = NumParamX
  gNumParamY         = NumParamY
  gNumParam          = NumParamX * NumParamY
  gNumElemX          = NumElemX
  gNumElemY          = NumElemY
  gnInvMethod        = InvMethod
  gnIPInvMethod      = IPInvMethod
  gnMaxNumIterInvCG  = MaxNumIterInvCG
  gnIPPosMethod      = IPPosMeth
  gfModResoFactor    = ModResoFactor
  gfEpsilonD         = abs(EpsilonD)
  gfEpsilonM         = abs(EpsilonM)
  gnRobustMethod     = 2  

  return
  end subroutine InitInvGlobals

