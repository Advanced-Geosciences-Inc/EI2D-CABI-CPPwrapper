!  InitForwGlobals.f90 
!
! Expose subroutine InitForwGlobals to users of this DLL
!
  subroutine InitForwGlobals(NumData, NumElectrodes, NumInfElectrodes, NumNodeX,  &
                             NumNodeY, ForwModMeth, ForwSolver, InvMethod,        &
                             ForwAccuracy, ForwCGIter, BCType, ForwCGResid,       &
                             MinTxRxSep, MaxTxRxSep)

!DEC$ ATTRIBUTES C, DLLEXPORT :: InitForwGlobals
!DEC$ ATTRIBUTES ALIAS: "initforwglobals" :: InitForwGlobals

  use GlobalForw
  implicit none

  integer, intent(in) :: NumData, NumElectrodes, NumNodeX, NumNodeY
  integer, intent(in) :: ForwModMeth, ForwSolver, InvMethod, ForwAccuracy
  integer, intent(in) :: ForwCGIter, BCType, NumInfElectrodes
  real(Rkind), intent(in) :: ForwCGResid, MinTxRxSep, MaxTxRxSep

! global variables: array dimensions
  gNumData          = NumData 
  gNumElectrodes    = NumElectrodes
  gNumInfElectrodes = NumInfElectrodes
  gNumNodeX         = NumNodeX
  gNumNodeY         = NumNodeY
  gnForwModMeth     = ForwModMeth
  gnForwSolver      = ForwSolver
  gnInvMethod       = InvMethod
  
  gNumNodes         = gNumNodeX * gNumNodeY
  gNumElemX         = gNumNodeX - 1
  gNumElemY         = gNumNodeY - 1
  gNumElem          = gNumElemX * gNumElemY
  gnForwCGIter      = ForwCGIter
  gnBCType          = BCType
  gfForwCGResid     = ForwCGResid
  
  if (gnForwModMeth == 0) then          !  Finite difference method
     gNumDiagonal = 3
     gnBand       = gNumNodeY + 1
  else if (gnForwModMeth == 1) then     !  Finite element method
     gNumDiagonal = 5
     gnBand       = gNumNodeY + 2
  end if

  call SetAbscissa(ForwAccuracy, MinTxRxSep)
!  call SetAbscissa1(ForwAccuracy, MinTxRxSep, MaxTxRxSep)
!  call FindAbscissa(ForwAccuracy, MinTxRxSep, MaxTxRxSep)

  return
  end subroutine InitForwGlobals

! -------------------------------------------------------------------------------
  subroutine SetNumParamForward(NumParamX, NumParamY)

!DEC$ ATTRIBUTES C, DLLEXPORT :: SetNumParamForward
!DEC$ ATTRIBUTES ALIAS: "setnumparamforward" :: SetNumParamForward

  use GlobalForw
  implicit none

  integer, intent(in) :: NumParamX, NumParamY

  gNumParamX = NumParamX
  gNumParamY = NumParamY
  gNumParam = NumParamX * NumParamY

  return
  end subroutine SetNumParamForward
