! =====================================================================
! ei2d_c_api.f90  —  C ABI wrappers for EI2D forward model (FD path)
! =====================================================================

subroutine ei2d_cb_stub(status)
  implicit none
  integer, intent(inout) :: status(0:1)
  status = 0
end subroutine ei2d_cb_stub

! ---------------------- C ABI: InitForwGlobals -----------------------
subroutine ei2d_InitForwGlobals(NumData, NumElectrodes, NumInfElectrodes, NumNodeX, NumNodeY, &
                                ForwModMeth, ForwSolver, InvMethod,                         &
                                ForwAccuracy, ForwCGIter, BCType,                           &
                                ForwCGResid, MinTxRxSep, MaxTxRxSep) bind(C, name='ei2d_InitForwGlobals')
  use iso_c_binding
  implicit none
  ! C ABI scalars
  integer(c_int),  value :: NumData, NumElectrodes, NumInfElectrodes, NumNodeX, NumNodeY
  integer(c_int),  value :: ForwModMeth, ForwSolver, InvMethod, ForwAccuracy, ForwCGIter, BCType
  real(c_double),  value :: ForwCGResid, MinTxRxSep, MaxTxRxSep
  ! engine interface (implicit okay)
  interface
    subroutine InitForwGlobals(NumData, NumElectrodes, NumInfElectrodes, NumNodeX, NumNodeY, &
                               ForwModMeth, ForwSolver, InvMethod, ForwAccuracy, ForwCGIter, &
                               BCType, ForwCGResid, MinTxRxSep, MaxTxRxSep)
      implicit none
      integer, intent(in) :: NumData, NumElectrodes, NumInfElectrodes, NumNodeX, NumNodeY
      integer, intent(in) :: ForwModMeth, ForwSolver, InvMethod, ForwAccuracy, ForwCGIter, BCType
      real   , intent(in) :: ForwCGResid, MinTxRxSep, MaxTxRxSep
    end subroutine
  end interface
  ! local kind-converted reals (engine uses real(Rkind) ~ double)
  real :: r_ForwCGResid, r_MinTxRxSep, r_MaxTxRxSep

  r_ForwCGResid = real(ForwCGResid, kind(r_ForwCGResid))
  r_MinTxRxSep  = real(MinTxRxSep , kind(r_MinTxRxSep ))
  r_MaxTxRxSep  = real(MaxTxRxSep , kind(r_MaxTxRxSep ))

  call InitForwGlobals(NumData, NumElectrodes, NumInfElectrodes, NumNodeX, NumNodeY, &
                       ForwModMeth, ForwSolver, InvMethod, ForwAccuracy, ForwCGIter, &
                       BCType, r_ForwCGResid, r_MinTxRxSep, r_MaxTxRxSep)
end subroutine ei2d_InitForwGlobals

! ---------------------- C ABI: SetNumParamForward --------------------
subroutine ei2d_SetNumParamForward(NumParamX, NumParamY) bind(C, name='ei2d_SetNumParamForward')
  use iso_c_binding
  implicit none
  integer(c_int), value :: NumParamX, NumParamY
  interface
    subroutine SetNumParamForward(NumParamX, NumParamY)
      implicit none
      integer, intent(in) :: NumParamX, NumParamY
    end subroutine
  end interface
  call SetNumParamForward(NumParamX, NumParamY)
end subroutine ei2d_SetNumParamForward

! ---------------------- C ABI: ForwardFD -----------------------------
subroutine ei2d_ForwardFD(NodeX, NodeY, Cond, VIcalc, Jacobian,                          &
                          ElecNodeID, StingCMD, ParamX1, ParamX2, ParamY1, ParamY2,       &
                          InfElec, GetJacobian, nNodes, nElem, nData) bind(C, name='ei2d_ForwardFD')
  use iso_c_binding
  use GlobalForw, only: gNumElectrodes, gNumParam, gNumParamX, gNumParamY, gNumNodeX, gNumNodeY, gNumInfElectrodes
  implicit none
  ! ------------ C ABI dummies ------------
  real(c_double),  intent(in)    :: NodeX(*), NodeY(*), Cond(*)
  real(c_double),  intent(inout) :: VIcalc(*), Jacobian(*)
  integer(c_int),  intent(in)    :: ElecNodeID(*), StingCMD(*), ParamX1(*), ParamX2(*), ParamY1(*), ParamY2(*), InfElec(*)
  integer(c_int),  value         :: GetJacobian, nNodes, nElem, nData

  ! ------------ engine explicit interface (only ForwardFD needs it because of procedure-arg) ------------
  interface
    subroutine ForwardFD(CallBackForw, NodeX, NodeY, Cond, VIcalc, Jacobian, ElecNodeID, &
                         StingCMD, ParamX1, ParamX2, ParamY1, ParamY2, InfElec, GetJacobian)
      use GlobalForw
      implicit none
      interface
        subroutine CallBackForw(CallBackStatus)
          implicit none
          integer, intent(inout), dimension(0:1) :: CallBackStatus
        end subroutine
      end interface
      real   , intent(in),  dimension(1:gNumNodes) :: NodeX, NodeY
      real   , intent(in),  dimension(1:gNumElem)  :: Cond
      real   , intent(out), dimension(1:gNumData)  :: VIcalc
      real   , intent(out), dimension(1:gNumData*gNumParam) :: Jacobian
      integer, intent(in),  dimension(1:gNumElectrodes)     :: ElecNodeID
      integer, intent(in),  dimension(1:gNumData*4)         :: StingCMD
      integer, intent(in),  dimension(1:gNumParamX) :: ParamX1, ParamX2
      integer, intent(in),  dimension(1:gNumParamY) :: ParamY1, ParamY2
      integer, intent(in),  dimension(1:gNumInfElectrodes) :: InfElec
      integer, intent(in) :: GetJacobian
    end subroutine
  end interface
  external :: ei2d_cb_stub

  ! ------------ locals (engine-kind working arrays) ------------
  integer :: i, oldNe, needJac, nparam, safeElec, maxCmd
  real    , allocatable :: x(:), y(:), c(:), vi(:), jac(:)
  integer , allocatable :: en(:), sc(:), px1(:), px2(:), py1(:), py2(:), inf(:)

  ! sizes from module (already set via InitForwGlobals / SetNumParamForward)
  nparam = max(1, gNumParam)

  ! allocate & copy primary arrays
  allocate(x(nNodes), y(nNodes), c(nElem), vi(nData))
  do i = 1, nNodes
    x(i) = real(NodeX(i), kind(x))
    y(i) = real(NodeY(i), kind(y))
  end do
  do i = 1, nElem
    c(i) = real(Cond(i), kind(c))
  end do
  vi = 0.0

  ! Jacobian workspace (always allocate so callee has a valid target; copy back only if requested)
  allocate(jac(nData*nparam))
  jac = 0.0

  ! copy StingCMD; compute highest referenced electrode index
  allocate(sc(nData*4))
  do i = 1, nData*4
    sc(i) = StingCMD(i)
  end do
  maxCmd = 0
  if (nData > 0) maxCmd = maxval(sc(1:nData*4))

  ! Expand gNumElectrodes if commands reference a higher index
  oldNe = gNumElectrodes
  if (maxCmd > gNumElectrodes) gNumElectrodes = maxCmd

  ! mirror ElecNodeID with the CURRENT module bound
  safeElec = gNumElectrodes
  allocate(en(safeElec))
  ! copy what we have
  do i = 1, min(oldNe, safeElec)
    en(i) = ElecNodeID(i)
  end do
  ! pad tail if expanded
  if (safeElec > oldNe) then
    do i = oldNe+1, safeElec
      en(i) = en(oldNe)
    end do
  end if

  ! parameter windows across full grid extents
  allocate(px1(gNumParamX), px2(gNumParamX))
  allocate(py1(gNumParamY), py2(gNumParamY))
  do i = 1, gNumParamX
    px1(i) = 1
    px2(i) = gNumNodeX
  end do
  do i = 1, gNumParamY
    py1(i) = 1
    py2(i) = gNumNodeY
  end do

  ! infinite electrodes (if none, pass a single zero to satisfy dummy shape)
  if (gNumInfElectrodes > 0) then
    allocate(inf(gNumInfElectrodes))
    do i = 1, gNumInfElectrodes
      inf(i) = InfElec(i)
    end do
  else
    allocate(inf(1))
    inf(1) = 0
  end if

  ! need Jacobian?
  needJac = 0
  if (GetJacobian /= 0) needJac = 1

  ! ---------------- call engine ----------------
  call ForwardFD(ei2d_cb_stub, x, y, c, vi, jac, en, sc, px1, px2, py1, py2, inf, needJac)

  ! ---------------- copy back ------------------
  do i = 1, nData
    VIcalc(i) = vi(i)
  end do
  if (needJac /= 0 .and. nparam > 0) then
    do i = 1, nData*nparam
      Jacobian(i) = jac(i)
    end do
  end if

  ! (no deallocation required; process ends—kept simple for debugging)
end subroutine ei2d_ForwardFD
