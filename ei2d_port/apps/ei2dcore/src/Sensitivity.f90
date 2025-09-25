
!  Calculate Jacobian matrix for resistivity inversion. Jacobian will be 
!  also used in IP forward modeling and inversion.
!  This file contains two subroutines: JacobFE and JacobFD

! Estimate Jacobian (Sensitivity) Matrix using electrical reciprocity.
! The source and receiver poles are interchangeable in the sense that 
! the voltage at a receiver due to a source point of strength s equals 
! the voltage that would be measured at the source point if a current of 
! strength s was placed at the receiver. Wannamaker, 1992 (pages 29-31).
! 
! The sensitivity is defined as derivative of V/I (data) with respect 
! to logarithm of electrical conductivities.

!  Old Version Declaration:

!  subroutine JacobFE(Jacobian, NodeX, NodeY, Cond, StingCMD, ParamX1, ParamX2,  &
!                     ParamY1, ParamY2, abscix, weight, CenterNodeX,             &
!                     CenterNodeY, ElemArea, ElemStiffAll)
!-------------------------------------------------------------------------------- 
  subroutine JacobFE(Jacobian, StingCMD, ParamX1, ParamX2, ParamY1, ParamY2,  &
                      weight, ElemStiffAll)

  use GlobalForw
  implicit none

  real(Rkind), intent(inout), dimension(1:gNumData*gNumParam) :: Jacobian
  integer, intent(in), dimension(1:gNumData*4)     :: StingCMD
  integer, intent(in), dimension(1:gNumParamX) :: ParamX1,ParamX2
  integer, intent(in), dimension(1:gNumParamY) :: ParamY1,ParamY2
  real(Rkind), intent(inout), dimension(1:gNumElem*16) :: ElemStiffAll
  real(Rkind), intent(in) ::  weight

  integer :: iParam, iData, iElemX, iElemY, iElem, iJacob0, iJacob
  integer :: i, j, k, iPos, A, B, M, N, ipx, ipy
  integer ::  LocalNodes(4)
  real(Rkind), allocatable, dimension(:, :) :: Vab, Vmn  
!----------------------------------------------------------------

  if (allocated(Vab)) deallocate(Vab)
  allocate(Vab(gNumData, gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Vmn)) deallocate(Vmn)
  allocate(Vmn(gNumData, gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

! Assuming a current electrode at A, B, M, N of iData. Below are the 
! beginning array index of corresponding solution in gfVoverI

  do iData = 1, gNumData               
     k = (iData-1) * 4                 
     A = StingCMD(k+1) * gNumNodes   
     B = StingCMD(k+2) * gNumNodes 
     M = StingCMD(k+3) * gNumNodes   
     N = StingCMD(k+4) * gNumNodes
     do j = 1, gNumNodes
        Vab(iData, j) = gfVoverI(A+j) - gfVoverI(B+j)
        Vmn(iData, j) = gfVoverI(M+j) - gfVoverI(N+j)
     end do
  end do
    
  ElemStiffAll = ElemStiffAll * weight


  do ipx = 1, gNumParamX
  do ipy = 1, gNumParamY
     iParam = (ipy - 1) * gNumParamX + ipx
     iJacob0 = iParam - gNumParam

!    Unscramble the parameter lumping scheme for element-based stiffness
!    ParamX1,ParamX2,ParamY1,ParamY2 are Elemental bounds of current parameter 
     do iElemY = ParamY1(ipy), ParamY2(ipy)
        do iElemX = ParamX1(ipx), ParamX2(ipx)

           iElem = iElemY + gNumElemY * (iElemX - 1)
           LocalNodes(1) = (iElemX-1) * gNumNodeY + iElemY
           LocalNodes(2) = LocalNodes(1) + 1
           LocalNodes(3) = LocalNodes(1) + gNumNodeY
           LocalNodes(4) = LocalNodes(2) + gNumNodeY

           iPos = 16 * (iElem - 1) - 4   

           ! Get local stiffness for current quadrilateral element (iElem)
           ! call StiffnessFE(NodeX, NodeY, Cond, ElemStiff, abscix, iElemX, &
           !            iElemY, LocalNodes, CenterNodeX, CenterNodeY, ElemArea)

           ! Loop through Reduced Stiffness matrix of current elements
           do iData = 1, gNumData
              iJacob = iJacob0 + gNumParam * iData
              do j = 1, 4
                 k = iPos + 4 * j
                 do i = 1, 4
                     k = k + 1
                     Jacobian(iJacob) = Jacobian(iJacob) - ElemStiffAll(k) *   &
                            Vab(iData, LocalNodes(j)) * Vmn(iData, LocalNodes(i))
                 end do   ! i
              end do      ! j
           end do         ! iData

        end do         ! iElemY
     end do            ! iElemX
  end do               ! ipy
  end do               ! ipx

  if (allocated(Vab)) deallocate(Vab)
  if (allocated(Vmn)) deallocate(Vmn)

  return
  end

!-------------------------------------------------------------------------------- 
  subroutine JacobFE2(Jacobian, StingCMD, ParamX1, ParamX2, ParamY1, ParamY2,  &
                      weight, ElemStiffAll)

  use GlobalForw
  implicit none

  real(Rkind), intent(inout), dimension(1:gNumData*gNumParam) :: Jacobian
  integer, intent(in), dimension(1:gNumData*4)     :: StingCMD
  integer, intent(in), dimension(1:gNumParamX) :: ParamX1,ParamX2
  integer, intent(in), dimension(1:gNumParamY) :: ParamY1,ParamY2
  real(Rkind), intent(inout), dimension(1:gNumElem*16) :: ElemStiffAll
  real(Rkind), intent(in) ::  weight

  integer :: iParam, iData, iElemX, iElemY, iElem, iJacob0, iJacob
  integer :: i, j, k, iPos, ipx, ipy
  integer ::  LocalNodes(4)
  real(Rkind), allocatable, dimension(:) :: Vab, Vmn  
  integer, allocatable, dimension(:) :: A, B, M, N  
!----------------------------------------------------------------

  if (allocated(Vab)) deallocate(Vab)
  allocate(Vab(gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Vmn)) deallocate(Vmn)
  allocate(Vmn(gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(A)) deallocate(A)
  allocate(A(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(B)) deallocate(B)
  allocate(B(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(M)) deallocate(M)
  allocate(M(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(N)) deallocate(N)
  allocate(N(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired


! Assuming a current electrode at A, B, M, N of iData. Below are the 
! beginning array index of corresponding solution in gfVoverI

  do iData = 1, gNumData               
     k = (iData-1) * 4                 
     A(iData) = StingCMD(k+1) * gNumNodes   
     B(iData) = StingCMD(k+2) * gNumNodes 
     M(iData) = StingCMD(k+3) * gNumNodes   
     N(iData) = StingCMD(k+4) * gNumNodes
  end do
    
  ElemStiffAll = ElemStiffAll * weight


 
  do iData = 1, gNumData

     do j = 1, gNumNodes
        Vab(j) = gfVoverI(A(iData)+j) - gfVoverI(B(iData)+j)
        Vmn(j) = gfVoverI(M(iData)+j) - gfVoverI(N(iData)+j)
     end do

     do ipx = 1, gNumParamX
     do ipy = 1, gNumParamY

        iParam = (ipy - 1) * gNumParamX + ipx
        iJacob = iParam - gNumParam + gNumParam * iData

!       Unscramble the parameter lumping scheme for element-based stiffness
!       ParamX1,ParamX2,ParamY1,ParamY2 are Elemental bounds of current parameter 
        do iElemY = ParamY1(ipy), ParamY2(ipy)
           do iElemX = ParamX1(ipx), ParamX2(ipx)

              iElem = iElemY + gNumElemY * (iElemX - 1)
              LocalNodes(1) = (iElemX-1) * gNumNodeY + iElemY
              LocalNodes(2) = LocalNodes(1) + 1
              LocalNodes(3) = LocalNodes(1) + gNumNodeY
              LocalNodes(4) = LocalNodes(2) + gNumNodeY

              iPos = 16 * (iElem - 1) - 4   

              ! Get local stiffness for current quadrilateral element (iElem)
              ! call StiffnessFE(NodeX, NodeY, Cond, ElemStiff, abscix, iElemX, &
              !      iElemY, LocalNodes, CenterNodeX, CenterNodeY, ElemArea)
              ! Loop through Reduced Stiffness matrix of current elements
              do j = 1, 4
                 k = iPos + 4 * j
                 do i = 1, 4
                     k = k + 1
                     Jacobian(iJacob) = Jacobian(iJacob) - ElemStiffAll(k) *   &
                            Vab(LocalNodes(j)) * Vmn(LocalNodes(i))
                 end do   ! i
              end do      ! j
           end do      ! iElemY
        end do         ! iElemX
     end do            ! ipy
     end do            ! ipx
  end do               ! iData

  if (allocated(Vab)) deallocate(Vab)
  if (allocated(Vmn)) deallocate(Vmn)
  if (allocated(A)) deallocate(A)
  if (allocated(B)) deallocate(B)
  if (allocated(M)) deallocate(M)
  if (allocated(N)) deallocate(N)

  return
  end

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  subroutine JacobFD(Jacobian, dx, dy, Cond, StingCMD, ParamX1, ParamX2,  &
                     ParamY1, ParamY2, abscissa, weight)

  use GlobalForw
  implicit none

  real(Rkind), intent(inout), dimension(1:gNumData*gNumParam) :: Jacobian
  real(Rkind), intent(in),  dimension(1:gNumNodes) :: dx, dy
  real(Rkind), intent(in),  dimension(1:gNumElem)  :: Cond
  integer, intent(in), dimension(1:gNumData*4)     :: StingCMD
  integer, intent(in), dimension(1:gNumParamX) :: ParamX1,ParamX2
  integer, intent(in), dimension(1:gNumParamY) :: ParamY1,ParamY2
  real(Rkind), intent(in) :: abscissa, weight

  integer :: iData, iElemX, iElemY, iElem, iJacob
  integer :: n2, n3, n4, i, k, ipx, ipy
  real(Rkind) :: c, HalfAbscAbsc
  real(Rkind), allocatable, dimension(:) :: Vab, Vmn
  real(Rkind), allocatable, dimension(:) :: Cx, Cy, Cself
  integer, allocatable, dimension(:) :: A, B, M, N, n1

!-----------------------------------------------------
  if (allocated(Vab)) deallocate(Vab)
  allocate(Vab(gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Vmn)) deallocate(Vmn)
  allocate(Vmn(gNumNodes), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Cx)) deallocate(Cx)
  allocate(Cx(gNumElem), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Cy)) deallocate(Cy)
  allocate(Cy(gNumElem), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(Cself)) deallocate(Cself)
  allocate(Cself(gNumElem), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(A)) deallocate(A)
  allocate(A(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(B)) deallocate(B)
  allocate(B(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(M)) deallocate(M)
  allocate(M(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(N)) deallocate(N)
  allocate(N(gNumData), stat=k)
  if (k /= 0) return                 ! ???: err message desired

  if (allocated(n1)) deallocate(n1)
  allocate(n1(gNumElem), stat=k)
  if (k /= 0) return                 ! ???: err message desired
!-----------------------------------------------------

  HalfAbscAbsc = 0.5D0 * abscissa * abscissa

!  Precalculate stiffness factors
!  Unscramble the parameter lumping scheme for element-based stiffness
!  ParamX1,ParamX2,ParamY1,ParamY2 are Elemental bounds of current parameter 
  do iElemX = 1, gNumElemX
     do iElemY = 1, gNumElemY
        iElem = iElemY + (iElemX-1) * gNumElemY


!                     n1            n3 
!                       ------------
!                      |            | 
!                      |   iElem    |
!                      |            |
!                       ------------
!                     n2            n4


        n1(iElem) = iElemY + (iElemX-1) * gNumNodeY
        c  = Cond(iElem)        
        Cx(iElem)    =  -c * dy(iElem) / dx(iElem)
        Cy(iElem)    =  -c * dx(iElem) / dy(iElem)             
        Cself(iElem) =  HalfAbscAbsc * c * dx(iElem) * dy(iElem) - Cx(iElem) - Cy(iElem)
	 end do
  end do

  do iData = 1, gNumData               
     ! Assuming a current electrode at A, B, M, N of iData. Below are the 
     ! beginning array index of corresponding solution in gfVoverI
     k = (iData-1) * 4                
     A(iData) = StingCMD(k+1) * gNumNodes   
     B(iData) = StingCMD(k+2) * gNumNodes 
     M(iData) = StingCMD(k+3) * gNumNodes   
     N(iData) = StingCMD(k+4) * gNumNodes  
  end do

!-----------------------------------------------------
! Loop through data

  do iData = 1, gNumData               

     do i = 1, gNumNodes
        Vab(i) = gfVoverI(A(iData)+i) - gfVoverI(B(iData)+i)
        Vmn(i) = gfVoverI(M(iData)+i) - gfVoverI(N(iData)+i)
     end do

     iJacob = (iData - 1) * gNumParam

     do ipy = 1, gNumParamY
     do ipx = 1, gNumParamX
        iJacob = iJacob + 1
!       Unscramble the parameter lumping scheme for element-based stiffness
!       ParamX1,ParamX2,ParamY1,ParamY2 are Elemental bounds of current parameter 
        do iElemY = ParamY1(ipy), ParamY2(ipy)
           do iElemX = ParamX1(ipx), ParamX2(ipx)

              iElem = iElemY + (iElemX-1) * gNumElemY
              n2 = n1(iElem) + 1
              n3 = n1(iElem) + gNumNodeY
              n4 = n3 + 1

              Jacobian(iJacob) = Jacobian(iJacob) - weight * (                   &
                Vab(n1(iElem)) * (Vmn(n1(iElem)) * Cself(iElem) +                & 
					      Vmn(n3) * Cx(iElem) + Vmn(n2) * Cy(iElem)) +           &
                Vab(n2) * (Vmn(n2) * Cself(iElem) +                              &
					      Vmn(n4) * Cx(iElem) + Vmn(n1(iElem)) * Cy(iElem)) +    &
                Vab(n3) * (Vmn(n3) * Cself(iElem) + Vmn(n1(iElem)) * Cx(iElem) + &
					      Vmn(n4) * Cy(iElem)) +                                 &
                Vab(n4) * (Vmn(n4) * Cself(iElem) + Vmn(n2) * Cx(iElem) +        &
					             Vmn(n3) * Cy(iElem)) )

            end do         ! iElemY
         end do            ! iElemX

     end do                ! ipx
     end do                ! ipy

  end do                   ! iData

  if (allocated(Vab)) deallocate(Vab)
  if (allocated(Vmn)) deallocate(Vmn)
  if (allocated(Cx)) deallocate(Cx)
  if (allocated(Cy)) deallocate(Cy)
  if (allocated(Cself)) deallocate(Cself)
  if (allocated(A)) deallocate(A)
  if (allocated(B)) deallocate(B)
  if (allocated(M)) deallocate(M)
  if (allocated(N)) deallocate(N)
  if (allocated(n1)) deallocate(n1)

  return
  end
