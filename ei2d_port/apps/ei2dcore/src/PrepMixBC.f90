
! Preparation for implementation of mixed boundary condition.
! The coefficient eta is calculated in this file
! Ref(1): Luo and Zhang (1987), pages 68, 82 and 83
! This subroutine should be called from GUI since it does NOT
! change during the inversion. This factor sought here is named
! as eta according to Luo and Zhang (1987).
! Ref(2): Dey and Morrison (1979), Geophysical Prospecting, 27, 106-136
!
! This subroutine calculates (where k = abscissa) 
!
!  k * K1(kr) * cos(theta) / K0(kr)

  subroutine PrepMixBC(NodeX, NodeY, ElecX, ElecY)

!DEC$ ATTRIBUTES C, DLLEXPORT :: PrepMixBC
!DEC$ ATTRIBUTES ALIAS: "prepmixbc" :: PrepMixBC

  use GlobalForw
  implicit none

  real(Rkind), intent(in), dimension(1:gNumNodes) :: NodeX, NodeY
  real(Rkind), intent(in), dimension(1:gNumElectrodes) :: ElecX, ElecY

  real(Rkind) :: BesselK0, BesselK1
  real(Rkind) :: Abscissa, AvgElecX, AvgElecY  ! Average electrode location
  real(Rkind) :: AbscDist, CosTheta
  integer :: i, j, k, iAbscissa


  AvgElecX = 10.0D0
  AvgElecY = 10.0D0
  AvgElecX = sum(ElecX) / gNumElectrodes
  AvgElecY = sum(ElecY) / gNumElectrodes

  if (allocated(gfEtaLeft)) deallocate(gfEtaLeft)
  allocate(gfEtaLeft(gNumNodeY, gNumAbscissa), stat=k)
  if (k /= 0) return                 ! ???: err message desired
  gfEtaLeft = 0D0

  if (allocated(gfEtaRight)) deallocate(gfEtaRight)
  allocate(gfEtaRight(gNumNodeY, gNumAbscissa), stat=k)
  if (k /= 0) return                 ! ???: err message desired
  gfEtaRight = 0D0

  if (allocated(gfEtaBottom)) deallocate(gfEtaBottom)
  allocate(gfEtaBottom(gNumNodeX, gNumAbscissa), stat=k)
  if (k /= 0) return                 ! ???: err message desired
  gfEtaBottom = 0D0

  do iAbscissa = gNumAbscissa-1, gNumAbscissa

     Abscissa = gfAbscissa(iAbscissa)
!  Left side boundary nodes
     do i = 1, gNumNodeY
        AbscDist = sqrt((NodeX(i)-AvgElecX)**2 + (NodeY(i)-AvgElecY)**2)
        CosTheta = abs(NodeX(i) - AvgElecX) / AbscDist
        AbscDist = Abscissa * AbscDist
        gfEtaLeft(i, iAbscissa) = Abscissa * BesselK1(AbscDist) * CosTheta / BesselK0(AbscDist)
     end do

!  Right side boundary nodes
     k = (gNumNodeX-1)*gNumNodeY
     do i = 1, gNumNodeY
        j = k + i
        AbscDist = sqrt((NodeX(j)-AvgElecX)**2 + (NodeY(j)-AvgElecY)**2)
        CosTheta = abs(NodeX(j) - AvgElecX) / AbscDist
        AbscDist = Abscissa * AbscDist
        gfEtaRight(i, iAbscissa) = Abscissa * BesselK1(AbscDist) * CosTheta / BesselK0(AbscDist)
     end do

!  Bottom boundary nodes
     do i = 1, gNumNodeX
        j = i * gNumNodeY
        AbscDist = sqrt((NodeX(j)-AvgElecX)**2 + (NodeY(j)-AvgElecY)**2)
        CosTheta = abs(NodeY(j) - AvgElecY) / AbscDist
        AbscDist = Abscissa * AbscDist
        gfEtaBottom(i, iAbscissa) = Abscissa * BesselK1(AbscDist) * CosTheta / BesselK0(AbscDist)
     end do

  end do      ! iNumAbscissa

  end subroutine PrepMixBC


!============================================================================
!  Calculate the 0-th order Bessel function
  function BesselK0(x)

  use GlobalForw, only : Rkind
  implicit none

  real(Rkind) :: BesselK0
  real(Rkind), intent(in) :: x
  real(Rkind) :: x1, x2, x3, xx

  x1 = 2.0D0 / x
  x2 = x*x / (3.75D0*3.75D0)
  x3 = x * x / 4.0D0
  
  if (x < 2.0) then
     BesselK0 = -log(x/20D0) * (1.0D0 + x2 * (3.5156229D0 + x2 * (3.0899424D0 + x2 *  &
                (1.2067492D0 + x2 * (0.2659732D0 + x2 * (0.0360768D0 + x2 *           & 
                0.0045813D0 ) ) ) ) ) ) - 0.57721566D0 + x3 * (0.42278420D0 + x3 *    &
                (0.23069756D0 + x3 * (0.03488590D0 + x3 * (0.00262698D0 + x3 *        &
                (0.00010750D0 + x3 * 0.00000740D0 ) ) ) ) )  
  else
     xx = x
     if (xx > 700.0D0) xx = 700.0D0
     BesselK0 = (1.25331414D0 + x1 * (-0.07832358D0 + x1 * (0.02189568D0 + x1  *    &
                (-0.01062446D0 + x1 * (0.00587872D0 + x1 * (-0.00251540D0 + x1 *    &
                0.00053208D0 ) ) ) ) ) )  / (sqrt(xx) * exp(xx))
  end if 

  end function BesselK0


!============================================================================
!  Calculate the first order Bessel function
  function BesselK1(x)
  
  use GlobalForw, only : Rkind
  implicit none

  real(Rkind) :: BesselK1
  real(Rkind), intent(in) :: x
  real(Rkind) :: x1, x2, x3, xx

  x1 = 2.0D0 / x
  x2 = x*x / (3.75D0*3.75D0)
  x3 = x*x / 4.0D0
  if (x < 2.0D0) then
     BesselK1 = x * (0.5D0 + x2 * (0.87890594D0 + x2 * (0.51498869D0 + x2      *    &
                (0.15084934D0 + x2 * (0.02658733D0 + x2 * (0.00301532D0 + x2   *    &
                0.00032411D0 ) ) ) ) ) ) * log(x/2D0) + 1.0/x * (1D0 + x3 *         &
                (0.15443144D0 + x3 * (-0.67278579D0 + x3 * (-0.18156897D0 + x3 *    &
                (-0.01919402D0 + x3 * (-0.00110404D0 - x3 * 0.00004686D0 ) ) ) ) ) )
  else
     xx = x
     if (xx > 700.0D0) xx = 700.0D0
     BesselK1 = (1.25331414D0 + x1 * (0.23498619D0 + x1 * (-0.03655620D0 + x1 *     &
                (0.01504268D0 + x1 * (-0.00780353D0 + x1 * (0.00325614D0 - x1 *     &
                 0.00068245D0 ) ) ) ) ) ) / (sqrt(xx) * exp(xx))              
  end if  

  end function BesselK1


