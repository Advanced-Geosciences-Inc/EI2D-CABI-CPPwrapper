  subroutine FreeMemForward

!DEC$ ATTRIBUTES C, DLLEXPORT :: FreeMemForward
!DEC$ ATTRIBUTES ALIAS: "freememforward" :: FreeMemForward

  use GlobalForw
  implicit none

  if (allocated(gfAbscissa)) deallocate(gfAbscissa)
  if (allocated(gfFTWeight)) deallocate(gfFTWeight)
  if (allocated(gfStiff)) deallocate(gfStiff)
  if (allocated(gfVoverI)) deallocate(gfVoverI)
!  if (allocated(gfVoverIold)) deallocate(gfVoverIold)
  if (allocated(gfEtaLeft)) deallocate(gfEtaLeft)
  if (allocated(gfEtaRight)) deallocate(gfEtaRight)
  if (allocated(gfEtaBottom)) deallocate(gfEtaBottom)

  end subroutine FreeMemForward