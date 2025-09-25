
!  This file contains utility functions and subroutines
!  for forward modeling

!----------------------------------------------------------
!  Get cooresponding array index of k in array IndexOffset
  function GetOffsetIndex(k,IndexOffset)
  use GlobalForw, only : gNumDiagonal
  implicit none

  integer :: GetOffsetIndex
  integer, intent(in) :: k
  integer, intent(in), dimension(1:gNumDiagonal) :: IndexOffset
  integer :: i

  i = 1
  do while ((i <= gNumDiagonal) .and. (k /= IndexOffset(i)))
     i = i + 1
  end do
  GetOffsetIndex = i

  end function GetOffsetIndex

!----------------------------------------------------------
! subroutine to write debug info to a file named debug.ech
! ???: This routine should be deleted later

  subroutine CheckPoint(k1, k2)
  implicit none

  integer, intent(in) :: k1, k2
  integer :: iech

  iech = 18
  open(unit=iech, file='C:\cvs\Projects\Volente\debug.ech', status="unknown")
  write(iech, *) 'Stop at ', k1, k2
  close(iech)

  end subroutine CheckPoint

!----------------------------------------------------------
!  Query if an electrode is an infinity electrode

  function IsInfElectrode(iElec, InfElec)
  use GlobalForw, only : gNumInfElectrodes
  implicit none

  logical :: IsInfElectrode
  integer, intent(in) :: iElec
  integer, intent(in), dimension(1:gNumInfElectrodes) :: InfElec
  integer :: i

  IsInfElectrode = .false.
  Loop1 : do i = 1, gNumInfElectrodes
     if (iElec == InfElec(i)) then
        IsInfElectrode = .true.
        exit Loop1
     end if
  end do Loop1

  end function IsInfElectrode

