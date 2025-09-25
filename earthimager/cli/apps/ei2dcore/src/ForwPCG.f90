!  Solve Ax=b using Conjugate gradient method with an SSOR 
!  preconditioner. There are three subroutines in this file
!
!  Ref: Golub and van Loan, 1996. Matrix Computation
!       Barrett et al., 1994. Temperates for the solution of ...
!-----------------------------------------------------------------
  subroutine ForwPCG(A, x, b, StiffLen)
  USE GlobalForw
  implicit none

  integer, intent(in) :: StiffLen
  real(Rkind), intent(in),    dimension(1:StiffLen)  :: A
  real(Rkind), intent(in),    dimension(1:gNumNodes) :: b
  real(Rkind), intent(inout), dimension(1:gNumNodes) :: x
!
  real(Rkind) :: dot1,dot2, pAp, alpha, beta
  real(Rkind), allocatable, dimension(:) :: r,p,precond,ApVec
  real(Rkind) :: omega, residual           
  integer :: iter

!----------------------------------------------------------
  allocate(r(1:gNumNodes), stat=iter)
  r = 0D0

  allocate(p(1:gNumNodes), stat=iter)
  p = 0D0

  allocate(precond(1:gNumNodes), stat=iter)
  precond = 0D0

  allocate(ApVec(1:gNumNodes), stat=iter)
  ApVec = 0D0

! SSOR relaxation factor
  omega = 1.5       
  if (gNumNodes < 10000) omega = 1.40D0
  if (gNumNodes > 50000) omega = 1.60D0

! stiffness matrix times a vector (initial guess /= 0) r = Ax

  if (gnForwModMeth == 0) then
    call StiffMatrixVecFD(A,x,r,StiffLen)
  else
    call StiffMatrixVecFE(A,x,r,StiffLen)
  end if

  r = b - r
  dot1 = 0D0

  do iter = 1, gnForwCGIter


     if (gnForwModMeth == 0) then
        call SSORFD(A, r, precond, omega, StiffLen)
     else
        call SSORFE(A, r, precond, omega, StiffLen)
     end if

     dot2 = dot1
     dot1 = dot_product(precond, r)
     if (iter == 1) then
        p = precond
     else
        beta = dot1 / dot2
        p = precond + beta * p
     end if

     ApVec = 0D0

     if (gnForwModMeth == 0) then
        call StiffMatrixVecFD(A,p,ApVec, StiffLen)
     else
        call StiffMatrixVecFE(A,p,ApVec, StiffLen)
     end if

     pAp = dot_product(p, ApVec)
     alpha = dot1 / pAp
     x = x + alpha * p
     r = r - alpha * ApVec
     residual = sqrt(dot_product(r, r)/gNumNodes)
     if(residual < gfForwCGResid) exit
  end do

  if (allocated(r)) deallocate(r)
  if (allocated(p)) deallocate(p)
  if (allocated(precond)) deallocate(precond)
  if (allocated(ApVec)) deallocate(ApVec)

  return
  end subroutine ForwPCG


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  SSOR -- Symmetric successive over-relaxation preconditioning.
!  Solve for a preconditioner from a L * transp(L) x = b
!  L is a scaled lower triangle of stiffness matrix (see also forward.f)
!
!  Ref. Weller, A., Seichter, M., Kampke, A., 1996. IP modeling using 
!       complex electrical conducitvities. Geoph. J. Int., 127, 387-398.
!
!       Barrett et al., 1994. Temperates for the solution of linear 
!       systems: building blocks for iterative methods, SIAM.
! ----------------------------------------------------------------------
  subroutine SSORFE(A, b, x, omega, StiffLen)   !  Finite element method
  use GlobalForw
  implicit none

  integer, intent(in) :: StiffLen
  real(Rkind), intent(in),  dimension(1:StiffLen)  :: A
  real(Rkind) ,intent(in),  dimension(1:gNumNodes) :: b
  real(Rkind) ,intent(out), dimension(1:gNumNodes) :: x
  real(Rkind), intent(in) :: omega

  integer :: i,i2, i3, i4, i5, n2, n3, n4, n5

  x = 0D0
  x(1) = b(1)

  do i = 2, gNumNodes                !  forward substitution

     i2 = i-1                        !  indexing x
     i3 = i-gNumNodeY+1
     i4 = i3 - 1
     i5 = i3 - 2

     ! indexing 1D stiffness matrix
     n2 = i2 * gNumDiagonal - 3 
     n3 = i3 * gNumDiagonal - 2
     n4 = i4 * gNumDiagonal - 1
     n5 = i5 * gNumDiagonal

     if ( i5 > 0) then
        x(i) = b(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3) + A(n4)*x(i4) + A(n5)*x(i5))
     else if (i4 > 0) then
        x(i) = b(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3) + A(n4)*x(i4))
     else if (i3 > 0) then
        x(i) = b(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3))
     else if (i2 > 0) then
        x(i) = b(i) - omega * A(n2) * x(i2)
     end if

  end do

  do i = gNumNodes-1, 1, -1        !  backward substitution

     i2 = i + 1                    !  indexing x and b
     i3 = i+gNumNodeY - 1
     i4 = i3 + 1
     i5 = i3 + 2

     ! indexing 1D stiffness matrix
     n2 = i * gNumDiagonal - 3  
     n3 = n2 + 1
     n4 = n2 + 2
     n5 = n2 + 3

     if (i5 <= gNumNodes) then
        x(i) = x(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3) + A(n4)*x(i4) + A(n5)*x(i5))
     else if (i4 <= gNumNodes) then
        x(i) = x(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3) + A(n4)*x(i4))
     else if (i3 <= gNumNodes) then
        x(i) = x(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3))
     else if (i2 <= gNumNodes) then
        x(i) = x(i) - omega * A(n2) * x(i2)
     end if

  end do

  return
  end subroutine SSORFE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  SSOR -- Symmetric successive over-relaxation preconditioning for
!  Finite difference method
  subroutine SSORFD(A, b, x, omega, StiffLen)   
  use GlobalForw
  implicit none

  integer, intent(in) :: StiffLen
  real(Rkind), intent(in),  dimension(1:StiffLen)  :: A
  real(Rkind) ,intent(in),  dimension(1:gNumNodes) :: b
  real(Rkind) ,intent(out), dimension(1:gNumNodes) :: x
  real(Rkind), intent(in) :: omega

  integer :: i,i2, i3, n2, n3

  x = 0D0
  x(1) = b(1)

  do i = 2, gNumNodes                !  forward substitution

     i2 = i-1                        !  indexing x
     i3 = i-gNumNodeY

     ! indexing 1D stiffness matrix
     n2 = i2 * gNumDiagonal - 1 
     n3 = i3 * gNumDiagonal

     if ( i3 > 0) then
        x(i) = b(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3))
     else if (i2 > 0) then
        x(i) = b(i) - omega * A(n2) * x(i2)
     end if

  end do

  do i = gNumNodes-1, 1, -1      !  backward substitution

     i2 = i+1                    !  indexing x and b
     i3 = i+gNumNodeY

     ! indexing 1D stiffness matrix
     n2 = i * gNumDiagonal - 1  
     n3 = i * gNumDiagonal

     if (i3 <= gNumNodes) then
        x(i) = x(i) - omega * (A(n2)*x(i2) + A(n3)*x(i3))
     else if (i2 <= gNumNodes) then
        x(i) = x(i) - omega * A(n2) * x(i2)
     end if

  end do

  return
  end subroutine SSORFD


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Compute stiffness matrix times a vector x. 
!
!       b = StiffMatrix * x
!
!  Note that StiffMatrix is stored in a 1D compact form of length
!  gNumDiagonal * gNumNodes.
!---------------------------------------------------------------------------
! For Finite element method
  subroutine StiffMatrixVecFE(A, x, b, StiffLen)   

  USE GlobalForw
  implicit none

  integer, intent(in) :: StiffLen
  real(Rkind), intent(in),  dimension(1:StiffLen)  :: A
  real(Rkind), intent(in),  dimension(1:gNumNodes) :: x
  real(Rkind), intent(out), dimension(1:gNumNodes) :: b

  integer :: i_2, i_3, i_4, i, i1, i2, i3, i4
  integer :: n_1, n_2, n_3, n_4, n

  b = 0D0

! First row, band1
  b(1) = A(1)*x(1) + A(2)*x(2) + A(3)*x(gNumNodeY) +    & 
         A(4)*x(gNumNodeY+1) + A(5)*x(gNumNodeY+2)
  
! Row 2 to Row gNumNodeY - 1
  do i = 2, gNumNodeY - 1
     i2 = i + gNumNodeY - 1
     i3  = i2 + 1 
     i4  = i2 + 2
     n = (i-1)*gNumDiagonal + 1
     n_1 = n - gNumDiagonal + 1

     b(i) = A(n_1)*x(i-1) + A(n)*x(i) + A(n+1)*x(i+1) + A(n+2)*x(i2) +   & 
            A(n+3)*x(i3) + A(n+4)*x(i4)
  end do 

! Row gNumNodeY
  i   = gNumNodeY
  i2  = i + gNumNodeY - 1
  i3  = i2 + 1 
  i4  = i2 + 2
  n   = (i-1)*gNumDiagonal + 1
  n_1 = n - gNumDiagonal + 1

  b(i) = A(3)*x(1) + A(n_1)*x(i-1) + A(n)*x(i) + A(n+1)*x(i+1) +         & 
                 A(n+2)*x(i2) + A(n+3)*x(i3) + A(n+4)*x(i4) 

! Row gNumNodeY + 1
  i   = gNumNodeY + 1
  i2  = i + gNumNodeY - 1
  i3  = i2 + 1 
  i4  = i2 + 2
  n   = (i-1)*gNumDiagonal + 1
  n_1 = n - gNumDiagonal + 1
  n_2 = gNumDiagonal + 3

  b(i) = A(4)*x(1) + A(n_2)*x(2) + A(n_1)*x(i-1) + A(n)*x(i) +           & 
         A(n+1)*x(i+1) + A(n+2)*x(i2) + A(n+3)*x(i3) + A(n+4)*x(i4) 

! Rows with a full band
  do i = gNumNodeY+2, gNumNodes-gNumNodeY-1
     i2  = i + gNumNodeY - 1
     i3  = i2 + 1 
     i4  = i2 + 2
     i_2 = i - gNumNodeY + 1
     i_3 = i_2 - 1
     i_4 = i_2 - 2
     n   = (i-1) * gNumDiagonal + 1
     n_1 = (i-2) * gNumDiagonal + 2
     n_2 = i_2 * gNumDiagonal - 2
     n_3 = i_3 * gNumDiagonal - 1
     n_4 = i_4 * gNumDiagonal

     b(i) = A(n_4)*x(i_4) + A(n_3)*x(i_3) + A(n_2)*x(i_2) + A(n_1)*x(i-1) +    & 
            A(n)*x(i) +A(n+1)*x(i+1) + A(n+2)*x(i2) + A(n+3)*x(i3) + A(n+4)*x(i4) 
  end do

! The rest of rows with a added triangle of zeros in stiffness
! ???: The triangular area outside the stiffness matrix is filled w/ 0s
  do i = gNumNodes-gNumNodeY, gNumNodes
     i1  = min(gNumNodes, i+1)
     i2  = min(gNumNodes, i+gNumNodeY-1)
     i3  = min(gNumNodes, i2+1) 
     i_2 = i - gNumNodeY + 1
     i_3 = i_2 - 1
     i_4 = i_2 - 2
     n   = (i-1) * gNumDiagonal + 1
     n_1 = (i-2) * gNumDiagonal + 2
     n_2 = i_2 * gNumDiagonal - 2
     n_3 = i_3 * gNumDiagonal - 1
     n_4 = i_4 * gNumDiagonal

     b(i) = A(n_4)*x(i_4) + A(n_3)*x(i_3) + A(n_2)*x(i_2) + A(n_1)*x(i-1) +    & 
            A(n)*x(i) +A(n+1)*x(i1) + A(n+2)*x(i2) + A(n+3)*x(i3) 
  end do

  return
  end subroutine StiffMatrixVecFE


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Stiffness matrix times a vector for Finite difference method

  subroutine StiffMatrixVecFD(A, x, b, StiffLen)   

  USE GlobalForw
  implicit none

  integer, intent(in) :: StiffLen
  real(Rkind), intent(in),  dimension(1:StiffLen)  :: A
  real(Rkind), intent(in),  dimension(1:gNumNodes) :: x
  real(Rkind), intent(out), dimension(1:gNumNodes) :: b

  integer :: i_2, i, i1, i2
  integer :: n_1, n_2, n

  b = 0D0

! First row
  b(1) = A(1)*x(1) + A(2)*x(2) + A(3)*x(gNumNodeY+1)
  
! Row 2 to Row gNumNodeY
  do i = 2, gNumNodeY
     i2 = i + gNumNodeY
     n = (i-1)*gNumDiagonal + 1
     n_1 = n - gNumDiagonal + 1
     b(i) = A(n_1)*x(i-1) + A(n)*x(i) + A(n+1)*x(i+1) + A(n+2)*x(i2) 
  end do 

! The rows with a full band
  do i = gNumNodeY+1, gNumNodes-gNumNodeY
     i2  = i + gNumNodeY
     i_2 = i - gNumNodeY
     n   = (i-1)*gNumDiagonal + 1
     n_1 = n - gNumDiagonal + 1
     n_2 = i_2 * gNumDiagonal
     b(i) = A(n_2)*x(i_2) + A(n_1)*x(i-1) + A(n)*x(i) + A(n+1)*x(i+1) + A(n+2)*x(i2) 
  end do

! The rest of rows
  do i = gNumNodes-gNumNodeY+1, gNumNodes
     i1  = min(gNumNodes, i+1)
     i_2 = i - gNumNodeY
     n   = (i-1)*gNumDiagonal + 1
     n_1 = n - gNumDiagonal + 1
     n_2 = i_2 * gNumDiagonal
     b(i) = A(n_2)*x(i_2) + A(n_1)*x(i-1) + A(n)*x(i) + A(n+1)*x(i1) 
  end do
!---------------------------------------------
  return
  end subroutine StiffMatrixVecFD
