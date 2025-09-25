!
!  There are two subroutines in this file for both finite element and
!  finite difference methods. 
!
!  Subroutine Decompose implements Cholesky decomposition of the 
!  stiffness matrix that is a symmetric, banded, positive, 
!  definite matrix.
!
!  This algorithm can be found in Goloub and van Loan (1996, p. 156),
!  Matrix Computation, and Dahlquist and Bjork (1974, p. 158), 
!  Numerical Methods.  
!
!  For efficient execution, Stiffness matrix is stored in a 1-D array 
!  StiffMatrix(gNumNodes * gnBand). Refer to ForwardFE.f90 and 
!  ForwardFD.f90. 
!  
!  A = L * transposed(L), where L is lower banded triangular matrix
!  
!  Decompose replaces the original StiffMatrix with L. So data in 
!  StiffMatrix disappears.
!-------------------------------------------------------------------------
  subroutine Decompose(A)

  use GlobalForw
  implicit none

  real(Rkind), intent(inout), dimension(1:gNumNodes*gnBand) :: A
  integer :: i, j, k, jBegin, jPos, kPos, band 
  real(Rkind) :: sum, tiny

  tiny = 1.0d-30
  band = gnBand - 1
  kPos = 0
  do k = 1, gNumNodes
     kPos = kPos + band             ! index the 1D StiffMatrix
     jBegin = max(k-gnBand+1, 1)
     jPos = (jBegin - 1) * band
     do j = jBegin, k-1
        jPos = jPos + band
        sum = 0D0
        do i = jBegin, j-1
            sum = sum + A(kPos+i) * A(jPos+i)
        end do
        A(kPos+j) = (A(kPos+j) - sum) / A(jPos+j)
     end do
     sum = 0D0
     do j = jBegin, k-1
        sum = sum + A(kPos+j) * A(kPos+j)
     end do
! ???:     if ((A(kPos+k) - sum) < 0.) call CheckPoint(k, A(kPos+k) - sum)
     A(kPos+k) = sqrt(A(kPos+k) - sum + tiny)
  end do

  return
  end subroutine Decompose


!======================================================================
!  Subroutine Solve solves decomposed matrix system for transformed  
!  electric potential distribution over entire grid by forward and then  
!  backward substitution. A now is the lower banded triangular matrix 
!  decomposed from StiffMatrix.  
!-------------------------------------------------------------------------
  subroutine Solve(A,V)

  use GlobalForw
  implicit none
  real(Rkind), intent(in), dimension(1:gNumNodes*gnBand) :: A
  real(Rkind), intent(inout), dimension(1:gNumNodes) ::  V
  integer :: i, j, iPos, jPos, band, jBegin, jEnd
  real(Rkind) :: sum

!  First solve L * V = V  for V, by forward substitution
!  V to the right is the source vector at the beginning

  band = gnBand - 1
  iPos = 0
  do i = 1, gNumNodes
     iPos = iPos + band
     jBegin = max(i-gnBand+1, 1)
     sum = 0D0
     do j = jBegin, i-1
        sum = sum + A(iPos+j) * V(j)
     end do
     V(i) = (V(i) - sum) / A(iPos+i)
  end do

!  Next solve  Lt * V = V  for V, by backward substitution 
!  Now V to the right is the solution of forward substituion  

  iPos = (gNumNodes + 1) * band
  do i = gNumNodes, 1, -1         ! backward indexing
     iPos = iPos - band
     sum = 0D0
     jEnd = min(i+gnBand-1, gNumNodes)
     jPos = i * band
     do j = i+1, jEnd
        jPos = jPos + band
        sum = sum + A(jPos+i) * V(J)
     end do
     V(i) = (V(i) - sum) / A(iPos+i)
  end do

  return
  end subroutine Solve
