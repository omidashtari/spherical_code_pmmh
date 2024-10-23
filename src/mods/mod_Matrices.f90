module mod_Matrices
 !$ USE OMP_LIB
  use mod_Globalvars
  use mod_PrecompSH
  
  implicit none

  interface operator(.dot.)
    module procedure MatrixProduct, MatVecProduct, VecVecProduct, VecMatProduct
  end interface

contains

  !############################################################################
  !                         Codes for the .dot. operator
  !############################################################################

  function MatrixProduct(A,B) result (AB)

    implicit none

    double precision, dimension(:,:), intent(in)      :: A,B
    double precision, dimension(size(A,1), size(B,2)) :: AB

    integer :: i,j,k

    AB = matmul(A,B)

  end function MatrixProduct

  !-----------------------------------------------------------------------------

  function MatVecProduct(A,X) result(AX)

    implicit none

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(:)  , intent(in) :: X
    double precision, dimension(size(A,1))       :: AX
    integer                                      :: i,j

    if (size(A,2) /= size(X)) error stop                                      &
          'Error .dot. product (Mat.Vec): Wrong shapes...'

    AX = matmul(A,X)

  end function MatVecProduct

  !-----------------------------------------------------------------------------

  function VecMatProduct(X,A) result(XA)

    implicit none

    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(:)  , intent(in) :: X
    double precision, dimension(size(A,2))       :: XA
    integer                                      :: i,j

    if (size(A,1) /= size(X)) error stop                                      &
          'Error .dot. product (Vec.Mat): Wrong shapes...'

    XA = matmul(X,A)

  end function VecMatProduct

  !-----------------------------------------------------------------------------

  function VecVecProduct(A,B) result(AB)

    implicit none

    double precision, dimension(:), intent(in) :: A,B
    double precision                           :: AB
    integer                                    :: i

    if (size(A) /= size(B)) error stop                                      &
          'Error .dot. product (Vec.Vec): Wrong shapes...'

    AB = dot_product(A,B)

  end function VecVecProduct

  !############################################################################
  !
  !                              LU decomposition
  !
  !############################################################################

  subroutine LU_decomp(A,LU,P)

    !---------------- Doolittle algorithm for LU decomposition -----------------

    implicit none

    double precision, dimension(:,:), intent(in) :: A ! Input square matrix
    double precision, dimension(size(A,1),size(A,2)), intent(out) :: LU
    integer, dimension(size(A,1)), intent(out) :: P ! Output permut.

    !--------------------------------------------------------------------------!
    ! L is a lower triangular matrix with ones on the main diagonal            !
    ! U is a upper triangular matrix                                           !
    ! At last, we want PA = LU. In most cases, P = id.                         !
    ! P corresponds to a permutation matrix. It's a vector containing column   !
    !         indexes where the permutation matrix has "1". For example :      !
    !             |1 0 0 0 0|           |1|                                    !
    !             |0 0 0 1 0|           |4|                                    !
    !             |0 0 1 0 0|     =>    |3|                                    !
    !             |0 1 0 0 0|           |2|                                    !
    !             |0 0 0 0 1|           |5|                                    !
    !         Plus, the first element of P is the number of permutations made  !
    !         (1 in our example)                                               !
    ! The values of L and U are stored in A in the form A = (L-id)+U           !
    !--------------------------------------------------------------------------!

    integer :: sizeA
    integer :: i,j,n
    integer :: nMax
    double precision :: max
    double precision :: tol = 1.d-8
    double precision, dimension(size(A,1)) :: piv

    sizeA = size(A,1)

    LU = A

    if (sizeA == size(A,2)) then ! Check if A is a square matrix

      do n = 1,sizeA    ! Initialising P
        P(n) = n
      end do
      do n = 1,sizeA

        nMax = n
        max = 0.d0

        ! To avoid A(i,i) = 0, we put the maximum value on the colum A(:,i) in
        ! the A(i,i) position using permutation. If this can't be done, then the
        ! matrix A is singular

        do j = n,sizeA
          if (abs(LU(j,n)) > max) then
            max = abs(LU(j,n))
            nMax = n
          end if
        end do

        if (max < tol) then
          Print*,"LU decomposition Error: Singular Matrix"
          exit
        end if

        if (nMax /= n) then
          ! Pivoting P
          j = P(n)
          P(n) = P(nMax)
          P(nMax) = j

          ! Pivoting the rows of A
          piv = LU(n,:)
          LU(n,:) = LU(nMax,:)
          LU(nMax,:) = piv

          ! Add a permutation to the counter
          !P(0) = P(0) + 1
        end if
        ! Compute the value of the L and U matrices
        do i = n+1,sizeA
          LU(i,n) = LU(i,n)/LU(n,n)
          do j = n+1,sizeA
            LU(i,j) = LU(i,j) - LU(i,n)*LU(n,j)
          end do
        end do
      end do
    else
      Print*,"LU decomposition Error: Not a square matrix"
    end if

  
    
  end subroutine LU_decomp

  !----------------------------------------------------------------------------

  subroutine giveExplicitLUP(A,L,U,P_vec,P)

    implicit none

    ! Function to get the matrices L, U and P from A and P computed by
    ! the LU_Decomp function

    double precision, dimension(:,:), intent(in) :: A
    integer, dimension(size(A,1)), intent(in) :: P_vec
    double precision, dimension(size(A,1),size(A,1)), intent(out) :: L
    double precision, dimension(size(A,1),size(A,1)), intent(out) :: U
    integer, dimension(size(A,1),size(A,1)), intent(out) :: P

    integer :: N
    integer :: i,j
    integer :: piv

    P = 0

    N = size(A,1)

    do i = 1,N
      piv = P_vec(i)
      P(piv,i) = 1
      do j = i,N
        if (i == j) then
          L(i,i) = 1
        else
          L(j,i) = A(j,i)
        end if
        U(i,j) = A(i,j)
      end do
    end do

  end subroutine giveExplicitLUP

  !*****************************************************************************

  function LU_solve(LU,P,b) result(x)

    implicit none

    ! Function to solve a linear system Ax = b
    ! LU and P are given by the LUP_Decomp subroutine

    double precision, dimension(:,:), intent(in) :: LU
    integer, dimension(size(LU,1)), intent(in) :: P
    double precision, dimension(size(LU,1)), intent(in) :: b
    double precision, dimension(size(LU,1)) :: x
    integer :: i,j
    integer :: N

    ! First solve Ly = Pb. Here, y is stored in x. Since L is lower triangular,
    ! solving the system is easy

    N = size(LU,1)

    do i = 1,N
      x(i) = b(P(i))
      do j = 1,i-1
        x(i) = x(i) - LU(i,j)*x(j)
      end do
    end do

    ! Then solve Ux = y. Again, since U is upper triangular, the procedure is
    ! straightforward

    do i = 0,N-1
      do j = 0,i-1
        x(N-i) = x(N-i) - LU(N-i,N-j)*x(N-j)
      end do
      x(N-i) = x(N-i) / LU(N-i,N-i)
    end do

  end function LU_solve

  !*****************************************************************************

  function LU_inv(LU,P) result(A_inv)

    implicit none

    ! Function returning the Inverse of the matrix A
    ! Lu and P are given from the LUP_decomp subroutine
    ! In fact, solving AX = id is equivalent to LUX = P and we can use the
    !  LUP_solve algorithm for each column

    double precision, dimension(:,:), intent(in) :: LU
    integer, dimension(size(LU,1))  , intent(in) :: P
    double precision, dimension(size(LU,1),size(LU,1)) :: A_inv

    integer :: i,j,k,N

    N = size(LU,1)

    do i = 1,N

      do j = 1,N
        if (P(j) == i) then
          A_inv(j,i) = 1.d0
        else
          A_inv(j,i) = 0.d0
        end if
        do k = 1,j-1
          A_inv(j,i) = A_inv(j,i) - LU(j,k)*A_inv(k,i)
        end do
      end do

      do j = 0,N-1
        do k = 0,j-1
          A_inv(N-j,i) = A_inv(N-j,i) - LU(N-j,N-k)*A_inv(N-k,i)
        end do
        A_inv(N-j,i) = A_inv(N-j,i) / LU(N-j,N-j)
      end do
    end do

  end function LU_inv

  !----------------------------------------------------------------------------

  function inv(A) result (Ai)

    implicit none

    double precision, dimension(:,:), intent(in)     :: A ! Input square matrix
    double precision, dimension(size(A,1),size(A,2)) :: LU
    double precision, dimension(size(A,1),size(A,2)) :: Ai
    integer, dimension(size(A,1))                    :: P ! Output permut.

    call LU_decomp(A,LU,P)
    Ai = LU_inv(LU,P)

  end function inv

  function invL(A) result(Ainv)
    ! Code taken from http://fortranwiki.org/fortran/show/Matrix+inversion
    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(size(A,1),size(A,2)) :: Ainv

    double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if



  end function invL

  !###########################################################################
  !            Thomas algorithm for tridiagonal matrices
  !               TODO: Allow block-tridiagonal matrices
  !###########################################################################

  function Thomas_solve(Diag,Sup,Sub,b) result(x)

    !--- Thomas algorithm for tridiagonal matrices

    implicit none

    double precision, dimension(:), intent(inout) :: Diag
    double precision, dimension(:), intent(in) :: Sup, Sub
    double precision, dimension(:), intent(inout) :: b
    double precision, dimension(size(b))       :: x
    integer :: n
    integer :: i
    double precision :: w ! Dummy variable for computations

    n = size(Diag)

    if (n /= (size(Sup)+1)) error stop                                        &
          'error Thomas algorithm: Wrong size of the Superdiagonal'
    if (n /= (size(Sub)+1)) error stop                                        &
          'error Thomas algorithm: Wrong size of the Subdiagonal'
    if (n /= size(b))       error stop 'Error Thomas algorithm: wrong size for b'

    do i = 2,n
      w = Sub(i-1)/Diag(i-1)
      Diag(i) = Diag(i) - w*Sup(i-1)
      b(i) = b(i) - w*b(i-1)
    end do

    ! Backward substitution
    x(n) = b(n)/Diag(n)
    do i = n-1,1,-1
      x(i) = (b(i) - Sup(i)*x(i+1))/Diag(i)
    end do

  end function Thomas_solve

function Banded_LU_solve(N,D,AB,LDAB, B,IPIV) result(X)
    
  implicit none
  integer, intent(in) :: N, LDAB, D
  double precision, dimension(LDAB, N), intent(in) :: AB
  integer :: INFO
  double precision, dimension(N,1), intent(in) :: B
  double precision, dimension(N,1) :: X
  integer, dimension(N) :: IPIV
  integer :: m

  external DGBTRS
  X=B
   
  call dgbtrs('N', N, D, D, 1, AB, LDAB, IPIV, X, N, INFO)
  ! print*, 'info', INFO

end function Banded_LU_solve


subroutine Banded_LU_decomp(N,D,LDAB,AB,IPIV)
  
  implicit none
  integer, intent(in) :: N, LDAB, D
  double precision, dimension(LDAB, N) :: AB
  integer :: INFO
  integer, dimension(N) :: IPIV
  integer :: m

  external DGBTRF
    
  call dgbtrf(N,N,D,D,AB,LDAB,IPIV,INFO)
  !print*, 'info', INFO

end subroutine Banded_LU_decomp

function Banded_Mult(N, D, A, LDA, X, B) result(Y)
  integer, intent(in) ::  N, LDA, D
  integer  :: INCX,INCY
  double precision, dimension(LDA, N), intent(in) :: A
  double precision, dimension(N,1) :: X, Y, B
  double precision :: alpha, beta
  
  external DGBMV

  alpha=1.d0
  beta=1.d0
  INCX=1
  INCY=1
  
  Y=B
  call dgbmv('N',N,N,D,D,alpha,A,LDA,X,INCX,beta,Y, INCY)
   
end function Banded_Mult

subroutine small2big(E, F, A)
    
  implicit none

  double complex, dimension(KK2, shtns%nlm), intent(in)                       :: E
  double complex, dimension(KK4, shtns%nlm), intent(in)                       :: F
  double precision, dimension(2 * LL * (KK2 + KK4), 0 : MM), intent(out) :: A
  integer :: step
  integer :: mid
  integer :: m, i
  integer :: t_odd, t_even, l_odd, l_even, n_even, n_odd
     
  t_odd = 0
  t_even = 0
  step = 2 * (KK2 + KK4)
  A = 0.0d0
  do m = 0, MM
    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    i = 1 - mod(max(m, 1), 2)
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    do l_odd = 1, n_odd
        A(2 * KK4 * i + (l_odd - 1) * step + 1 : &
          & 2 * KK4 * i + (l_odd-1) * step + KK2, m) = &
          & real(E(:, idx_odd(t_odd + l_odd)))
        A(2 * KK4 * i + (l_odd - 1) * step + KK2 + 1 : & 
          & 2 * KK4 * i + (l_odd - 1) * step + 2 * KK2, m) = &
          & aimag(E(:, idx_odd(t_odd + l_odd)))
        A(2 * KK2 * i + mid + (l_odd - 1) * step + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + KK4, m) = &
          & real(F(:, idx_odd(t_odd + l_odd)))
        A(2 * KK2 * i + mid + (l_odd - 1) * step + KK4 + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + 2 * KK4, m) = &
          & aimag(F(:, idx_odd(t_odd + l_odd)))
    end do
    do l_even = 1, n_even
        A(2 * KK2 * (1 - i) + (l_even - 1) * step + 1 : & 
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + KK4, m) = &
          & real(F(:, idx_even(t_even + l_even)))
        A(2 * KK2 * (1 - i) + (l_even - 1) * step + KK4 + 1 : &
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + 2 * KK4, m) = &
          & aimag(F(:, idx_even(t_even + l_even)))
        A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2, m) = &
          & real(E(:, idx_even(t_even + l_even)))
        A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2 + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 2 * KK2, m) = &
          & aimag(E(:, idx_even(t_even + l_even)))
    end do
    t_even = t_even + n_even
    t_odd = t_odd + n_odd
  end do

end subroutine small2big

subroutine big2small(A, E, F)
  
  implicit none
  
  double precision, dimension(2 * LL * (KK2 + KK4), 0 : MM), intent(in) :: A
  double complex, dimension(KK2, shtns%nlm), intent(out)                     :: E
  double complex, dimension(KK4, shtns%nlm), intent(out)                     :: F
  integer :: step
  integer :: mid
  integer :: m, i
  integer :: t_odd, t_even, l_odd, l_even, n_even, n_odd
     
  t_odd = 0
  t_even = 0
  step = 2 * (KK2 + KK4)
  E = 0.0d0
  F = 0.0d0
  do m = 0, MM
    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    i = 1 - mod(max(m, 1), 2)
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    do l_odd = 1, n_odd
        E(:, idx_odd(t_odd + l_odd))%re = &
          & A(2 * KK4 * i + (l_odd - 1) * step + 1 : &
          & 2 * KK4 * i + (l_odd-1)*step + KK2, m)
        E(:, idx_odd(t_odd + l_odd))%im = & 
          & A(2 * KK4 * i + (l_odd - 1) * step + KK2 + 1 : & 
          & 2 * KK4 * i + (l_odd - 1) * step + 2 * KK2, m)
        F(:, idx_odd(t_odd + l_odd))%re = & 
          & A(2 * KK2 * i + mid + (l_odd - 1) * step + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + KK4, m)
        F(:, idx_odd(t_odd + l_odd))%im = & 
          & A(2 * KK2 * i + mid + (l_odd - 1) * step + KK4 + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + 2 * KK4, m)
    end do
    do l_even = 1, n_even
        F(:, idx_even(t_even + l_even))%re = & 
          & A(2 * KK2 * (1 - i) + (l_even - 1) * step + 1 : & 
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + KK4, m)
        F(:, idx_even(t_even + l_even))%im = & 
          & A(2 * KK2 * (1 - i) + (l_even - 1) * step + KK4 + 1 : &
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + 2 * KK4, m)
        E(:, idx_even(t_even + l_even))%re = & 
          & A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2, m)
        E(:, idx_even(t_even + l_even))%im = & 
          & A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2 + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 2 * KK2, m)
    end do
    t_even = t_even + n_even
    t_odd = t_odd + n_odd
  end do

end subroutine big2small  

end module mod_Matrices
