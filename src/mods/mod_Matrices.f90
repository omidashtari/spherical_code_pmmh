! spherical_code_pmmh: solve the Boussinesq equations in a rotating spherical shell
! Copyright (C) 2025 Alan Riquier, Camille Rambert, Juan Cruz Gonzalez Sembla and Laurette Tuckerman

! Contact: juan-cruz.gonzalez-sembla@espci.fr, laurette.tuckerman@espci.fr

! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

! This module was written by A. Riquier, C. Rambert and J. C. Gonzalez Sembla.
! It contains the many algebraic functions used in the code. Functions:
! - MatrixProduct: product of 2 matrices.
! - MatVecProduct: product between matrix and vector.
! - VecMatProduct: product between vector and matrix.
! - VecVecProduct: product between vectors.
! - invL: compute inverse of matrix.
! - Banded_LU_solve: backsolve for banded LU decomposition.
! - Banded_LU_decomp: banded LU decomposition.
! - Banded_Mult: banded matrix multiplication.
! Subroutines:
! - small2big: turn E and F arrays into big array that respects symmetry clases for implicit Coriolis.
! - big2small: turn big array respecting symmetry clases for implicit Coriolis into E and F.

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
  integer :: m, i, m_idx
  integer :: t_odd, t_even, l_odd, l_even, n_even, n_odd
     
  t_odd = 0
  t_even = 0
  m_idx = 0
  step = 2 * (KK2 + KK4)
  A = 0.0d0

  do m = 0, MM*mres, mres
    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    i = 1 - mod(max(m, 1), 2)
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    do l_odd = 1, n_odd
        A(2 * KK4 * i + (l_odd - 1) * step + 1 : &
          & 2 * KK4 * i + (l_odd-1) * step + KK2, m_idx) = &
          & real(E(:, idx_odd(t_odd + l_odd)))
        A(2 * KK4 * i + (l_odd - 1) * step + KK2 + 1 : & 
          & 2 * KK4 * i + (l_odd - 1) * step + 2 * KK2, m_idx) = &
          & aimag(E(:, idx_odd(t_odd + l_odd)))
        A(2 * KK2 * i + mid + (l_odd - 1) * step + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + KK4, m_idx) = &
          & real(F(:, idx_odd(t_odd + l_odd)))
        A(2 * KK2 * i + mid + (l_odd - 1) * step + KK4 + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + 2 * KK4, m_idx) = &
          & aimag(F(:, idx_odd(t_odd + l_odd)))
    end do
    do l_even = 1, n_even
        A(2 * KK2 * (1 - i) + (l_even - 1) * step + 1 : & 
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + KK4, m_idx) = &
          & real(F(:, idx_even(t_even + l_even)))
        A(2 * KK2 * (1 - i) + (l_even - 1) * step + KK4 + 1 : &
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + 2 * KK4, m_idx) = &
          & aimag(F(:, idx_even(t_even + l_even)))
        A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2, m_idx) = &
          & real(E(:, idx_even(t_even + l_even)))
        A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2 + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 2 * KK2, m_idx) = &
          & aimag(E(:, idx_even(t_even + l_even)))
    end do
    t_even = t_even + n_even
    t_odd = t_odd + n_odd
    m_idx = m_idx + 1
  end do

end subroutine small2big

subroutine big2small(A, E, F)
  
  implicit none
  
  double precision, dimension(2 * LL * (KK2 + KK4), 0 : MM), intent(in) :: A
  double complex, dimension(KK2, shtns%nlm), intent(out)                     :: E
  double complex, dimension(KK4, shtns%nlm), intent(out)                     :: F
  integer :: step
  integer :: mid
  integer :: m, i, m_idx
  integer :: t_odd, t_even, l_odd, l_even, n_even, n_odd
     
  t_odd = 0
  t_even = 0
  m_idx = 0
  step = 2 * (KK2 + KK4)
  E = 0.0d0
  F = 0.0d0

  do m = 0, MM*mres, mres
    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    i = 1 - mod(max(m, 1), 2)
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    do l_odd = 1, n_odd
        E(:, idx_odd(t_odd + l_odd))%re = &
          & A(2 * KK4 * i + (l_odd - 1) * step + 1 : &
          & 2 * KK4 * i + (l_odd-1)*step + KK2, m_idx)
        E(:, idx_odd(t_odd + l_odd))%im = & 
          & A(2 * KK4 * i + (l_odd - 1) * step + KK2 + 1 : & 
          & 2 * KK4 * i + (l_odd - 1) * step + 2 * KK2, m_idx)
        F(:, idx_odd(t_odd + l_odd))%re = & 
          & A(2 * KK2 * i + mid + (l_odd - 1) * step + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + KK4, m_idx)
        F(:, idx_odd(t_odd + l_odd))%im = & 
          & A(2 * KK2 * i + mid + (l_odd - 1) * step + KK4 + 1 : & 
          & 2 * KK2 * i + mid + (l_odd - 1) * step + 2 * KK4, m_idx)
    end do
    do l_even = 1, n_even
        F(:, idx_even(t_even + l_even))%re = & 
          & A(2 * KK2 * (1 - i) + (l_even - 1) * step + 1 : & 
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + KK4, m_idx)
        F(:, idx_even(t_even + l_even))%im = & 
          & A(2 * KK2 * (1 - i) + (l_even - 1) * step + KK4 + 1 : &
          & 2 * KK2 * (1 - i) + (l_even - 1) * step + 2 * KK4, m_idx)
        E(:, idx_even(t_even + l_even))%re = & 
          & A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2, m_idx)
        E(:, idx_even(t_even + l_even))%im = & 
          & A(mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + KK2 + 1 : &
          & mid + 2 * KK4 * (1 - i) + (l_even - 1) * step + 2 * KK2, m_idx)
    end do
    t_even = t_even + n_even
    t_odd = t_odd + n_odd
    m_idx = m_idx + 1
  end do

end subroutine big2small  

end module mod_Matrices
