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

! This module was written by J. C. Gonzalez Sembla.
! It contains the iterative solvers used in Newton's method. Subroutines:
! - gmresm: GMRESm solver.
! - hookstep: performs hookstep.

module mod_IterativeSolvers

  use mod_Globalvars
  use mod_Matrices
  use mod_ActNewtonSolver

  implicit none
  
contains

subroutine gmresm(time_stepping_sub, m, n, itmax, tol, b, ubest, gmres_iters)

  !----------------------------------------------------------------------
  ! This solver was taken from https://openpipeflow.org/images/d/df/GMRESm.f90
  ! All credits to Willis, A. 
  ! See Willis, A. (2017) SoftwareX 6, 124-127.
  ! https://doi.org/10.1016/j.softx.2017.05.003
  !----------------------------------------------------------------------
  ! Solves Ax = b using the GMRES algorithm with reinitialisation.
  ! requires lapack routines dgelsy, dgesvd.
  !----------------------------------------------------------------------
  ! Inputs
  !   m            - Integer: GMRES dimension
  !   n            - Integer: The size of the system (number of unknowns).
  !   itmax        - Integer: The maximum number of iterations allowed.
  !   tol          - Double precision: The convergence criterion.
  !   b            - Double precision array: The right-hand-side vector.
  ! Outputs
  !   ubest        - Double precision array: Solution that fullfils criterion.
  !   gmres_iters  - Amounts of iterations until reaching criterion
  ! Work Arrays:
  !   h   Hessian matrix,  size (m+1)*m
  !   v   Krylov subspace, size n*(m+1)
  !   w, z, y, p   - double precision arrays: Working arrays for the gmres algorithm.
  !   work, piv    - double precision and integer arrays: for the least-squares algorithm.
  !----------------------------------------------------------------------

  implicit none

  procedure(), pointer, intent(in) :: time_stepping_sub

  integer, intent(in) :: m, n

  integer, intent(in) :: itmax

  double precision, intent(in) :: tol

  double precision, dimension(n), intent(in) :: b

  double precision, dimension(n), intent(out) :: ubest

  integer, intent(out) :: gmres_iters

  double precision, dimension(n) :: x

  double precision, dimension(m+1, m) :: h, h_

  double precision, dimension(n, m+1) :: v

  double precision :: res, res_, stgn

  double precision, dimension(n) :: w, z

  double precision, dimension(m + 1) :: y, p

  double precision, dimension(4*m+1) :: work  ! see LAPACK user manual for DGELSY

  integer :: its, rank, i, j

  integer, dimension(m) :: piv

  double precision, save :: beta
  
  logical :: done

  ! Initialize parameters
  its = 0
  x = 0.0
  v = 0.0

  outer_loop: do while (.true.)

    ! Initialize residual and stagnation values
    res_ = 1d99
    stgn = 1d0 - 1d-14

    ! Initialize beta and the Krylov basis

    call dot(n, x, x, beta)
    beta = sqrt(beta)

    if (beta /= 0.0) then
      call subA(time_stepping_sub, x, w)
    else
      w = 0.0
    end if

    w = b - w

    call dot(n, w, w, beta)
    beta = sqrt(beta)

    v(:, 1) = w / beta

    ! Initialize Hessenberg matrix h to zero
    h = 0.0

    inner_loop: do j = 1, m
      its = its + 1
      z = v(:, j)

      call subA(time_stepping_sub, z, w)
      
      ! Orthogonalize w against previous Krylov vectors
      do i = 1, j
        call dot(n, w, v(:, i), h(i, j))
        w = w - h(i, j) * v(:, i)
      end do

      call dot(n, w, w, h(j+1, j))
      h(j+1, j) = sqrt(h(j+1, j))
      
      v(:, j+1) = w / h(j+1, j)

      ! Solve the least squares problem using LAPACK's dgelsy
      p(1) = beta
      p(2:j+1) = 0.0
      h_(1:j+1, 1:j) = h(1:j+1, 1:j)

      call dgelsy(j+1, j, 1, h_, m+1, p, m+1, piv, m, rank, work, 4*m+1, i)
      ! dgelsy arguments:
      !   j+1, j   top-left block size of h_ used as the coeff. matrix
      !   1        num. of columns in the RHS (simultaneous problems)
      !   h_       coeff. matrix. Will be overwritten on exit
      !   m+1      total rows in h_; coeff. matrix is a block within h_
      !   p        RHS vector. The solution will be stored in p on exit
      !   m+1      total rows in p; the RHS is the j+1 top block of p
      !   piv      vector storing the permuitations in h_ on exit
      !   m        parameter related to rank determination threshold
      !   rank     rank of the coeff. matrix
      !   work     temporary workspace allocated for dgelsy calculations
      !   4*m+1    length of the work vector
      !   i        exit flag; 0 means successful

      if (i /= 0) stop 'gmresm: dgelsy failed'
      y = p

      ! Calculate the residual
      p(1:j+1) = - matmul(h(1:j+1, 1:j), y(1:j))
      p(1) = p(1) + beta

      res = sqrt(p(1:j+1) .dot. p(1:j+1))

      if (res > 1.0e5) then
        print*, "Residual too big, stopping simulation..."
        stop
      end if

      print*, 'gmresm: iteration = ', its, ' residual = ', res

      ! done = (res <= tol .or. its == imx .or. res > res_) ! From original routine by Willis, A.
      done = (res <= tol .or. its == itmax)

      if (done .or. j == m) then

        z = matmul(v(:, 1:j), y(1:j))
        x = x + z

        if (done) exit outer_loop  ! Exit the outer loop as well

        exit inner_loop  ! Exit only the inner loop if needed

      end if

      res_ = res * stgn

    end do inner_loop

  end do outer_loop

  gmres_iters = its

  ubest = x

end subroutine gmresm

!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
!-----------------------------------------------------------------
subroutine hookstep(j, h, m, beta, del, y, info)
  !----------------------------------------------------------------------
  ! This solver was taken from https://openpipeflow.org/images/d/df/GMRESm.f90
  ! All credits to Willis, A. 
  ! See Willis, A. (2017) SoftwareX 6, 124-127.
  ! https://doi.org/10.1016/j.softx.2017.05.003
  !----------------------------------------------------------------------

  ! Generate the hookstep as per Viswanath (2008).

  implicit none

  integer, intent(in)    :: j, m
  double precision, intent(in)   :: h(m+1, j), beta
  double precision, intent(inout) :: del
  double precision, intent(out)  :: y(j)
  double precision :: a(j+1, j), s(j), u(j+1, j+1), vt(j, j), work(5*(j+1))
  double precision :: p(j+1), q(j), mu, qn
  integer, intent(in) :: info

  a = h(1:j+1, 1:j)
  
  call dgesvd('A', 'A', j+1, j, a, j+1, s, u, j+1, vt, j, work, 5*(j+1), info)
  if (info /= 0) stop 'hookstep: dgesvd failed'

  p(1:j) = beta * u(1, 1:j)

  mu = max(s(j)*s(j)*1.0e-6, 1.0e-99)
  qn = 1.0e99

  do while (qn > del)
    mu = mu * 1.1

    q = p(1:j) * s / (mu + s * s)

    call dot(j, q, q, qn)
    qn = sqrt(qn)
  end do

  y = matmul(q, vt)
  
  p = - matmul(h(1:j+1, 1:j), y(1:j))
  p(1) = p(1) + beta
  
  call dot(j+1, p, p, del)
  del = sqrt(del)
end subroutine hookstep

subroutine dot(n, x, y, s)

  !--------------------------------------------------------------------
  !  Subroutine: dot
  !  Purpose: Computes the dot product of two vectors and normalizes it.
  !  Inputs:
  !    n  - Integer: Size of the vectors.
  !    x  - double precision array: The first vector.
  !    y  - double precision array: The second vector.
  !  Output:
  !    s  - double precision: The normalized dot product of x and y.
  !--------------------------------------------------------------------

  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n)
  double precision, intent(out) :: s
  integer :: i

  s = 0.0

  do i = 1, n
    s = s + x(i) * y(i)
  end do

  s = s / n

end subroutine dot

! Subroutine: suba - Performs y = A*x where A is the matrix defined above
subroutine suba_test(n, x, y)

  implicit none

  integer, intent(in) :: n
  
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y

  ! Define the matrix A (3x3 for this example)
  double precision, dimension(3, 3) :: A

  ! Initialize the matrix A
  A(1, :) = (/ 4.0d0, 1.0d0, 0.0d0 /)
  A(2, :) = (/ 1.0d0, 3.0d0, 1.0d0 /)
  A(3, :) = (/ 0.0d0, 1.0d0, 2.0d0 /)

  ! Perform the matrix-vector multiplication y = A*x
  y(:) = 0.0d0
  y(1) = A(1, 1) * x(1) + A(1, 2) * x(2) + A(1, 3) * x(3)
  y(2) = A(2, 1) * x(1) + A(2, 2) * x(2) + A(2, 3) * x(3)
  y(3) = A(3, 1) * x(1) + A(3, 2) * x(2) + A(3, 3) * x(3)

end subroutine suba_test

end module mod_IterativeSolvers
