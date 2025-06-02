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

! This code was developped progressively by Alan Riquier, Camille Rambert 
! and Juan Cruz Gonzalez Sembla at the Laboratoire Physique et mécanique 
! des milieux Hétérogènes (PMMH), under the supervision of Laurette Tuckerman.

program iterative_solver_test
  
  use mod_Globalvars
  use mod_CommandLineParser
  use mod_PrecompSH
  use mod_Matrices
  use mod_Tchebyshev
  use mod_precompXY
  use mod_read
  use mod_init
  use mod_TimeSteppingSolver
  use mod_NewtonSolver
  use mod_TimeStep
  use mod_IterativeSolvers
  use mod_Output

  implicit none

  ! definition of program inputs and parameters
  integer, parameter :: m = 11                   ! the size of the linear system
  character(len = 100) :: coeff_dir = './A.txt'  ! coefficient matrix directory (column-major order)
  character(len = 100) :: rhs_dir = './b.txt'    ! RHS vector directory
  
  ! internal variables
  double precision, dimension(m, m) :: A_     ! the coefficient matrix
  double precision, dimension(m) :: b_        ! the RHS vector
  integer :: unit = 0                         ! unit number for file input and output
  integer :: ios                              ! input/output status


  print*, "################################################################"
  print*, "#                                                              #"
  print*, "#             Testing different iterative solvers              #"
  print*, "#                                                              #"
  print*, "# Solver 1: GMRES                                              #"
  print*, "# Solver 2: BiCGSTAB                                           #"
  print*, "# Solver 3: IDR                                                #"
  print*, "#                                                              #"
  print*, "################################################################"


  ! ---------- read coefficient matrix from file
  open(unit, file = coeff_dir, status = 'old', action = 'read', iostat = ios)
    
  if (ios /= 0) then
    print *, "Error opening file."
    stop
  end if

  read(unit, *) A_
  close(unit)
  unit = unit + 1
  print*, "...Coefficient matrix successfully read from file."

  
  ! ---------- read RHS vector from file
  open(unit, file = rhs_dir, status = 'old', action = 'read', iostat = ios)
    
  if (ios /= 0) then
    print *, "Error opening file."
    stop
  end if

  read(unit, *) b_
  close(unit)
  unit = unit + 1
  print*, "...RHS vector successfully read from file."
  
end program iterative_solver_test
