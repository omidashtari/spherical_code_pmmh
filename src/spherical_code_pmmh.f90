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

program spherical_code_pmmh
  
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


  print*, "################################################################"
  print*, "#                                                              #"
  print*, "#                   Spherical code PMMH                        #"
  print*, "#                                                              #"
  print*, "# This code solves the Boussinesq equations in a rotating      #"
  print*, "# spherical shell using no slip boundary conditions.           #"
  print*, "#                                                              #"
  print*, "# Avalable solvers:                                            #"
  print*, "#   - Timestepping solver (Predictor-corrector,                #"
  print*, "#   Crank-Nicolson, Forward-Backward Euler, BDF2)              #"
  print*, "#   using explicit or implicit Coriolis.                       #"
  print*, "#   - Newton solver using CN or FBE with explicit or           #"
  print*, "#   implicit Coriolis.                                         #"
  print*, "#   - Continuation solver using CN or FBE with explicit or     #"
  print*, "#   implicit Coriolis. Path following performed in Rayleigh    #"
  print*, "#   or Ekman (keeping Ra * Ek ** (1. / 3.) * Pr constant).     #"
  print*, "#                                                              #"
  print*, "# This code uses SHTns for fast SH transforms.                 #"
  print*, "#                                                              #"
  print*, "################################################################"
  print*
  print*, "--------------- Reading the Simulation parameters --------------"
  print*, "#                                                              #"
  call parse_command_line_arguments()
  print*, "#                                                              #"
  print*, "------------------- Precomputing the matrices ------------------"
  print*
  print*, "Precomputing the Chebyshev terms..."
  call PrecompTchebyshev()
  print*, "Precomputing the spherical harmonics grid..."
  call PrecompSH()
  call output_coordinates()

  if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
    ! We call the covective solver with explicit Coriolis
    call convective_solver()
  else if ((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit")) then
    ! We call the Newton covective solver
    call convective_newton_solver()
  else if ((solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) then
    ! We call the continuation covective solver
    call continuation_convective_solver()
  end if

end program spherical_code_pmmh
