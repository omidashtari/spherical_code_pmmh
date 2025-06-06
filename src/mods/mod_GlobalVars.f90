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
! It contains the declaration of global variables and some subroutines:
! - fact: compute the factorial of a number.
! - LinfNormArray: compute the L_ifty norm of a 3D array.
! - LinfNormMat: compute the L_ifty norm of a 2D array.
! - LinfNormVec: compute the L_ifty norm of a 1D array.

!------------------------------- All variables --------------------------------!

module mod_Globalvars
  !$ use OMP_LIB
  use ieee_arithmetic

  implicit none

  double precision, parameter :: pi = 4*atan(1.0) ! Always useful

  !--- Time Stepping parameters

  double precision :: time         ! Current time
  double precision :: delta_t = 0. ! time step
  double precision :: t0 = 0.d0    ! initial time
  integer :: step_0 = 0            ! initial time step
  integer :: NTS = 0               ! Number of Time Steps
  integer :: save_every = 0      ! Save output files every n TS
  integer :: step = 0

  character(len=100) :: solver, time_step, restart, directory, dealiasing, init
  character(len=100) :: restart_filename, dim_filename

  !--- Truncation parameters
  integer :: LL, MM, KK
  integer :: mres = 0
  integer :: kN     ! Number of r points in physical space
  integer :: lN, mN ! Numer of theta and phi points in physical space

  integer :: KK2, KK4 ! These are KK+2 and KK+4

  ! Truncation parameters of restart file
  integer :: KKp, LLp, MMp, nlmp, mresp

  !--- Size of the domain

  double precision :: Rin = 7. / 13.
  double precision :: Rout = 20. / 13.       ! Inner and outer radii

  !--- Prarameters

  double precision :: Pr = 0.   ! Prandtl number
  double precision :: Ek = 0.   ! Ekman number
  double precision :: Ro = 0.   ! Rossby number
  double precision :: Ra = 0.   ! Rayleigh number
  double precision :: IER = -1. ! Ratio beetwen the implicit and explicit steps as defined in eq (21)
                                ! We set it to be -1. to detect if the user has no set it
  double precision :: init_amp  ! Amplitude for symmetric initial condition
  integer          :: sym       ! M azymuthal symmetry for initial condition

  !--- Arrays containing the computed data

  integer, dimension(:), allocatable :: ell    ! Auxiliary array for the solution of the linear system
  
  double complex, dimension(:, :), allocatable, target :: E, F, T ! Spectral arrays for e, f and theta - real part for cosine and imaginary for sine

  double complex, dimension(:, :), allocatable :: DE, DF, DT              ! Spectral arrays for the RHS
  double complex, dimension(:, :), allocatable :: DEp, DFp, DTp           ! Spectral arrays for the RHS - predictor step in case of PC
  double complex, dimension(:, :), allocatable :: DEr, DFr                ! Arrays for the RHS in intermidiate space (real in r spectral in theta an dphi) - corrector step
  double complex, dimension(:, :), allocatable :: DEpr, DFpr              ! Arrays for the RHS in intermidiate space (real in r spectral in theta an dphi) - predictor step

  !--- For BDF2 timestepping scheme
  double complex, dimension(:, :), allocatable :: E_tm1, F_tm1, T_tm1          ! Spectral arrays for e, f and theta at time t-1
  double complex, dimension(:, :), allocatable :: E_tm2, F_tm2, T_tm2          ! Spectral arrays for e, f and theta at time t-2
  double complex, dimension(:, :), allocatable :: DE_tm1, DF_tm1, DT_tm1       ! Spectral arrays for the RHS at time t-1
  double complex, dimension(:, :), allocatable :: DE_tm2, DF_tm2, DT_tm2       ! Spectral arrays for the RHS at time t-2
  double complex, dimension(:, :), allocatable :: E_tm1_p, F_tm1_p, T_tm1_p    ! Spectral arrays to restart e, f and theta at time t-1
  double complex, dimension(:, :), allocatable :: DE_tm1_p, DF_tm1_p, DT_tm1_p ! Spectral arrays to restart the RHS at time t-1
  logical :: read_tm1 = .false. ! Flag to check if we have read files from time t-1

  !--- Arrays to contain the solution of the linear systems
  double complex, dimension(:, :), allocatable :: wE, wF, wT
  double complex, dimension(:, :), allocatable :: wDE, wDF, wDT

  !--- Arrays containing real velocity components and temperature
  double precision, dimension(:, :, :), allocatable :: Ur, Ut, Up, T_real

  !--- Auxiliary arrays to compute U, curl(U) and grad(T)
  double complex, dimension(:, :), allocatable :: Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St

  !--- Intermediate arrays for SH transform
  double precision, dimension(:, :), allocatable :: Sh1, Sh2, Sh3

  !--- Pointers to select RHS computation for timestep: These will be used only if timestepping simulation with CN is selected
  ! or if Newton or Continuation solvers are selected.
  ! Explicit_RHS_ptr => comp_ExplicitRHS if user chooses CN timestep for timestepping simulation
  ! Implicit_RHS_ptr => comp_ImplicitRHS if user chooses CN timestep for timestepping simulation
  ! Explicit_RHS_ptr => comp_RHS_with_rot and Implicit_RHS_ptr => comp_RHS_with_rot if user chooses Newton solver
  procedure(), pointer :: Explicit_RHS_ptr, Implicit_RHS_ptr

  !--- To compute the solve matrices
  double precision, dimension(:, :), allocatable :: Chb0XY, ChbD1XY, ChbD2XY, ChbD4XY, divR

  !--- To impose the boundary conditions in the solve matrices
  double precision, dimension(:, :), allocatable :: BC_mat

  !--- Inverse matrix X⁻¹ and X⁻¹*Y for e, f and Theta:

  double precision, dimension(:, :, :), allocatable :: Xe_inv, Xe_invYe
  double precision, dimension(:, :, :), allocatable :: Xf_inv, Xf_invYf
  double precision, dimension(:, :, :), allocatable :: XT_inv, XT_invYT

  !--- Matrices for implicit Coriolis
  double precision, dimension(:, :, :), allocatable :: Xef, Yef ! Coupled Yef is only used if timestep is CN or PC
  double precision, dimension(:, :, :), allocatable :: Ye_mat, Yf_mat ! These are de-coupled and are only used if timestep is BDF
  integer, dimension(:, :), allocatable :: PIVOT
  double precision, dimension(:, :), allocatable :: A, Ap, A1, DA, DA1

  !--- Newton solver and GMRES
  double complex, dimension(:, :), allocatable, target :: E_base, F_base, T_base ! Spectral arrays to contain base values of U in Newton solver
  double complex, dimension(:, :), pointer :: E_per, F_per, T_per ! Spectral arrays to contain the perturbation of U in Newton solver
  double precision, dimension(:), pointer :: E_ptr, F_ptr, T_ptr ! Pointers for E, F and T
  double precision, dimension(:), pointer :: E_base_ptr, F_base_ptr, T_base_ptr ! Pointers for E_base, F_base and T_base
  integer :: max_newt, max_gmres, restart_gmres, M_wave, wavespeed_loc
  double precision :: newt_eps, newt_delta, tol_gmres
  double precision :: C_base, c_per

  !--- Continuation method
  double precision :: Ek_final = 0. ! Final value in Ekman continuation
  double precision :: Ra_final = 0. ! Final value in Rayleigh continuation 
  double precision, dimension(:), allocatable, target :: E_nm1, F_nm1, T_nm1 ! Arrays to store the state in iteration n minus 1
  double precision, dimension(:), allocatable, target :: E_nm2, F_nm2, T_nm2 ! Arrays to store the state in iteration n minus 2
  double precision, dimension(:), allocatable, target :: E_nm3, F_nm3, T_nm3 ! Arrays to store the state in iteration n minus 3
  double precision, dimension(:), pointer :: S, S_nm1, S_nm2, S_nm3 ! Pointers to perform extrapolation on turning point
  double precision :: Ra_nm1, Ra_nm2, Ra_nm3 ! To store Ra in iteration n minus 1, n minus 2 and n minus 3
  double precision :: Ek_nm1, Ek_nm2, Ek_nm3 ! To store Ra in iteration n minus 1, n minus 2 and n minus 3
  double precision :: C_base_nm1, C_base_nm2, C_base_nm3 ! To store C_base in iteration n minus 1, n minus 2 and n minus 3
  double precision :: dRa, dEk, dsmax ! To store differentials
  double precision :: delta_param = 0. ! To store delta_Ekman or delta_Rayleigh
  double precision :: gamma = 0. ! For turning point detection
  double precision :: gr_threshold = 1. ! Grid refinement threshold
  logical :: adapt_param ! Flag to adapt or not the delta_param using Nopt
  character(len=10) :: adapt_param_string ! To read adapt_param variable from console
  logical :: grid_refine ! Flag to perform grid refinement in continuation
  character(len=10) :: grid_refine_string ! To read grid_refine variable from console
  integer :: Nopt = 0 ! Optimal Newton steps per continuation iteration
  logical :: max_flag ! Flag to compute location of max in E, F and T for turning point detection
  integer :: idsmax, loc_dE, loc_dF, loc_dT ! To store location of max of dE, dF or dT for turning point detection
  logical :: ur_SH, ur_Cheb ! Flags to check if we are underresolved in SH or Chebyshev when doing continuation in Ekman
  integer :: refinement_count ! To store the amount of refinements done throughout continuation in Ekman.
  integer :: threshold_count ! To know how many times in a row the threshold is being overpassed

contains

  function fact(n)

    !***************************************************************************
    ! Knowing n, compute the factorial of n, i.e. n!
    !***************************************************************************

    integer, intent(in) :: n
    double precision :: fact
    integer :: i

    if (n < 0) error stop 'factorial is singular for negative integers'
    fact = 1.0
    do i = 2, n
      fact = fact * i
    end do

  end function fact

  !-----------------------------------------------------------------------------

  subroutine LinfNormArray(X)

    implicit none

    double precision, dimension(:,:,:), intent(in) :: X
    double precision :: Xmax
    integer :: i,j,k,imax,jmax,kmax,N1,N2,N3

    N1 = size(X,1)
    N2 = size(X,2)
    N3 = size(X,3)

    Xmax = 0.d0

    do i = 1,N1
      do j = 1,N2
        do k = 1,N3
          if (ieee_is_nan(X(i,j,k))) then
            print*, 'Found NaN'
            stop
          end if
          if (abs(X(i,j,k)) > Xmax) then
            Xmax = abs(X(i,j,k))
            imax = i
            jmax = j
            kmax = k
          end if
        end do
      end do
    end do

    write(*,*) "Maximum value :", Xmax, "at (i,j,k) =", imax,jmax,kmax, "size=", size(X,1),size(X,2),size(X,3)

  end subroutine LinfNormArray

  !----------------------------------------------------------------------------

  subroutine write_mat(M,imax,jmax,unit)

    integer, intent(in)                                        :: imax,jmax
    double precision, dimension(1:imax,1:jmax), intent(in)     :: M
    integer, intent(in)                                        :: unit

    integer :: i,j

    do j = 1,jmax
      do i = 1,imax
        write(unit, '(*(g0))', advance='no') M(i,j),','
      end do
      write(unit,*) ''
    end do

  end subroutine write_mat

  !-----------------------------------------------------------------------------

  subroutine LinfNormVec(X)

    implicit none

    double precision, dimension(:), intent(in) :: X
    double precision :: Xmax
    integer :: i,imax,N

    Xmax = 0.d0

    N = size(X)

    do i = 1,N
      if (ieee_is_nan(X(i))) then
        print*, 'Found NaN'
        stop
      end if
      if (abs(X(i)) > Xmax) then
        Xmax = abs(X(i))
        imax = i
      end if
    end do

    write(*,*) "Maximum value :", Xmax

  end subroutine LinfNormVec

  !-----------------------------------------------------------------------------

  subroutine LinfNormMat(X)

    implicit none

    double precision, dimension(:,:), intent(in) :: X
    double precision :: Xmax
    integer :: i,j,imax,jmax,N1,N2

    N1 = size(X,1)
    N2 = size(X,2)

    Xmax = 0.d0

    do i = 1,N1
      do j = 1,N2
        if (ieee_is_nan(X(i, j))) then
          print*, 'Found NaN at (i, j) = ', i, j
          stop
        end if
        if (abs(X(i,j)) > Xmax) then
          Xmax = abs(X(i,j))
          imax = i
          jmax = j
        end if
      end do
    end do

    write(*,*) "Maximum value :", Xmax, "at (i,j)=",imax,jmax, "size=",N1,N2

  end subroutine LinfNormMat

end module mod_GlobalVars
