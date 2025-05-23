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
! It contains the subroutines to initialise fields for timestepping simulations. 
! Subroutines:
! - init_TEMP: initialise temperature field according to the selected initial condition.
! - init_EF: intialise E and F with 0.

module mod_init
!$ use OMP_LIB
  use mod_Globalvars
  use mod_PrecompSH
  use mod_Matrices
  use mod_Tchebyshev
  use mod_ExplicitTerms

  implicit none

contains

subroutine init_TEMP()

  !##########################################################################
  !                    Initial value of the temperature
  !##########################################################################
  
  implicit none

  double precision, dimension(kN, lN, mN) :: Ti, Ti_test  ! Real initial field
  double precision, dimension(kN) :: x
  double precision :: A = 0.1d0
  integer :: k, l, m, lm

  x = 2 * rN - Rin - Rout

  select case (init)

    case("christensen")
      ! Initial temperature corresponds to the article "A numerical dynamo benchmark" by Christensen et al. - 
      ! 2001, Physics of the Earth and planetary interiors
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = 210 * A / sqrt(17920 * pi) &
                          & * (1 - 3 * x(k) ** 2 + 3 * x(k) ** 4 - x(k) ** 6) * SinTh(l) ** 4 &
                          & * cos(4 * phi(m)) ! + Rout*Rin / rN(k) - Rin ! To add if solving for total temperature
          end do
        end do
      end do

    case("symmetric")
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = init_amp * cos(sym * phi(m))
          end do
        end do
      end do

  end select

  call BackToSpectral(Ti, T, KK2)

  call ToReal(T, Ti_test, KK2)

  write(*,*) "Init Check:"
  call LinfNormArray(Ti - Ti_test)

end subroutine init_TEMP

!----------------------------------------------------------------------------

subroutine init_EF()

  implicit none

  E = 0.

  F = 0.

end subroutine init_EF

end module mod_init
