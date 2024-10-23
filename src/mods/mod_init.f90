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

    case("Christensen")
      ! Initial temperature corresponds to the article "A numerical dynamo benchmark" by Christensen et al. - 
      ! 2001, Physics of the Earth and planetary interiors
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = Rout*Rin / rN(k) - Rin + 210 * A / sqrt(17920 * pi) &
                          & * (1 - 3 * x(k) ** 2 + 3 * x(k) ** 4 - x(k) ** 6) * SinTh(l) ** 4 &
                          & * cos(4 * phi(m))
          end do
        end do
      end do

    case("sym_2")
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = init_amp * cos(2 * phi(m))
          end do
        end do
      end do

    case("sym_3")
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = init_amp * cos(3 * phi(m))
          end do
        end do
      end do

    case("sym_4")
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = init_amp * cos(4 * phi(m))
          end do
        end do
      end do

    case("sym_5")
      do k = 1, kN
        do l = 1, lN
          do m = 1, mN
            Ti(k, l, m) = init_amp * cos(5 * phi(m))
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
