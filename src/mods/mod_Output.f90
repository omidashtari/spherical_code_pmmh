module mod_Output
  !$ USE OMP_LIB
  use mod_Globalvars
  use mod_Tchebyshev
  use mod_PrecompSH
  use mod_ExplicitTerms
  use mod_Matrices

  implicit none

  double precision :: Ekin ! Kinetic energy

contains

subroutine comp_KineticEnergy(Ur,Ut,Up)

  !##########################################################################
  ! This subroutine compute the kinetic energy of the system provided:
  !    - Ur: the radial velocity
  !    - Ut: the velocity along theta
  !    - Up: the velocity along phi
  ! using the following formula:
  !    Ek = 1/(2*V) * integral_V u² dx
  ! with V the volume of the spherical shell
  !##########################################################################

  implicit none

  double precision, dimension(kN, lN, mN), intent(in) :: Ur, Ut, Up
  double precision, dimension(kN, lN, mN)             :: U2 ! U²
  integer :: k, l, m

  !--- Computing U²
  do k = 1, kN
    do l = 1, lN
      U2(k, l, :) = Ur(k, l, :)**2 + (Ut(k, l, :)**2 + Up(k, l, :)**2) / SinTh(l)**2 ! We recall that the theta and phi components of U are multiplied by sin(theta)
    end do
  end do

  Ekin = 0.

  do k = 1, kN
    do l = 1, lN
      do m = 1, mN
        Ekin =  Ekin + U2(k, l, m) * wts(l) * rN(k) ** 2 / dxdr   &
             &  * pi/kN * 2 * pi/mN * sqrt(1. - xN(k) ** 2)
      end do
    end do
  end do

  Ekin = Ekin / (2 * Vol)

end subroutine comp_KineticEnergy

!-------------------------------------------------------------------------------------

subroutine comp_spectral_KE(Ekin_spec)

  double precision, dimension(0:MM), intent(out) :: Ekin_spec

  double complex, dimension(kN, lN, 0:MM) :: Ur_m, Ut_m, Up_m ! Velocity components in real r and theta and spectral phi
  double complex, dimension(shtns%nlm) :: Slm_r, Slm_t, Slm_p ! Intermediate arrays for SH transform
  double complex, dimension(lN) :: Sm                         ! Intermediate array for SH transform
  double precision, dimension(lN, mN) :: Sh                   ! Intermediate array for SH transform

  double precision :: U2 ! To store modulus of the velocity components

  integer :: k, l, m, lm

  ! Ititialise arrays
  Ur_m = 0. ; Ut_m = 0. ; Up_m = 0. ;
  Slm_r = 0. ; Slm_t = 0. ; Slm_p = 0. ; 
  Sm = 0. ; Sh = 0. ;
  Ekin_spec = 0.

  ! Turn to SH and then Legendre to real
  do k = 1, kN
    Sh = Ur(k, :, :)
    call spat_to_SH(shtns_c, Sh, Slm_r)

    Sh = Ut(k, :, :)
    call spat_to_SH(shtns_c, Sh, Slm_t)

    Sh = Up(k, :, :)
    call spat_to_SH(shtns_c, Sh, Slm_p)

    do m = 0, MM
      lm = shtns_lmidx(shtns_c, m, m)

      call SH_to_spat_ml(shtns_c, m, Slm_r(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Ur_m(k, :, m) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0

      call SH_to_spat_ml(shtns_c, m, Slm_t(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Ut_m(k, :, m) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0

      call SH_to_spat_ml(shtns_c, m, Slm_p(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Up_m(k, :, m) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0
    end do
  end do

  do k = 1, kN
    do l = 1, lN
      do m = 0, MM
        U2 = real(Ur_m(k, l, m))**2 + aimag(Ur_m(k, l, m))**2 + &
                  & (real(Ut_m(k, l, m))**2 + real(Up_m(k, l, m))**2) / SinTh(l)**2 + &
                  & (aimag(Ut_m(k, l, m))**2 + aimag(Up_m(k, l, m))**2) / SinTh(l)**2 ! Modulus
        Ekin_spec(m) =  Ekin_spec(m) + U2 * wts(l) * rN(k) ** 2 / dxdr   &
                      &  * pi/kN * sqrt(1. - xN(k) ** 2) * pi / (2 * Vol) ! Integration in r and theta
      end do
    end do
  end do

end subroutine comp_spectral_KE

!-------------------------------------------------------------------------------------

subroutine writeRestart(Ra, Ek)

  implicit none

  double precision, optional, intent(in) :: Ra
  double precision, optional, intent(in) :: Ek

  integer :: Ra_int  ! To store the rounded Ra as integer
  character(len = 20) :: Ek_str  ! To store Ek in scientific notation

  character(len = 10) :: nom
  character(len = 50) :: suffix

  integer :: i

  if (present(Ra)) then
    Ra_int = int(Ra + 0.5)  ! Round to the nearest integer
    write(nom, "(I0)") Ra_int  ! Use I0 format to remove leading spaces
    suffix = "_Ra_" // trim(adjustl(nom))
    open(unit=11, file=trim(directory)//"/Restart"//trim(adjustl(suffix))//".b", form="unformatted")
  end if

  if (present(Ek)) then
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
      if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
      end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
    open(unit=11, file=trim(directory)//"/Restart"//trim(adjustl(suffix))//".b", form="unformatted")
  end if

  if ((.not. present(Ra)) .and. (.not. present(Ek))) then
    open(unit=11, file=trim(directory)//"/Restart.b", form="unformatted")
  end if

  write(11) E
  write(11) F
  write(11) T
  write(11) time
  write(11) step
  close(11)

end subroutine writeRestart

!-------------------------------------------------------------------------------------

subroutine writeDim(Ra, Ek)

  implicit none

  double precision, optional, intent(in) :: Ra
  double precision, optional, intent(in) :: Ek

  integer :: Ra_int  ! To store the rounded Ra as integer
  character(len = 20) :: Ek_str  ! To store Ek in scientific notation

  character(len = 10) :: nom
  character(len = 50) :: suffix

  integer :: i

  if (present(Ra)) then
    Ra_int = int(Ra + 0.5)  ! Round to the nearest integer
    write(nom, "(I0)") Ra_int  ! Use I0 format to remove leading spaces
    suffix = "_Ra_" // trim(adjustl(nom))
    open(unit=12, file=trim(directory)//"/Dim"//trim(adjustl(suffix))//".b", form="unformatted")
  end if

  if (present(Ek)) then
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
      if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
      end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
    open(unit=12, file=trim(directory)//"/Dim"//trim(adjustl(suffix))//".b", form="unformatted")
  end if

  if ((.not. present(Ra)) .and. (.not. present(Ek))) then
    open(unit=12, file=trim(directory)//"/Dim.b", form="unformatted")
  end if

  write(12) KK
  write(12) LL
  write(12) MM
  write(12) shtns%nlm
  close(12)
end subroutine writeDim

!-------------------------------------------------------------------------------------

subroutine Output_files(Ur, Up, T_real, step, Ra, Ek)

  double precision, dimension(kN, lN, mN), intent(in) :: Ur, Up, T_real
  
  ! Optional arguments for output filename
  integer, optional, intent(in) :: step
  double precision, optional, intent(in) :: Ra
  double precision, optional, intent(in) :: Ek

  double precision, dimension(0:MM) :: Ekin_spec

  character(len = 10) :: file
  character(len = 10) :: nom
  character(len = 50) :: suffix
  integer :: Ra_int  ! To store the rounded Ra as integer
  character(len = 20) :: Ek_str  ! To store Ek in scientific notation
  integer :: i, k, l, m

  ! Determine if we use Ra, Ek, or the step in the file name
  if (present(Ra)) then
    Ra_int = int(Ra + 0.5)  ! Round to the nearest integer
    write(nom, "(I0)") Ra_int  ! Use I0 format to remove leading spaces
    suffix = "_Ra_" // trim(adjustl(nom))
  end if

  if (present(Ek)) then
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
      if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
      end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
  end if

  if (present(step)) then
    write(nom, "(I10)") step
    suffix = "_" // trim(adjustl(nom))
  end if

  open(unit=14, file=trim(directory)//"/Ur_eq_plane"//trim(adjustl(suffix))//".dat", form="formatted")
  open(unit=15, file=trim(directory)//"/T_mid_gap"//trim(adjustl(suffix))//".dat"  , form="formatted")
  open(unit=16, file=trim(directory)//"/Ur_mid_gap"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=17, file=trim(directory)//"/T_r_profile"//trim(adjustl(suffix))//".dat", form="formatted")
  open(unit=18, file=trim(directory)//"/T_eq_plane"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=19, file=trim(directory)//"/Ekin_spectral"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=20, file=trim(directory)//"/Up_mer_plane"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=21, file=trim(directory)//"/T_k_spectral"//trim(adjustl(suffix))//".dat" , form="formatted")

  do k = 1, kN
    write(17,"(E16.6)", advance="no") rN(k)
    write(17,"(E16.6)") T_real(k, int(lN/2), mN)
    do m  = 1, mN
      write(18,"(E16.6)", advance="no") T_real(k, int(lN/2), m)
      write(18,"(A1)", advance="no") ','
      write(14,"(E16.6)", advance="no") Ur(k, int(lN/2), m)
      write(14,"(A1)", advance="no") ','
    end do
    do l = 1, lN
      write(20,"(E16.6)", advance="no") Up(k, l, 1) / SinTh(l) ! We recall Ut and Up are multiplied by sin(theta)
      write(20,"(A1)", advance="no") ','
    end do
    write(14,*)
    write(18,*)
    write(20,*)
  end do

  do l = 1, lN
    do m = 1, mN
      write(15,"(E16.6)", advance="no") T_real(int(kN/2), l, m)
      write(15,"(A1)", advance="no") ','
      write(16,"(E16.6)", advance="no") Ur(int(kN/2), l, m)
      write(16,"(A1)", advance="no") ','
    end do
    write(15,*)
    write(16,*)
  end do

  ! Compute spectral Kinetic Energy
  call comp_spectral_KE(Ekin_spec)

  do m = 0, MM
    write(19, '(I5, 3x, E16.6)') m, Ekin_spec(m)
  end do

  do k = 1, KK2
    write(21, '(I5, 3x, E16.6)') k, dot_product(real(T(k, :)), real(T(k, :))) + dot_product(aimag(T(k, :)), aimag(T(k, :)))
  end do

  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)

end subroutine Output_files

subroutine output_coordinates()

  integer :: k, l, m

  open(unit=11, file=trim(directory)//"/r.dat"    , form="formatted")
  open(unit=12, file=trim(directory)//"/theta.dat", form="formatted")
  open(unit=13, file=trim(directory)//"/phi.dat" , form="formatted")

  do k = 1, kN
    write(11,"(E16.6)") rN(k)
  end do

  do l = 1, lN
    write(12,"(E16.6)") acos(CosTh(l))
  end do

  do m = 1, mN
    write(13,"(E16.6)") phi(m)
  end do

  close(11)
  close(12)
  close(13)

end subroutine output_coordinates

end module mod_Output
