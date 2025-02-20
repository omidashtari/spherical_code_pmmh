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

  integer :: k, l, m, lm, m_idx

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

    m_idx = 0

    do m = 0, MM*mres, mres
      lm = shtns_lmidx(shtns_c, m, m)

      call SH_to_spat_ml(shtns_c, m_idx, Slm_r(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Ur_m(k, :, m_idx) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0

      call SH_to_spat_ml(shtns_c, m_idx, Slm_t(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Ut_m(k, :, m_idx) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0

      call SH_to_spat_ml(shtns_c, m_idx, Slm_p(lm:lm + (LL + 1 - m)), Sm, LL+1)
      Up_m(k, :, m_idx) = Sm * (1. + dble(min(1, max(0, m)))) ! We need to multiply Fourier by 2 for every m >= 0

      m_idx = m_idx + 1
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

subroutine writeRestart(step, Ra, Ek)

  implicit none

  integer, optional, intent(in) :: step
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
  else if (present(Ek)) then
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
      if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
      end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
    open(unit=11, file=trim(directory)//"/Restart"//trim(adjustl(suffix))//".b", form="unformatted")
  else if (present(step)) then
    write(nom, "(I10)") step
    suffix = "_" // trim(adjustl(nom))
    open(unit=11, file=trim(directory)//"/Restart"//trim(adjustl(suffix))//".b", form="unformatted")
  else
    open(unit=11, file=trim(directory)//"/Restart.b", form="unformatted")
  end if

  write(11) E
  write(11) F
  write(11) T
  if (present(step)) then
    write(11) time
    write(11) step
  end if
  close(11)

end subroutine writeRestart

!-------------------------------------------------------------------------------------

subroutine writeDim(step, Ra, Ek)

  implicit none

  integer, optional, intent(in) :: step
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
  else if (present(Ek)) then
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
      if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
      end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
    open(unit=12, file=trim(directory)//"/Dim"//trim(adjustl(suffix))//".b", form="unformatted")
  else if (present(step)) then
    write(nom, "(I10)") step
    suffix = "_" // trim(adjustl(nom))
    open(unit=12, file=trim(directory)//"/Dim"//trim(adjustl(suffix))//".b", form="unformatted")
  else
    open(unit=12, file=trim(directory)//"/Dim.b", form="unformatted")
  end if

  write(12) KK
  write(12) LL
  write(12) MM
  write(12) shtns%nlm
  write(12) mres
  close(12)
end subroutine writeDim

!-------------------------------------------------------------------------------------

subroutine Output_files(Ur, Up, T_real, step, Ra, Ek, ur_SH, ur_Cheb)

  double precision, dimension(kN, lN, mN), intent(in) :: Ur, Up, T_real
  
  ! Optional arguments for output filename
  integer, optional, intent(in) :: step
  double precision, optional, intent(in) :: Ra
  double precision, optional, intent(in) :: Ek
  ! Optional arguments for underresolution checks
  logical, optional, intent(inout) :: ur_SH, ur_Cheb

  double precision, dimension(0:MM) :: Ekin_spec
  double precision :: T1_mod, TKK2_mod ! Auxiliary variables

  character(len = 10) :: file
  character(len = 10) :: nom
  character(len = 50) :: suffix
  integer :: Ra_int  ! To store the rounded Ra as integer
  character(len = 20) :: Ek_str  ! To store Ek in scientific notation
  integer :: i, k, l, m, merid_plane

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

  if ((.not. present(Ra)) .and. (.not. present(Ek)) .and.(.not. present(step))) then
    suffix = ""
  end if

  open(unit=14, file=trim(directory)//"/Data_eq_plane"//trim(adjustl(suffix))//".dat", form="formatted")
  open(unit=15, file=trim(directory)//"/Data_mer_plane"//trim(adjustl(suffix))//".dat"  , form="formatted")
  open(unit=16, file=trim(directory)//"/Data_mid_gap"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=17, file=trim(directory)//"/Ekin_spectral"//trim(adjustl(suffix))//".dat" , form="formatted")
  open(unit=18, file=trim(directory)//"/T_k_spectral"//trim(adjustl(suffix))//".dat" , form="formatted")

  ! Write the header
  write(14, "(A10, 3x, A16, 3x, A16, 3x, A16)") "Ur", "Ut", "Up", "T"
  write(15, "(A10, 3x, A16, 3x, A16, 3x, A16)") "Ur", "Ut", "Up", "T"
  write(16, "(A10, 3x, A16, 3x, A16, 3x, A16)") "Ur", "Ut", "Up", "T"
  write(17, '(A5, 3x, A15)') "m", "Ekin_spec"
  write(18, '(A5, 3x, A12)') "k", "T_k"

  merid_plane = maxloc(abs(Up(kN/2, lN/2, :) / SinTh(lN/2)), dim=1) ! Choosing a meridional plane

  do k = 1, kN
    do m = 1, mN
      write(14, "(E16.6, 3x, E16.6, 3x, E16.6, 3x, E16.6)") Ur(k, int(lN/2), m), Ut(k, int(lN/2), m) / SinTh(int(lN/2)), & 
            & Up(k, int(lN/2), m) / SinTh(int(lN/2)), T_real(k, int(lN/2), m) ! We recall Ut and Up are multiplied by sin(theta)
    end do
    do l = 1, lN
      write(15, "(E16.6, 3x, E16.6, 3x, E16.6, 3x, E16.6)") Ur(k, l, merid_plane), Ut(k, l, merid_plane) / SinTh(l), &
            & Up(k, l, merid_plane) / SinTh(l), T_real(k, l, merid_plane) ! We recall Ut and Up are multiplied by sin(theta)
    end do
  end do

  do l = 1, lN
    do m = 1, mN
      write(16, "(E16.6, 3x, E16.6, 3x, E16.6, 3x, E16.6)") Ur(int(kN/2), l, m), Ut(int(kN/2), l, m)  / SinTh(l), &
            & Up(int(kN/2), l, m)  / SinTh(l),  T_real(int(kN/2), l, m)  ! We recall Ut and Up are multiplied by sin(theta)
    end do
  end do

  ! Compute spectral Kinetic Energy
  call comp_spectral_KE(Ekin_spec)

  do m = 0, MM*mres, mres
    write(17, '(I5, 3x, E16.6)') m, Ekin_spec(m / mres)
  end do

  ! Check for underresolution in SH for continuation in Ekman
  if (present(ur_SH) .and. (sqrt(Ekin_spec(MM) / maxval(Ekin_spec)) > gr_threshold)) then
    print*, "This is sqrt(Ekin_spec(MM) / maxval(Ekin_spec)) = ", sqrt(Ekin_spec(MM) / maxval(Ekin_spec))
    ur_SH = .true.
  end if

  do k = 1, KK2
    write(18, '(I5, 3x, E16.6)') k, dot_product(real(T(k, :)), real(T(k, :))) + dot_product(aimag(T(k, :)), aimag(T(k, :)))
  end do

  ! Check for underresolution in Chebyshev for continuation in Ekman
  if (present(ur_Cheb)) then
    T1_mod = dot_product(real(T(1, :)), real(T(1, :))) + dot_product(aimag(T(1, :)), aimag(T(1, :)))
    TKK2_mod = dot_product(real(T(KK2, :)), real(T(KK2, :))) + dot_product(aimag(T(KK2, :)), aimag(T(KK2, :)))
    if (sqrt(TKK2_mod / T1_mod) > gr_threshold) then
      print*, "This is sqrt(TKK2_mod / T1_mod) = ", sqrt(TKK2_mod / T1_mod)
      ur_Cheb = .true.
    end if
  end if

  close(14)
  close(15)
  close(16)
  close(17)
  close(17)
  close(18)

end subroutine Output_files

subroutine output_coordinates(step)

  integer :: k, l, m
  integer, optional :: step
  character(len = 10) :: suffix, nom

  if (present(step)) then
    write(nom, "(I10)") step
    suffix = "_" // trim(adjustl(nom))
    open(unit=11, file=trim(directory)//"/r"//trim(adjustl(suffix))//".dat", form="formatted")
    open(unit=12, file=trim(directory)//"/theta"//trim(adjustl(suffix))//".dat", form="formatted")
    open(unit=13, file=trim(directory)//"/phi"//trim(adjustl(suffix))//".dat", form="formatted")
  else
    open(unit=11, file=trim(directory)//"/r.dat", form="formatted")
    open(unit=12, file=trim(directory)//"/theta.dat", form="formatted")
    open(unit=13, file=trim(directory)//"/phi.dat", form="formatted")
  end if

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
