module mod_read

  use mod_Globalvars
  use mod_PrecompSH
  use mod_Matrices

  implicit none

contains

subroutine readDim()

  implicit none

  integer :: ios  ! Status of the I/O operation

  open(unit=12, file=trim(directory)//"/"//trim(dim_filename),form="unformatted")
  read(12) KKp
  read(12) LLp
  read(12) MMp
  read(12) nlmp

  ! The lines that follows are due to the fact that in previous versions mres could not be controlled by the user.
  ! Attempt to read mresp
  read(12, IOSTAT=ios) mresp
  
  if (ios /= 0) then
    ! If reading mresp fails (ios is non-zero), set a default value for mresp
    mresp = 1
  end if
  
  ! Close the file
  close(12)
  
end subroutine readDim

subroutine readRestart()

  implicit none

  integer :: ios  ! Status of the I/O operation

  integer :: i, l, m, lm, lmp, lNp, mNp
  double complex, dimension(:, :), allocatable :: E_p, F_p, T_p

  type(c_ptr) :: shtns_c_p
  type(shtns_info), pointer :: shtns_p

  allocate(E_p(KKp + 2, nlmp), F_p(KKp + 4, nlmp), T_p(KKp + 2, nlmp))

  if (time_step == "bdf2") then
    allocate(E_tm1_p(KKp + 2, nlmp), F_tm1_p(KKp + 4, nlmp), T_tm1_p(KKp + 2, nlmp))
    allocate(DE_tm1_p(KKp + 2, nlmp), DF_tm1_p(KKp + 4, nlmp), DT_tm1_p(KKp + 2, nlmp))
  end if
  
  open(unit=11, file=trim(directory)//"/"//trim(restart_filename), form="unformatted")
  read(11) E_p
  read(11) F_p
  read(11) T_p

  ! Attempt to read t0 and step_0 from timestepping simulation binary files because they could come from continuation.
  read(11, IOSTAT=ios) t0
  read(11, IOSTAT=ios) step_0

  ! Attempt to read DE_tm1, DF_tm1, DT_tm1, E_tm1, F_tm1 and T_tm1 in case restart file comes from a timestepping simulation using BDF2
  read(11, IOSTAT=ios) DE_tm1_p
  read(11, IOSTAT=ios) DF_tm1_p
  read(11, IOSTAT=ios) DT_tm1_p
  read(11, IOSTAT=ios) E_tm1_p
  read(11, IOSTAT=ios) F_tm1_p
  read(11, IOSTAT=ios) T_tm1_p

  if ((time_step /= "bdf2") .or. ((ios /= 0) .and. time_step == "bdf2")) then
  
    ! Initialise fields
    E = (0., 0.) ; F = (0., 0.) ; T = (0., 0.) ;

    print*, 'Allocating spectral fields from restart files...'

    if (shtns%nlm == nlmp) then
      E(1:min(KK, KKp) + 2, :) = E_p(:min(KK, KKp) + 2, :)
      F(1:min(KK, KKp) + 4, :) = F_p(:min(KK, KKp) + 4, :)
      T(1:min(KK, KKp) + 2, :) = T_p(:min(KK, KKp) + 2, :)
    else
      ! We recreate previous grid
      lNp = 0
      mNp = 0
      shtns_c_p = shtns_create(LLp + 1, MMp, mresp, norm)
      call shtns_set_grid_auto(shtns_c_p, layout, eps_polar, 2, lNp, mNp)

      !-- C/Fortran pointer mapping
      call c_f_pointer(cptr=shtns_c_p, fptr=shtns_p)

      ! Tracking and initialising
      do lmp = 1, nlmp
        l = shtns_lm2l(shtns_c_p, lmp)
        m = shtns_lm2m(shtns_c_p, lmp)
        if ((mod(m, mres) == 0) .and. (m <= MM*mres)) then
          lm = shtns_lmidx(shtns_c, l, m)
          E(1:min(KK, KKp) + 2, lm) = E_p(1:min(KK, KKp) + 2, lmp)
          F(1:min(KK, KKp) + 4, lm) = F_p(1:min(KK, KKp) + 4, lmp)
          T(1:min(KK, KKp) + 2, lm) = T_p(1:min(KK, KKp) + 2, lmp)
        end if
      end do

      ! Destroy previous grid
      call shtns_unset_grid(shtns_c_p)
      call shtns_destroy(shtns_c_p)

    end if
  
  else ! We are restarting a BDF2 with files for the second to last timestep

    ! Initialise fields
    E = (0., 0.) ; F = (0., 0.) ; T = (0., 0.) ;
    DE_tm1 = (0., 0.) ; DF_tm1 = (0., 0.) ; DT_tm1 = (0., 0.) ;
    E_tm1 = (0., 0.) ; F_tm1 = (0., 0.) ; T_tm1 = (0., 0.) ;

    print*, 'Allocating spectral fields from restart files...'

    if (shtns%nlm == nlmp) then
      ! Allocate state
      E(1:min(KK, KKp) + 2, :) = E_p(:min(KK, KKp) + 2, :)
      F(1:min(KK, KKp) + 4, :) = F_p(:min(KK, KKp) + 4, :)
      T(1:min(KK, KKp) + 2, :) = T_p(:min(KK, KKp) + 2, :)

      ! Allocate previous RHS
      DE_tm1(1:min(KK, KKp) + 2, :) = DE_tm1_p(:min(KK, KKp) + 2, :)
      DF_tm1(1:min(KK, KKp) + 4, :) = DF_tm1_p(:min(KK, KKp) + 4, :)
      DT_tm1(1:min(KK, KKp) + 2, :) = DT_tm1_p(:min(KK, KKp) + 2, :)

      ! Allocate previous state
      E_tm1(1:min(KK, KKp) + 2, :) = E_tm1_p(:min(KK, KKp) + 2, :)
      F_tm1(1:min(KK, KKp) + 4, :) = F_tm1_p(:min(KK, KKp) + 4, :)
      T_tm1(1:min(KK, KKp) + 2, :) = T_tm1_p(:min(KK, KKp) + 2, :)
    else
      ! We recreate previous grid
      lNp = 0
      mNp = 0
      shtns_c_p = shtns_create(LLp + 1, MMp, mresp, norm)
      call shtns_set_grid_auto(shtns_c_p, layout, eps_polar, 2, lNp, mNp)

      !-- C/Fortran pointer mapping
      call c_f_pointer(cptr=shtns_c_p, fptr=shtns_p)

      ! Tracking and initialising
      do lmp = 1, nlmp
        l = shtns_lm2l(shtns_c_p, lmp)
        m = shtns_lm2m(shtns_c_p, lmp)
        if ((mod(m, mres) == 0) .and. (m <= MM*mres)) then
          lm = shtns_lmidx(shtns_c, l, m)

          ! Allocate state
          E(1:min(KK, KKp) + 2, lm) = E_p(1:min(KK, KKp) + 2, lmp)
          F(1:min(KK, KKp) + 4, lm) = F_p(1:min(KK, KKp) + 4, lmp)
          T(1:min(KK, KKp) + 2, lm) = T_p(1:min(KK, KKp) + 2, lmp)

          ! Allocate previous RHS
          DE_tm1(1:min(KK, KKp) + 2, lm) = DE_tm1_p(1:min(KK, KKp) + 2, lmp)
          DF_tm1(1:min(KK, KKp) + 4, lm) = DF_tm1_p(1:min(KK, KKp) + 4, lmp)
          DT_tm1(1:min(KK, KKp) + 2, lm) = DT_tm1_p(1:min(KK, KKp) + 2, lmp)

          ! Allocate previous state
          E_tm1(1:min(KK, KKp) + 2, lm) = E_tm1_p(1:min(KK, KKp) + 2, lmp)
          F_tm1(1:min(KK, KKp) + 4, lm) = F_tm1_p(1:min(KK, KKp) + 4, lmp)
          T_tm1(1:min(KK, KKp) + 2, lm) = T_tm1_p(1:min(KK, KKp) + 2, lmp)
        end if
      end do

      ! Destroy previous grid
      call shtns_unset_grid(shtns_c_p)
      call shtns_destroy(shtns_c_p)

    end if

    read_tm1 = .true.

  end if

  deallocate(E_p, F_p, T_p)

  if (time_step == "bdf2") then
    deallocate(E_tm1_p, F_tm1_p, T_tm1_p)
    deallocate(DE_tm1_p, DF_tm1_p, DT_tm1_p)
  end if

  print*, 'Allocation completed.'

  close(11)

end subroutine readRestart

end module mod_read
