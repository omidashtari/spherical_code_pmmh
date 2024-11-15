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

  allocate(E_p(KKp + 2, nlmp))
  allocate(F_p(KKp + 4, nlmp))
  allocate(T_p(KKp + 2, nlmp))  
  
  open(unit=11, file=trim(directory)//"/"//trim(restart_filename), form="unformatted")
  read(11) E_p
  read(11) F_p
  read(11) T_p

  ! Attempt to read t0 and step_0 from timestepping simulation binary files.
  read(11, IOSTAT=ios) t0
  read(11, IOSTAT=ios) step_0
  
  ! Initialise fields
  E = 0.
  F = 0.
  T = 0.

  print*, 'Allocating spectral fields from restart files...'

  if (shtns%nlm == nlmp) then
    E(1:min(KK, KKp) + 2, :) = E_p(:min(KK, KKp) + 2, :)
    F(1:min(KK, KKp) + 4, :) = F_p(:min(KK, KKp) + 4, :)
    T(1:min(KK, KKp) + 2, :) = T_p(:min(KK, KKp) + 2, :)
  else
    ! We recreate previous grid
    lNp = 0
    mNp = 0
    shtns_c_p = shtns_create(LLp, MMp, mresp, norm)
    call shtns_set_grid_auto(shtns_c_p, layout, eps_polar, 2, lNp, mNp)

    !-- C/Fortran pointer mapping
    call c_f_pointer(cptr=shtns_c_p, fptr=shtns_p)

    ! Tracking and initialising
    i = 1
    do m = 0, MM*mres, mres
      do l = m, LL
        lm = shtns_lmidx(shtns_c, l, m) ! indexing of actual grid
        lmp = shtns_lmidx(shtns_c_p, l, m) ! indexing of previous grid
        E(1:min(KK, KKp) + 2, lm) = E_p(1:min(KK, KKp) + 2, lmp)
        F(1:min(KK, KKp) + 4, lm) = F_p(1:min(KK, KKp) + 4, lmp)
        T(1:min(KK, KKp) + 2, lm) = T_p(1:min(KK, KKp) + 2, lmp)
        i = i + 1
      end do
    end do

    ! Destroy previous grid
    call shtns_unset_grid(shtns_c_p)
    call shtns_destroy(shtns_c_p)

  end if

  print*, 'Allocation completed.'

  close(11)

end subroutine readRestart

end module mod_read
