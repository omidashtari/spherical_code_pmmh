module mod_PrecompSH

    use iso_c_binding
    use mod_Globalvars

    implicit none

    include 'shtns.f03' ! Including the interface between our fortran code and SHTns

    !--- shtns_info is a struct that contains the information of our SH transforms
    ! for more information see https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/structshtns__info.html
    type(shtns_info), pointer :: shtns
    type(c_ptr) :: shtns_c

    double precision, pointer :: CosTh(:), SinTh(:)
    double precision, dimension(:), allocatable :: phi
    double precision :: eps_polar                             ! Polar optimisation threshold
    integer :: nthreads, norm, layout, nlm
    double precision, dimension(:), allocatable :: ll1, dphi  ! Arrays for the Laplacian and d/dphi
    double precision, dimension(:), allocatable :: sth_dth    ! Array to perform sin(theta).d/dtheta
    double precision, dimension(:), allocatable :: costh_mul  ! Array to perform multiplication by cos(theta)
    double precision, dimension(:), allocatable :: wts        ! Array to contain the Gauss weights
    double precision, dimension(:, :), allocatable :: SH_norm ! Array to contain the spherical harmonics norm
    integer, dimension(:), allocatable :: idx_even, idx_odd   ! Arrays to contain the odd and even indexes for the mapping of coefficients in the implicit Coriolis case

contains

subroutine PrecompSH()

    implicit none

    integer :: i, l, m, lm, i_even, i_odd

    call shtns_verbose(0)           ! 0 for SHTns not to print anything, see SHTns documentation
    nthreads = shtns_use_threads(0) ! Enables multi-thread transform using OpenMP with num_threads (if available). 

    !--- Polar optimisation threshold polar values of Legendre Polynomials below that threshold are 
    ! neglected (for high m), leading to increased performance (a few percents). 0 = no polar optimization; 
    ! 1.e-14 = VERY safe; 1.e-10 = safe; 1.e-6 = aggresive, but still good accuracy.
    eps_polar = 0.

    !--- Norm and layout for SHTs
    ! For the norm we set the real norm and we omit the 4*pi - 
    ! see https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/spec.html
    norm  = SHT_FOURPI + SHT_REAL_NORM

    !--- For the layout we use Gaussian grid and quadrature - 
    ! see https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/shtns_8h.html#abdccbb8fbce176dbe189d494c94f0f7b
    ! and contiguous latitudes to have A(it, ip) - 
    ! see https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/spat.html
    layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS

    ! The series is truncated at degree LL and order MM*mres - 
    ! see https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/spec.html

    !--- Define the size of the spectral description
    shtns_c = shtns_create(LL + 1, MM, mres, norm)

    !--- Creation of the grid
    if ((dealiasing == 'yes') .or. (dealiasing == 'y')) then
        lN = 0; mN = 0 ! They will be set to optimal values by SHTns
        call shtns_set_grid_auto(shtns_c, layout, eps_polar, 2, lN, mN)
    else
        lN = LL
        mN = 2 * MM
        call shtns_set_grid(shtns_c, layout, eps_polar, lN, mN)
    end if

    !--- C/Fortran pointer mapping
    call c_f_pointer(cptr=shtns_c, fptr=shtns)
    call c_f_pointer(cptr=shtns%ct, fptr=CosTh, shape=[shtns%nlat])
    call c_f_pointer(cptr=shtns%st, fptr=SinTh, shape=[shtns%nlat])

    print*, 'The total number of (l, m) spherical harmonics components is: ', shtns%nlm

    ! Auxiliary array
    ! Allocation
    allocate(ell(shtns%nlm))

    i = 1
    do m=0, MM*mres, mres
        do l = m, LL + 1
            ell(i) = l
            i = i + 1
        end do
    end do

    ! Allocation of spectral arrays
    ! For the poloidal and toroidal field and the spectral temperature field
    allocate(E(KK2, shtns%nlm))
    allocate(F(KK4, shtns%nlm))
    allocate(T(KK2, shtns%nlm))

    ! For the time-stepping
    ! For the RHS
    allocate(DE(KK2, shtns%nlm))
    allocate(DF(KK4, shtns%nlm))
    allocate(DT(KK2, shtns%nlm))

    allocate(DEp(KK2, shtns%nlm))
    allocate(DFp(KK4, shtns%nlm))
    allocate(DTp(KK2, shtns%nlm))

    if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_implicit")) then
        allocate(DEr(KK2, shtns%nlm))
        allocate(DFr(KK4, shtns%nlm))
        allocate(DEpr(KK2, shtns%nlm))
        allocate(DFpr(KK4, shtns%nlm))
    end if

    ! Allocate state and RHS in t-1 if using BDF2
    if (time_step == 'bdf2') then 
        allocate(E_tm1(KK2, shtns%nlm), F_tm1(KK4, shtns%nlm), T_tm1(KK2, shtns%nlm))
        allocate(E_tm2(KK2, shtns%nlm), F_tm2(KK4, shtns%nlm), T_tm2(KK2, shtns%nlm))
        allocate(DE_tm1(KK2, shtns%nlm), DF_tm1(KK4, shtns%nlm), DT_tm1(KK2, shtns%nlm))
        allocate(DE_tm2(KK2, shtns%nlm), DF_tm2(KK4, shtns%nlm), DT_tm2(KK2, shtns%nlm))
    end if

    ! For the Newton solver
    if ((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. & 
        & (solver == "continuation_convective_implicit")) then
        allocate(E_base(KK2, shtns%nlm))
        allocate(F_base(KK4, shtns%nlm))
        allocate(T_base(KK2, shtns%nlm))
    end if

    ! For the solution of the linear system (explicit)
    if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
        & .or. (solver == "continuation_convective_explicit")) then
        allocate(wE(KK2, shtns%nlm), wDE(KK2, shtns%nlm))
        allocate(wF(KK4, shtns%nlm), wDF(KK4, shtns%nlm))
    end if
    allocate(wT(KK2, shtns%nlm), wDT(KK2, shtns%nlm))

    ! For the solution of the linear system (implicit)
    if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_implicit")) then
        allocate(A(2 * LL * (KK2 + KK4), 0:MM), Ap(2 * LL * (KK2 + KK4), 0:MM), A1(2 * LL * (KK2 + KK4), 1))
        allocate(DA(2 * LL * (KK2 + KK4), 0:MM), DA1(2 * LL * (KK2 + KK4), 1))
    end if

    ! To store states if continuation method selected
    if ((solver == "continuation_convective_explicit") .or. &
        & (solver == "continuation_convective_implicit")) then
        allocate(E_nm1(2 * shtns%nlm * KK2), F_nm1(2 * shtns%nlm * KK4), T_nm1(2 * shtns%nlm * KK2))
        allocate(E_nm2(2 * shtns%nlm * KK2), F_nm2(2 * shtns%nlm * KK4), T_nm2(2 * shtns%nlm * KK2))
        allocate(E_nm3(2 * shtns%nlm * KK2), F_nm3(2 * shtns%nlm * KK4), T_nm3(2 * shtns%nlm * KK2))
    end if

    ! Allocate velocity components and temperature
    allocate(Ur(kN, lN, mN), Ut(kN, lN, mN), Up(kN, lN, mN), T_real(kN, lN, mN))

    ! Allocation and computation of spherical harmonics norm
    allocate(SH_norm(0 : MM*mres, 0 : LL))
    do m = 0, MM*mres
        do l = m, LL
            if ((l == 0) .and. (m == 0)) then
                SH_norm(m, l) = 1.
            else if ((m == 0) .and. (l /= 0)) then
                SH_norm(m, l) = sqrt(dble(2 * l + 1))
            else
                SH_norm(m, l) = sqrt(2. * dble(2 * l + 1)) * sqrt(fact(l - m) / fact(l + m))
            end if
        end do
    end do

    ! Allocate matrices for solving time step
    if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
        & .or. (solver == "continuation_convective_explicit")) then
        !--- Inverse matrix X⁻¹ and X⁻¹*Y for e
        allocate(Xe_inv(KK2, KK2, LL + 1), Xe_invYe(KK2, KK2, LL + 1))
        !--- Inverse matrix X⁻¹ and X⁻¹*Y for f
        allocate(Xf_inv(KK4, KK4, LL + 1), Xf_invYf(KK4, KK4, LL + 1))
        !--- Inverse matrix X⁻¹ and X⁻¹*Y for the temperature
        allocate(XT_inv(KK2, KK2, 0:LL + 1), XT_invYT(KK2, KK2, 0:LL + 1)) 
    end if

    if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_implicit")) then
        allocate(Xef(6 * (KK2 + KK4) - 2, 2 * LL * (KK2 + KK4), 0 : MM))
        allocate(PIVOT(2 * LL * (KK2 + KK4), 0 : MM))
        if (time_step == "cn") then
            allocate(Yef(4 * (KK2 + KK4) - 1, 2 * LL * (KK2 + KK4), 0 : MM))
        else if ((time_step == "iee") .or. (time_step == "bdf2")) then
            allocate(Ye_mat(KK2, KK2, LL + 1), Yf_mat(KK4, KK4, LL + 1))
        end if
        !--- Inverse matrix X⁻¹ and X⁻¹*Y for the temperature
        allocate(XT_inv(KK2, KK2, 0:LL + 1), XT_invYT(KK2, KK2, 0:LL + 1))
    end if


    ! Allocation and computation of even and odd indexes for implicit Coriolis
    if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_implicit")) then
        allocate(idx_odd(sum([(((LL + mod(LL, 2)) + 1 - m) / 2, m=0, MM*mres, mres)])))
        allocate(idx_even(sum([(((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2), m=0, MM*mres, mres)])))
        i_even = 1
        i_odd = 1
        i = 1
        do m = 0, MM*mres, mres
            do l = max(m, 1), LL
                lm = shtns_lmidx(shtns_c, l, m)
                if (mod(l, 2) /= 0) then ! odd case
                    idx_odd(i_odd) = lm
                    i_odd = i_odd + 1
                else ! even case
                    idx_even(i_even) = lm
                    i_even = i_even + 1
                end if
                i = i + 1
            end do
        end do
    end if

    ! Allocation and discretisation of phi
    allocate(phi(mN))
    do i = 0, mN-1 
        phi(i + 1) = 2 * pi * i / (mN * mres)
    end do

    ! Allocation and initialisation of Laplacian operator and d/dphi
    allocate(ll1(shtns%nlm), dphi(shtns%nlm))
    ll1 = 0.
    dphi = 0.
    do m = 0, MM*mres, mres
        do l = m, LL + 1
            lm = shtns_lmidx(shtns_c, l, m)
            ll1(lm) = dble(l * (l + 1))
            dphi(lm) = dble(m)
        end do
    end do

    ! Allocation and initialisation of the array for the sin(theta).d/dtheta operation
    allocate(sth_dth(2 * shtns%nlm))
    call st_dt_matrix(shtns_c, sth_dth)

    ! Allocation and initialisation of the array for the multiplication by cos(theta)
    allocate(costh_mul(2 * shtns%nlm))
    call mul_ct_matrix(shtns_c, costh_mul)

    ! Allocation and initialisation of the array containing the Gauss weights
    allocate(wts(shtns%nlat))
    call shtns_gauss_wts(shtns_c, wts(:shtns%nlat / 2)) ! Only fills half of the vector. The second half is a mirror of the first one
    wts(shtns%nlat / 2 + 1:) = wts(shtns%nlat/2:1:-1)
    
end subroutine PrecompSH

end module mod_PrecompSH