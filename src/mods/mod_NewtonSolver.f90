module mod_NewtonSolver

    use iso_c_binding

    use mod_GlobalVars
    use mod_PrecompSH
    use mod_read
    use mod_ActNewtonSolver
    use mod_IterativeSolvers
    use mod_TimeStep
    use mod_Output

    implicit none

contains

subroutine convective_newton_solver()

    implicit none

    procedure(), pointer :: NonLinTimeStep_ptr, LinNonLinTimeStep_ptr

    ! Check for solver type
    if (solver == "newton_convective_explicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xe and Ye matrices..."
        call precompXeYe()
        print*, "Precomputing the Xf and Yf matrices..."
        call precompXfYf()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        NonLinTimeStep_ptr => compute_time_step_convective_explicit_IEE
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit
        Explicit_RHS_ptr => comp_RHS_with_rot

    else if (solver == "newton_convective_implicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xef and Yef matrices..."
        call PrecompimplicitXY()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        NonLinTimeStep_ptr => compute_time_step_convective_implicit_IEE
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit
        Implicit_RHS_ptr => comp_RHS_with_rot

    end if

    print*
    print*, "--------------------- Initial values ---------------------------"
    print*

    print*, "Volume of the shell = ", Vol

    print*, "Reading restart files"
    call readDim()
    call readRestart()

    ! Set pointers for E, F and T
    call c_f_pointer(c_loc(E), E_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F), F_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T), T_ptr, [2 * KK2 * shtns%nlm])

    ! Set pointers for E_base, F_base and T_base
    call c_f_pointer(c_loc(E_base), E_base_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F_base), F_base_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T_base), T_base_ptr, [2 * KK2 * shtns%nlm])

    ! Initialise wavespeed
    C_base = 0.

    ! Now we call the Newton solver
    call newton_solver(NonLinTimeStep_ptr, LinNonLinTimeStep_ptr, C_base)

end subroutine convective_newton_solver

subroutine continuation_convective_solver()
    
    implicit none

    procedure(), pointer :: NonLinTimeStep_ptr, LinNonLinTimeStep_ptr

    double precision :: delta_Ra, Ra_max, Ra_min ! For continuation in Ra

    double precision :: delta_Ek, Ek_max, Ek_min ! For continuation in Ek

    double precision :: Ur_output

    integer :: count, newt_steps, gmres_its
    
    logical ::  condition

    logical :: Ra_max_flag, adapt_Ra  ! For continuation in Ra

    logical :: Ek_min_flag  ! For continuation in Ek

    ! Check for solver type
    if (solver == "continuation_convective_explicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xe and Ye matrices..."
        call precompXeYe()
        print*, "Precomputing the Xf and Yf matrices..."
        call precompXfYf()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        NonLinTimeStep_ptr => compute_time_step_convective_explicit_IEE
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit
        Explicit_RHS_ptr => comp_RHS_with_rot

    else if (solver == "continuation_convective_implicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xef and Yef matrices..."
        call PrecompimplicitXY()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        NonLinTimeStep_ptr => compute_time_step_convective_implicit_IEE
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit
        Implicit_RHS_ptr => comp_RHS_with_rot

    end if

    ! Open file for KE continuation
    open(51,file=trim(directory)//"/Continuation_params.dat", status='unknown', position='append')
    write(51, "(A10,3x,A10,3x,A11,3x,A16,3x,A16,3x,A16,3x,A16)") "Cont. step", "Newt steps", & 
        & "Total GMRES", "Cont. param", "Ekin", "Drift f", "Ur"

    print*
    print*, "--------------------- Initial values ---------------------------"
    print*

    print*, "Volume of the shell = ", Vol

    print*, "Reading restart files"
    call readDim()
    call readRestart()

    ! Set pointers for E, F and T
    call c_f_pointer(c_loc(E), E_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F), F_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T), T_ptr, [2 * KK2 * shtns%nlm])

    ! Set pointers for E_base, F_base and T_base
    call c_f_pointer(c_loc(E_base), E_base_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F_base), F_base_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T_base), T_base_ptr, [2 * KK2 * shtns%nlm])

    ! Set parameters for continuation method
    gamma = 40.

    ! Initialise max_flag
    max_flag = .false.

    ! Set parameters for continuation in Ra
    Ra_max = 140.
    Ra_min = 88.
    delta_Ra = 2.5
    Ra_max_flag = .false.
    adapt_Ra = .true.

    ! Set parameters for continuation in Ek
    Ek_max = 1.0e-2 
    Ek_min = 1.0e-5
    delta_Ek = 1 / (10. ** (1./30.))
    Ek_min_flag = .false.
    
    ! Initialise dsmax
    dsmax = 0.

    ! Initialise C_base
    C_base = 0.

    ! Initialise count
    count = 0

    ! IDEAS ON HOW TO CODE ALL CONTINUATION METHODS IN ONE
    ! We do a condition for the do while => do while (condition)
    ! And that condition will depend on the user's choice for continuation parameter
    ! Then we'll have different assign_new_value functions that will be called using a pointer

    ! Set condition
    condition = Ra <= Ra_max ! For continuation in Ra
    ! condition = Ra >= Ra_min ! For continuation in Ra
    ! condition = Ek >= Ek_min ! For continuation in Ek
    ! condition = Ek <= Ek_max ! For continuation in Ek

    do while (condition)

        count = count + 1

        print*, "############## Starting continuation step n°", count
        print*, "Rayleigh number in current step Ra = ", Ra  ! For continuation in Ra
        print*
        ! print*, "Ekman number in current step Ek = ", Ek  ! For continuation in Ek
        ! print*

        call newton_solver(NonLinTimeStep_ptr, LinNonLinTimeStep_ptr, C_base, &
                          & newt_steps=newt_steps, gmres_its=gmres_its)

        ! Check sign of delta_Ra for Ur output
        if (delta_Ra > 0.) then
            Ur_output = maxval(Ur(kN / 2, lN / 2, :))
        else
            Ur_output = minval(Ur(kN / 2, lN / 2, :))
        end if

        write(51,"(I5,8x,I5,10x,I5,7x,E16.6,7x,E16.6,2x,E16.6,5x,E16.6)") count, newt_steps, & 
                & gmres_its, Ra, Ekin, C_base, Ur_output ! For continuation in Ra
        ! write(51,"(I5,8x,I5,10x,I5,7x,E16.6,7x,E16.6,2x,E16.6,5x,E16.6)") count, newt_steps, & 
        !     & gmres_its, Ek, Ekin, C_base, maxval(Ur(kN / 2, lN / 2, :)) ! For continuation in Ek

        ! Flush the file buffer to ensure data is written
        flush(51)

        call max_search()

        dRa = abs(delta_Ra) / Ra
        ! dEk = abs(delta_Ek) / Ek ! Not checked. 15/10/24

        if (count /= 1) print*, "This is S(idsmax) = ", S(idsmax)
        print*, 'This is dsmax = ', dsmax
        print*, 'This is dRa = ', dRa ! For continuation in Ra

        call assign_new_value_Ra(count, newt_steps, delta_Ra, Ra_max, Ra_max_flag, adapt_Ra) ! For continuation in Ra
        ! call assign_new_value_Ek(count, newt_steps, delta_Ek, Ek_min, Ek_min_flag)         ! For continuation in Ek

        print*
        print*, "############## End continuation step n°", count
        print*

        ! Update condition
        condition = Ra <= Ra_max ! For continuation in Ra
        ! condition = Ra >= Ra_min ! For continuation in Ra
        ! condition = Ek >= Ek_min ! For continuation in Ek
        ! condition = Ek <= Ek_max ! For continuation in Ek

    end do

    print*
    print*, "Path following completed"
    print*

end subroutine continuation_convective_solver

subroutine assign_new_value_Ra(count, newt_steps, delta_Ra, Ra_max, Ra_max_flag, adapt_Ra)

    implicit none
    
    integer, intent(in) :: count, newt_steps

    double precision, intent(inout) :: delta_Ra ! This may not be necessary as input

    double precision, intent(in) :: Ra_max ! This may not be necessary as input

    logical, intent(inout) ::  Ra_max_flag

    logical, intent(in) :: adapt_Ra

    double precision :: a, b, c, x, delta_x

    integer :: N_opt, i

    N_opt = 6

    if ((dsmax < gamma * dRa) .or. (count <= 3)) then

        if (max_flag) max_flag = .false. ! Reset max_flag
    
        if (count == 1) then
            Ra_nm1 = Ra                                     ! Save Ra from first iteration
            Ra = Ra_nm1 + delta_Ra                          ! Compute new Ra
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from first iteration
            C_base_nm1 = C_base                             ! Save wavespeed from first iteration
        end if

        if (count == 2) then
            Ra_nm2 = Ra_nm1                                                           ! Save Ra from previous iteration
            Ra_nm1 = Ra                                                               ! Save Ra from current iteration
            if (adapt_Ra) then
                delta_Ra = dble(N_opt + 1) / dble(newt_steps + 1) * (Ra_nm1 - Ra_nm2) ! Compute delta_Ra
            end if
            Ra = Ra_nm1 + delta_Ra                                                    ! Compute new Ra

            E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
            C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
            C_base_nm1 = C_base                             ! Save wavespeed from current iteration

            ! Linear interpolation for new state
            a = (Ra - Ra_nm2) / (Ra_nm1 - Ra_nm2) ! Compute slope
            E_ptr = E_nm2 + a * (E_nm1 - E_nm2)
            F_ptr = F_nm2 + a * (F_nm1 - F_nm2)
            T_ptr = T_nm2 + a * (T_nm1 - T_nm2)
            C_base = C_base_nm2 + a * (C_base_nm1 - C_base_nm2)
        end if

        if (count >= 3) then
            Ra_nm3 = Ra_nm2                                                           ! Save Ra from previous to previous iteration
            Ra_nm2 = Ra_nm1                                                           ! Save Ra from previous iteration
            Ra_nm1 = Ra                                                               ! Save Ra from current iteration
            if (adapt_Ra) then
                delta_Ra = dble(N_opt + 1) / dble(newt_steps + 1) * (Ra_nm1 - Ra_nm2) ! Compute delta_Ra
            end if
            Ra = Ra_nm1 + delta_Ra                                                    ! Compute new Ra

            ! Check that we are not going over the limit
            if ((Ra > Ra_max) .and. (Ra_max_flag .eqv. .false.)) then
                Ra = Ra_max
                Ra_max_flag = .true.
            end if

            E_nm3 = E_nm2 ; F_nm3 = F_nm2 ; T_nm3 = T_nm2 ; ! Save state from previous to previous iteration
            E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
            C_base_nm3 = C_base_nm2                         ! Save wavespeed from previous to previous iteration
            C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
            C_base_nm1 = C_base                             ! Save wavespeed from current iteration

            ! Quadratic extrapolation for new state
            a = (Ra - Ra_nm3) * (Ra - Ra_nm2) / ((Ra_nm1 - Ra_nm3) * (Ra_nm1 - Ra_nm2))
            b = (Ra - Ra_nm3) * (Ra - Ra_nm1) / ((Ra_nm2 - Ra_nm3) * (Ra_nm2 - Ra_nm1))
            c = (Ra - Ra_nm2) * (Ra - Ra_nm1) / ((Ra_nm3 - Ra_nm2) * (Ra_nm3 - Ra_nm1))

            E_ptr = a * E_nm1 + b * E_nm2 + c * E_nm3
            F_ptr = a * F_nm1 + b * F_nm2 + c * F_nm3
            T_ptr = a * T_nm1 + b * T_nm2 + c * T_nm3
            C_base = a * C_base_nm1 + b * C_base_nm2 + c * C_base_nm3
        end if
    else ! We have detected a turning point, we will extrapolate using S(idsmax)
        print*, '###################################'
        print*, '####### Extrapolating in Ui #######'
        print*, '###################################'

        ! First we compute the new guess for S(idsmax)
        delta_x = S(idsmax) - S_nm1(idsmax)
        x = S(idsmax) + delta_x

        ! Saving
        Ra_nm3 = Ra_nm2                                 ! Save Ra from previous to previous iteration
        Ra_nm2 = Ra_nm1                                 ! Save Ra from previous iteration
        Ra_nm1 = Ra                                     ! Save Ra from current iteration
        E_nm3 = E_nm2 ; F_nm3 = F_nm2 ; T_nm3 = T_nm2 ; ! Save state from previous to previous iteration
        E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
        E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
        C_base_nm3 = C_base_nm2                         ! Save wavespeed from previous to previous iteration
        C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
        C_base_nm1 = C_base                             ! Save wavespeed from current iteration

        ! Compute the extrapolation coefficients
        a = (x - S_nm3(idsmax)) * (x - S_nm2(idsmax)) / ((S_nm1(idsmax) - S_nm3(idsmax)) * (S_nm1(idsmax) - S_nm2(idsmax)))
        b = (x - S_nm3(idsmax)) * (x - S_nm1(idsmax)) / ((S_nm2(idsmax) - S_nm3(idsmax)) * (S_nm2(idsmax) - S_nm1(idsmax)))
        c = (x - S_nm2(idsmax)) * (x - S_nm1(idsmax)) / ((S_nm3(idsmax) - S_nm2(idsmax)) * (S_nm3(idsmax) - S_nm1(idsmax)))

        ! Perform extrapolation
        E_ptr = a * E_nm1 + b * E_nm2 + c * E_nm3
        F_ptr = a * F_nm1 + b * F_nm2 + c * F_nm3
        T_ptr = a * T_nm1 + b * T_nm2 + c * T_nm3
        C_base = a * C_base_nm1 + b * C_base_nm2 + c * C_base_nm3
        Ra = a * Ra_nm1 + b * Ra_nm2 + c * Ra_nm3

        ! Set S(idsmax)
        S(idsmax) = x

        ! Compute delta_Ra
        delta_Ra = Ra - Ra_nm1

        ! Set max_flag
        max_flag = .true.
    end if

end subroutine assign_new_value_Ra

subroutine assign_new_value_Ek(count, newt_steps, delta_Ek, Ek_min, Ek_min_flag)

    implicit none
    
    integer, intent(in) :: count, newt_steps

    double precision, intent(in) :: delta_Ek ! This may not be necessary as input

    double precision, intent(in) :: Ek_min ! This may not be necessary as input

    logical, intent(inout) ::  Ek_min_flag

    double precision :: a, b, c, x, delta_x

    integer :: N_opt, i

    N_opt = 6
    
    if (count == 1) then
        Ek_nm1 = Ek                                     ! Save Ra from first iteration
        Ek = Ek_nm1 * delta_Ek                          ! Compute new Ek
        E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from first iteration
        C_base_nm1 = C_base                             ! Save wavespeed from first iteration
    end if

    if (count == 2) then
        Ek_nm2 = Ek_nm1                                                     ! Save Ek from previous iteration
        Ek_nm1 = Ek                                                         ! Save Ek from current iteration
        Ek = Ek_nm1 * delta_Ek                                              ! Compute new Ek

        E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
        E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
        C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
        C_base_nm1 = C_base                             ! Save wavespeed from current iteration

        ! Linear interpolation for new state
        a = (Ek - Ek_nm2) / (Ek_nm1 - Ek_nm2) ! Compute slope
        E_ptr = E_nm2 + a * (E_nm1 - E_nm2)
        F_ptr = F_nm2 + a * (F_nm1 - F_nm2)
        T_ptr = T_nm2 + a * (T_nm1 - T_nm2)
        C_base = C_base_nm2 + a * (C_base_nm1 - C_base_nm2)
    end if

    if (count >= 3) then
        Ek_nm3 = Ek_nm2                                                     ! Save Ek from previous to previous iteration
        Ek_nm2 = Ek_nm1                                                     ! Save Ek from previous iteration
        Ek_nm1 = Ek                                                         ! Save Ek from current iteration
        Ek = Ek_nm1 * delta_Ek                                              ! Compute new Ek

        ! Check that we are not going over the limit
        if ((Ek < Ek_min) .and. (Ek_min_flag .eqv. .false.)) then
            Ek = Ek_min
            Ek_min_flag = .true.
        end if

        E_nm3 = E_nm2 ; F_nm3 = F_nm2 ; T_nm3 = T_nm2 ; ! Save state from previous to previous iteration
        E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
        E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
        C_base_nm3 = C_base_nm2                         ! Save wavespeed from previous to previous iteration
        C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
        C_base_nm1 = C_base                             ! Save wavespeed from current iteration

        ! Quadratic extrapolation for new state
        a = (Ek - Ek_nm3) * (Ek - Ek_nm2) / ((Ek_nm1 - Ek_nm3) * (Ek_nm1 - Ek_nm2))
        b = (Ek - Ek_nm3) * (Ek - Ek_nm1) / ((Ek_nm2 - Ek_nm3) * (Ek_nm2 - Ek_nm1))
        c = (Ek - Ek_nm2) * (Ek - Ek_nm1) / ((Ek_nm3 - Ek_nm2) * (Ek_nm3 - Ek_nm1))

        E_ptr = a * E_nm1 + b * E_nm2 + c * E_nm3
        F_ptr = a * F_nm1 + b * F_nm2 + c * F_nm3
        T_ptr = a * T_nm1 + b * T_nm2 + c * T_nm3
        C_base = a * C_base_nm1 + b * C_base_nm2 + c * C_base_nm3
    end if

    ! Now we re-compute the matrices
    if (solver == "continuation_convective_explicit") then
        print*, "Precomputing the Xe and Ye matrices..."
        call precompXeYe()
        print*, "Precomputing the Xf and Yf matrices..."
        call precompXfYf()
    else if (solver == "continuation_convective_implicit") then
        print*, "Precomputing the Xef and Yef matrices..."
        call PrecompimplicitXY()
    end if

end subroutine assign_new_value_Ek

subroutine max_search()

    implicit none

    double precision :: max_E, max_F, max_T

    if (.not. max_flag) then

        ! Find maximum absolute values and their locations
        max_E = maxval(merge(abs((E_ptr - E_nm1) / E_ptr), 0.0, E_ptr > 1e-2))
        loc_dE = maxloc(merge(abs((E_ptr - E_nm1) / E_ptr), 0.0, E_ptr > 1e-2), dim=1)

        max_F = maxval(merge(abs((F_ptr - F_nm1) / F_ptr), 0.0, F_ptr > 1e-2))
        loc_dF = maxloc(merge(abs((F_ptr - F_nm1) / F_ptr), 0.0, F_ptr > 1e-2), dim=1)

        max_T = maxval(merge(abs((T_ptr - T_nm1) / T_ptr), 0.0, T_ptr > 1e-2))
        loc_dT = maxloc(merge(abs((T_ptr - T_nm1) / T_ptr), 0.0, T_ptr > 1e-2), dim=1)

        ! Check max value among them all and set pointers S, S_nm1, S_nm2 and S_nm3
        if (max_E > max(max_F, max_T)) then
            dsmax = max_E
            idsmax = loc_dE
            call c_f_pointer(c_loc(E), S, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(E_nm1), S_nm1, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(E_nm2), S_nm2, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(E_nm3), S_nm3, [2 * KK2 * shtns%nlm])
        else if (max_F > max(max_T, max_E)) then
            dsmax = max_F
            idsmax = loc_dF
            call c_f_pointer(c_loc(F), S, [2 * KK4 * shtns%nlm])
            call c_f_pointer(c_loc(F_nm1), S_nm1, [2 * KK4 * shtns%nlm])
            call c_f_pointer(c_loc(F_nm2), S_nm2, [2 * KK4 * shtns%nlm])
            call c_f_pointer(c_loc(F_nm3), S_nm3, [2 * KK4 * shtns%nlm])
        else if (max_T > max(max_E, max_F)) then
            dsmax = max_T
            idsmax = loc_dT
            call c_f_pointer(c_loc(T), S, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(T_nm1), S_nm1, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(T_nm2), S_nm2, [2 * KK2 * shtns%nlm])
            call c_f_pointer(c_loc(T_nm3), S_nm3, [2 * KK2 * shtns%nlm])
        end if
    else
        dsmax = abs((S(idsmax) - S_nm1(idsmax)) / S(idsmax))
    end if

    ! if (.not. max_flag) print*, "This is S(idsmax) = ", S(idsmax)

end subroutine max_search

subroutine newton_solver(NonLinTimeStep_ptr, LinNonLinTimeStep_ptr, C_base, newt_steps, gmres_its)

    implicit none

    procedure(), pointer, intent(in) :: NonLinTimeStep_ptr, LinNonLinTimeStep_ptr

    double precision, intent(inout) :: C_base

    integer, optional, intent(out) :: newt_steps, gmres_its

    double precision, dimension(2 * shtns%nlm * (2 * KK2 + KK4)) :: EFT_RHS, EFT_best

    double precision :: max_real, norm_FU, norm_EFT_new, norm_EFT_best

    integer :: i_newt, k, j, lm, k_max, lm_max, lm_init, gmres_iters

    double complex, dimension(:), pointer :: E_ptr_cplx, F_ptr_cplx, T_ptr_cplx
    double precision :: max_E, max_F, max_T
    integer :: loc_max_E, loc_max_F, loc_max_T

    call c_f_pointer(c_loc(E), E_ptr_cplx, [KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F), F_ptr_cplx, [KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T), T_ptr_cplx, [KK2 * shtns%nlm])

    ! Compute wavespeed location - search for biggest real part of the k, l, m = M_wave coefficient 
    max_real = 0.
    lm_init = shtns_lmidx(shtns_c, M_wave, M_wave) ! lm index for l = 1, m = M_wave
    lm_max = shtns_lmidx(shtns_c, LL + 1, M_wave)  ! Maximum lm for m = M_wave

    max_E = maxval(real(E_ptr_cplx(KK2 * (lm_init - 1) + 1 : KK2 * lm_max)))
    loc_max_E = maxloc(real(E_ptr_cplx(KK2 * (lm_init - 1) + 1 : KK2 * lm_max)), dim=1)
    max_F = maxval(real(F_ptr_cplx(KK4 * (lm_init - 1) + 1 : KK4 * lm_max)))
    loc_max_F = maxloc(real(F_ptr_cplx(KK4 * (lm_init - 1) + 1 : KK4 * lm_max)), dim=1)
    max_T = maxval(real(T_ptr_cplx(KK2 * (lm_init - 1) + 1 : KK2 * lm_max)))
    loc_max_T = maxloc(real(T_ptr_cplx(KK2 * (lm_init - 1) + 1 : KK2 * lm_max)), dim=1)

    if (max_E > max(max_F, max_T)) then
        wavespeed_loc = 2 * (KK2 * (lm_init - 1) + loc_max_E)
    else if (max_F > max(max_T, max_E)) then
        wavespeed_loc = 2 * shtns%nlm * KK2 + 2 * (KK4 * (lm_init - 1) + loc_max_F)
    else if (max_T > max(max_E, max_F)) then
        wavespeed_loc = 2 * shtns%nlm * (KK2 + KK4) + 2 * (KK2 * (lm_init - 1) + loc_max_T)
    end if

    print*
    print*, "------------------ Beginning Newton iteration ------------------"
    print*

    if (present(gmres_its)) gmres_its = 0

    do i_newt = 1, max_newt

        print*, "############## Starting Newton iteration n°", i_newt
        print*

        ! Compute F(U)
        call act_on_RHS(NonLinTimeStep_ptr, EFT_RHS(1), EFT_RHS(2 * shtns%nlm * KK2 + 1), &
                        & EFT_RHS(2 * shtns%nlm * (KK2 + KK4) + 1))

        ! Compute ||F(U)|| and check if it is lower than newt_eps
        norm_FU = sqrt(EFT_RHS .dot. EFT_RHS)
        print*, 'This is ||F(U)|| = ', norm_FU
        if (norm_FU <= newt_eps) exit

        call gmresm(LinNonLinTimeStep_ptr, restart_gmres, 2 * shtns%nlm * (2 * KK2 + KK4), max_gmres, &
                    & tol_gmres, EFT_RHS, EFT_best, gmres_iters)

        ! Save total amount of GMRES iterations
        if (present(gmres_its)) gmres_its = gmres_its + gmres_iters

        ! Compute the norm of EFT_best (with c_per inside)
        norm_EFT_best = sqrt(EFT_best .dot. EFT_best)

        ! Estract c_per from EFT_best
        c_per = EFT_best(wavespeed_loc)
        EFT_best(wavespeed_loc) = 0.

        ! Update C_base
        C_base = C_base - c_per

        ! Update EFT and compute norm - Here we should include C
        call update_EFT(EFT_best(1), EFT_best(2 * shtns%nlm * KK2 + 1), &
                        & EFT_best(2 * shtns%nlm * (KK2 + KK4) + 1), norm_EFT_new)

        print*, 'This is norm_EFT_best / norm_EFT_new = ', norm_EFT_best / norm_EFT_new

        print*, 'This is the wavespeed: ', C_base

        print*
        print*, "############## End Newton iteration n°", i_newt
        print*
        
        if (norm_EFT_best / norm_EFT_new <= newt_delta) then
            print*, "Newton method did not converge, delta criterion attained"
            print*, "Stopping simulation..."
            stop
        end if

    end do

    print*
    print*, "Exiting Newton solver and saving data"
    print*

    if ((i_newt == max_newt + 1) .and. (norm_FU >= newt_eps)) then
        print*, "Newton method did not converge, maximum number of iterations attained" 
        print*, "Stopping simulation..."
        stop
    end if

    if (present(newt_steps)) newt_steps = i_newt

    ! Saving restart files
    if ((solver == "continuation_convective_explicit") .or. &
        & (solver == "continuation_convective_implicit")) then
        call writeRestart(Ra=Ra) ! For continuation in Ra
        call writeDim(Ra=Ra) ! For continuation in Ra
        ! call writeRestart(Ek=Ek) ! For continuation in Ek
        ! call writeDim(Ek=Ek) ! For continuation in Ek
    else
        call writeRestart()
        call writeDim()
    end if

    ! Compute velocity components and temperature for output
    call comp_U(E, F, Ur, Ut, Up) ! theta and phi components are multiplied by sin(theta)
    call ToReal(T, T_real, KK2)
    call Output_files(Ur, Up, T_real, Ra=Ra) ! For continuation in Ra
    ! call Output_files(Ur, Up, T_real, Ek=Ek) ! For continuation in Ek

    call comp_KineticEnergy(Ur, Ut, Up) 
    print*, "Kinetic energy:", Ekin

    print*, 'Simulation completed.'

end subroutine newton_solver

end module mod_NewtonSolver