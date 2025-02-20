module mod_NewtonSolver

    use iso_c_binding

    use mod_GlobalVars
    use mod_PrecompSH
    use mod_read
    use mod_ActNewtonSolver
    use mod_IterativeSolvers
    use mod_TimeStep
    use mod_Output
    use mod_tools

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

    double precision :: Ra_tilde ! For continuation in Ek

    integer :: count, newt_steps, gmres_its

    integer :: ios
    
    logical ::  condition

    logical :: final_flag ! Flag to check if we have gone over the limit

    character(len=6) :: cont_type

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

    ! Open file for continuation, checking if file already exists
    open(51, file=trim(directory)//"/Continuation_params.dat", status='unknown', action='read', iostat=ios)
    close(51)
    if (ios /= 0) then 
        ! File does not exist
        open(51,file=trim(directory)//"/Continuation_params.dat", status='unknown', position='append')
        write(51, "(A10,3x,A10,3x,A11,3x,A19,3x,A24,3x,A25,3x,A22)") "Cont. step", "Newt steps", & 
            & "Total GMRES", "Cont. param", "Ekin", "Drift f", "Ur"
    else
        ! File exists
        open(51,file=trim(directory)//"/Continuation_params.dat", status='unknown', position='append')
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

    ! Check continuation parameter and direction of continuation
    if (Ek_final /= 0.) then ! We are doing continuation in Ekman
        if (Ek_final > Ek) then
            cont_type = "Ek_max"
        else 
            cont_type = "Ek_min"
        end if
    else ! We are doing continuation in Rayleigh
        if (Ra_final > Ra) then
            cont_type = "Ra_max"
        else 
            cont_type = "Ra_min"
        end if
    end if

    ! Initialise flags
    max_flag = .false.
    final_flag = .false.
    ur_SH = .false.
    ur_Cheb = .false.
    
    ! Initialise dsmax
    dsmax = 0.

    ! Initialise C_base
    C_base = 0.

    ! Initialise count
    count = 0

    ! Initialise counts for grid refinement
    refinement_count = 0 ! How many times we have done grid refinement
    threshold_count = 0  ! How many times we have gone over the threshold

    ! Compute \tilde{Ra} = Ra_{Th} * EK^{4/3} = Ra_{rot} * Ek^{1/3} * Pr for continuation in Ek
    Ra_tilde = Ra * Ek ** (1. / 3.) * Pr

    ! Set condition according to cont_type
    if (cont_type == "Ra_max") then
        condition = Ra <= Ra_final
    else if (cont_type == "Ra_min") then
        condition = Ra >= Ra_final
    else if (cont_type == "Ek_min") then
        condition = Ek >= Ek_final
    else
        condition = Ek <= Ek_final
    end if

    do while (condition)

        count = count + 1

        print*, "############## Starting continuation step n°", count
        if ((cont_type ==  "Ra_max") .or. (cont_type ==  "Ra_min")) then
            print*, "Rayleigh number in current step Ra = ", Ra
            print*
        else
            print*, "Ekman number in current step Ek = ", Ek
            print*, "Rayleigh number in current step Ra = ", Ra
            print*
        end if

        call newton_solver(NonLinTimeStep_ptr, LinNonLinTimeStep_ptr, C_base, &
                          & newt_steps=newt_steps, gmres_its=gmres_its, cont_type=cont_type)

        if ((cont_type ==  "Ra_max") .or. (cont_type ==  "Ra_min")) then
            write(51,"(I5,8x,I5,10x,I5,7x,E24.16,7x,E24.16,2x,E24.16,5x,E24.16)") count, newt_steps, & 
                    & gmres_its, Ra, Ekin, C_base, maxval(Ur(kN / 2, lN / 2, :))
        else
            write(51,"(I5,8x,I5,10x,I5,7x,E24.16,7x,E24.16,2x,E24.16,5x,E24.16)") count, newt_steps, & 
                    & gmres_its, Ek, Ekin, C_base, maxval(Ur(kN / 2, lN / 2, :))
        end if

        ! Flush the file buffer to ensure data is written
        flush(51)

        if (grid_refine) then
            if (ur_SH .or. ur_Cheb) then ! For continuation in Ek
                if (threshold_count == 2) then 
                    print*, "The fields are underresolved. Increasing resolution..."
                    call grid_refinement()
                    ! Reset flags
                    ur_SH = .false.
                    ur_Cheb = .false.
                    ! Reset count
                    count = 0
                    ! Reset threshold_count
                    threshold_count = 0
                else
                    ! Increase threshold_count
                    threshold_count = threshold_count + 1
                end if
            else
                ! Reset threshold_count
                threshold_count = 0
            end if
        end if

        if (count /= 0) then

            call max_search() ! dsmax and idsmax

            if (count /= 1) print*, "This is S(idsmax) = ", S(idsmax)
            print*, 'This is dsmax = ', dsmax

            if ((cont_type ==  "Ra_max") .or. (cont_type ==  "Ra_min")) then
                dRa = abs(delta_param) / Ra
                print*, 'This is dRa = ', dRa
                print*, 'This is dsmax / dRa = ', dsmax / dRa
                call assign_new_value_Ra(count, newt_steps, cont_type, final_flag)
            else
                dEk = abs(delta_param / log10(Ek))
                print*, 'This is dEk = ', dEk
                print*, 'This is dsmax / dEk = ', dsmax / dEk
                call assign_new_value_Ek(count, newt_steps, cont_type, final_flag)
                Ra = Ra_tilde * Ek ** (- 1. / 3.) / Pr
            end if

            print*
            print*, "############## End continuation step n°", count
            print*

        end if

        ! Update condition
        if (cont_type == "Ra_max") then
            condition = Ra <= Ra_final
        else if (cont_type == "Ra_min") then
            condition = Ra >= Ra_final
        else if (cont_type == "Ek_min") then
            condition = Ek >= Ek_final
        else
            condition = Ek <= Ek_final
        end if

    end do

    print*
    print*, "Path following completed"
    print*

end subroutine continuation_convective_solver

subroutine assign_new_value_Ra(count, newt_steps, cont_type, final_flag)

    implicit none
    
    integer, intent(in) :: count, newt_steps

    character(len=6), intent(in) :: cont_type

    logical, intent(inout) ::  final_flag

    double precision :: a, b, c, x, delta_x

    integer :: i

    if ((dsmax < gamma * dRa) .or. (count <= 3)) then

        if (max_flag) max_flag = .false. ! Reset max_flag
    
        if (count == 1) then
            Ra_nm1 = Ra                                     ! Save Ra from first iteration
            Ra = Ra_nm1 + delta_param                          ! Compute new Ra
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from first iteration
            C_base_nm1 = C_base                             ! Save wavespeed from first iteration
        end if

        if (count == 2) then
            Ra_nm2 = Ra_nm1                                                           ! Save Ra from previous iteration
            Ra_nm1 = Ra                                                               ! Save Ra from current iteration
            if (adapt_param) then
                delta_param = dble(Nopt + 1) / dble(newt_steps + 1) * (Ra_nm1 - Ra_nm2) ! Compute delta_param
            end if
            Ra = Ra_nm1 + delta_param                                                    ! Compute new Ra

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
            if (adapt_param) then
                delta_param = dble(Nopt + 1) / dble(newt_steps + 1) * (Ra_nm1 - Ra_nm2) ! Compute delta_param
            end if
            Ra = Ra_nm1 + delta_param                                                    ! Compute new Ra

            ! Check that we are not going over the limit
            if ((cont_type == "Ra_max") .and. (Ra > Ra_final) .and. (final_flag .eqv. .false.)) then
                Ra = Ra_final
                final_flag = .true.
            else if ((cont_type == "Ra_min") .and. (Ra < Ra_final) .and. (final_flag .eqv. .false.)) then
                Ra = Ra_final
                final_flag = .true.
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

        ! Compute delta_param
        delta_param = Ra - Ra_nm1

        ! Set max_flag
        max_flag = .true.
    end if

end subroutine assign_new_value_Ra

subroutine assign_new_value_Ek(count, newt_steps, cont_type, final_flag)

    implicit none
    
    integer, intent(in) :: count, newt_steps

    character(len=6), intent(in) :: cont_type

    logical, intent(inout) ::  final_flag

    double precision :: a, b, c, x, delta_x

    integer :: i

    if ((dsmax < gamma * dEk) .or. (count <= 3)) then

        if (max_flag) max_flag = .false. ! Reset max_flag
    
        if (count == 1) then
            Ek_nm1 = Ek                                     ! Save Ra from first iteration
            Ek = Ek_nm1 * 10 ** delta_param                    ! Compute new Ek
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from first iteration
            C_base_nm1 = C_base                             ! Save wavespeed from first iteration
        end if

        if (count == 2) then
            Ek_nm2 = Ek_nm1                                                                ! Save Ek from previous iteration
            Ek_nm1 = Ek                                                                    ! Save Ek from current iteration
            if (adapt_param) then
                delta_param = dble(Nopt + 1) / dble(newt_steps + 1) * log10(Ek_nm1 / Ek_nm2) ! Compute delta_param
            end if
            Ek = Ek_nm1 * 10 ** delta_param                                                   ! Compute new Ek

            E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
            C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
            C_base_nm1 = C_base                             ! Save wavespeed from current iteration

            ! Linear interpolation for new state
            a = log10(Ek / Ek_nm2) / log10(Ek_nm1 / Ek_nm2) ! Compute slope
            E_ptr = E_nm2 + a * (E_nm1 - E_nm2)
            F_ptr = F_nm2 + a * (F_nm1 - F_nm2)
            T_ptr = T_nm2 + a * (T_nm1 - T_nm2)
            C_base = C_base_nm2 + a * (C_base_nm1 - C_base_nm2)
        end if

        if (count >= 3) then
            Ek_nm3 = Ek_nm2                                                                ! Save Ek from previous to previous iteration
            Ek_nm2 = Ek_nm1                                                                ! Save Ek from previous iteration
            Ek_nm1 = Ek                                                                    ! Save Ek from current iteration
            if (adapt_param) then
                delta_param = dble(Nopt + 1) / dble(newt_steps + 1) * log10(Ek_nm1 / Ek_nm2) ! Compute delta_param
            end if
            Ek = Ek_nm1 * 10 ** delta_param                                                   ! Compute new Ek

            ! Check that we are not going over the limit
            if ((cont_type == "Ek_max") .and. (Ek > Ek_final) .and. (final_flag .eqv. .false.)) then
                Ek = Ek_final
                final_flag = .true.
            else if ((cont_type == "Ek_min") .and. (Ek < Ek_final) .and. (final_flag .eqv. .false.)) then
                Ek = Ek_final
                final_flag = .true.
            end if

            E_nm3 = E_nm2 ; F_nm3 = F_nm2 ; T_nm3 = T_nm2 ; ! Save state from previous to previous iteration
            E_nm2 = E_nm1 ; F_nm2 = F_nm1 ; T_nm2 = T_nm1 ; ! Save state from previous iteration
            E_nm1 = E_ptr ; F_nm1 = F_ptr ; T_nm1 = T_ptr ; ! Save state from current iteration
            C_base_nm3 = C_base_nm2                         ! Save wavespeed from previous to previous iteration
            C_base_nm2 = C_base_nm1                         ! Save wavespeed from previous iteration
            C_base_nm1 = C_base                             ! Save wavespeed from current iteration

            ! Quadratic extrapolation for new state
            a = log10(Ek / Ek_nm3) * log10(Ek / Ek_nm2) / (log10(Ek_nm1 / Ek_nm3) * log10(Ek_nm1 / Ek_nm2))
            b = log10(Ek / Ek_nm3) * log10(Ek / Ek_nm1) / (log10(Ek_nm2 / Ek_nm3) * log10(Ek_nm2 / Ek_nm1))
            c = log10(Ek / Ek_nm2) * log10(Ek / Ek_nm1) / (log10(Ek_nm3 / Ek_nm2) * log10(Ek_nm3 / Ek_nm1))

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
        Ek_nm3 = Ek_nm2                                 ! Save Ek from previous to previous iteration
        Ek_nm2 = Ek_nm1                                 ! Save Ek from previous iteration
        Ek_nm1 = Ek                                     ! Save Ek from current iteration
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
        Ek = 10 ** log10(Ek_nm1 ** a * Ek_nm2 ** b * Ek_nm3 ** c)

        ! Set S(idsmax)
        S(idsmax) = x

        ! Compute delta_param
        delta_param = log10(Ek / Ek_nm1)

        ! Set max_flag
        max_flag = .true.
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

end subroutine max_search

subroutine newton_solver(NonLinTimeStep_ptr, LinNonLinTimeStep_ptr, C_base, newt_steps, gmres_its, cont_type)

    implicit none

    procedure(), pointer, intent(in) :: NonLinTimeStep_ptr, LinNonLinTimeStep_ptr

    double precision, intent(inout) :: C_base

    integer, optional, intent(out) :: newt_steps, gmres_its

    character(len=6), intent(in), optional :: cont_type

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
        if ((cont_type == "Ra_max") .or. (cont_type == "Ra_min")) then
            call writeRestart(Ra=Ra)
            call writeDim(Ra=Ra)
        else
            call writeRestart(Ek=Ek)
            call writeDim(Ek=Ek)
        end if
    else
        call writeRestart()
        call writeDim()
    end if

    ! Compute velocity components and temperature for output
    call comp_U(E, F, Ur, Ut, Up) ! theta and phi components are multiplied by sin(theta)
    call ToReal(T, T_real, KK2)
    if (present(cont_type)) then
        if ((cont_type == "Ra_max") .or. (cont_type == "Ra_min")) then
            call Output_files(Ur, Up, T_real, Ra=Ra, ur_SH=ur_SH, ur_Cheb=ur_Cheb)
        else
            call Output_files(Ur, Up, T_real, Ek=Ek, ur_SH=ur_SH, ur_Cheb=ur_Cheb)
        end if
    else
        call Output_files(Ur, Up, T_real)
    end if

    call comp_KineticEnergy(Ur, Ut, Up) 
    print*, "Kinetic energy:", Ekin

    print*, 'Simulation completed.'

end subroutine newton_solver

end module mod_NewtonSolver