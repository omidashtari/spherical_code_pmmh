module mod_CommandLineParser

    use mod_Globalvars
    
    implicit none

    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

contains

subroutine parse_command_line_arguments()

    ! Declare local variables
    character(len=100) :: arg, param, param_lc, value
    character(len=1024) :: command_line_string
    character(len=100), dimension(41) :: params_string_array
    integer :: i, lc, ic

    ! Set params_string_array
    params_string_array = [character(len=100) :: &
        "-delta_t", "-nts", "-save_every", "-save_restart", &
        "-solver", "-restart", "-dealiasing", "-directory", &
        "-kk", "-ll", "-mm", "-mres", "-rin", "-rout", "-pr", "-ek", &
        "-ro", "-ra", "-ier", "-save_ur_mgep_t", "-init", &
        "-init_amp", "-sym", "-max_newt", "-max_gmres", "-restart_gmres", &
        "-newt_eps", "-newt_delta", "-tol_gmres", "-m_wave", &
        "-time_step", "-restart_filename", "-dim_filename", &
        "-ek_final", "-ra_final", "-delta_param", "-adapt_param", &
        "-gamma", "-nopt", "-grid_refine", "-gr_threshold" &
    ]

    ! Set strings to default values
    restart = 'no'
    solver = 'convective_explicit'
    time_step = 'no'
    directory = '.'
    dealiasing = 'yes'
    init = 'no'
    restart_filename = "no"
    dim_filename = "no"
    adapt_param_string = "no"

    ! Parse command line arguments
    i = 1
    do while (i <= command_argument_count())
       
        call get_command_argument(i, arg)
       
        ! Extract parameter key
       param = trim(adjustl(arg))
       
       ! Check if the next argument exists
       if (i < command_argument_count()) then

            call get_command_argument(i+1, value)

            ! Turn to lower case
            param_lc = param
            do lc = 1, LEN_TRIM(param)
                ic = INDEX(cap, param(lc:lc))
                if (ic > 0) param_lc(lc:lc) = low(ic:ic)
            end do

            if (any(params_string_array == param_lc) .eqv. .false.) then
                print*, 'Simulation stopped, invalid parameter : ', param_lc
                stop
            end if

            select case (param_lc)

            case ("-delta_t")
                read(value, *) delta_t

            case ("-nts")
                read(value, *) NTS

            case ("-save_every")
                read(value, *) save_every
                
            case ("-save_restart")
                read(value, *) save_restart

            case ("-save_ur_mgep_t")
                read(value, *) save_Ur_mgep_t
            
            case ("-solver")
                read(value, *) solver
                ! Turn to lower case
                do lc = 1, LEN_TRIM(solver)
                    ic = INDEX(cap, solver(lc:lc))
                    if (ic > 0) solver(lc:lc) = low(ic:ic)
                end do

            case ("-time_step")
                read(value, *) time_step
                ! Turn to lower case
                do lc = 1, LEN_TRIM(time_step)
                    ic = INDEX(cap, time_step(lc:lc))
                    if (ic > 0) time_step(lc:lc) = low(ic:ic)
                end do

            case ("-restart")
                read(value, *) restart
                ! Turn to lower case
                do lc = 1, LEN_TRIM(restart)
                    ic = INDEX(cap, restart(lc:lc))
                    if (ic > 0) restart(lc:lc) = low(ic:ic)
                end do

            case ("-restart_filename")
                read(value, *) restart_filename

            case ("-dim_filename")
                read(value, *) dim_filename

            case ("-dealiasing")
                read(value, *) dealiasing
                ! Turn to lower case
                do lc = 1, LEN_TRIM(dealiasing)
                    ic = INDEX(cap, dealiasing(lc:lc))
                    if (ic > 0) dealiasing(lc:lc) = low(ic:ic)
                end do

            case ("-directory")
                directory = value

            case ("-init")
                read(value, *) init
                ! Turn to lower case
                do lc = 1, LEN_TRIM(init)
                    ic = INDEX(cap, init(lc:lc))
                    if (ic > 0) init(lc:lc) = low(ic:ic)
                end do

            case ("-init_amp")
                read(value, *) init_amp

            case ("-sym")
                read(value, *) sym

            case ("-kk")
                read(value, *) KK

            case ("-ll")
                read(value, *) LL

            case ("-mm")
                read(value, *) MM

            case ("-mres")
                read(value, *) mres

            case ("-rin")
                read(value, *) Rin

            case ("-rout")
                read(value, *) Rout

            case ("-pr")
                read(value, *) Pr

            case ("-ek")
                read(value, *) Ek

            case ("-ro")
                read(value, *) Ro

            case ("-ra")
                read(value, *) Ra

            case ("-ek_final")
                read(value, *) Ek_final

            case ("-ra_final")
                read(value, *) Ra_final

            case ("-ier")
                read(value, *) IER

            case ("-max_newt")
                read(value, *) max_newt
            
            case ("-max_gmres")
                read(value, *) max_gmres

            case ("-restart_gmres")
                read(value, *) restart_gmres

            case ("-newt_eps")
                read(value, *) newt_eps

            case ("-newt_delta")
                read(value, *) newt_delta

            case ("-tol_gmres")
                read(value, *) tol_gmres

            case ("-m_wave")
                read(value, *) M_wave

            case ("-adapt_param")
                read(value, *) adapt_param_string
                ! Turn to lower case
                do lc = 1, LEN_TRIM(adapt_param_string)
                    ic = INDEX(cap, adapt_param_string(lc:lc))
                    if (ic > 0) adapt_param_string(lc:lc) = low(ic:ic)
                end do

            case ("-delta_param")
                read(value, *) delta_param

            case ("-nopt")
                read(value, *) Nopt

            case ("-gamma")
                read(value, *) gamma

            case ("-grid_refine")
                read(value, *) grid_refine_string
                ! Turn to lower case
                do lc = 1, LEN_TRIM(grid_refine_string)
                    ic = INDEX(cap, grid_refine_string(lc:lc))
                    if (ic > 0) grid_refine_string(lc:lc) = low(ic:ic)
                end do

            case ("-gr_threshold")
                read(value, *) gr_threshold

          end select
          ! Increment i by 2 to skip both parameter key and value
          i = i + 2

       else

          print *, "Invalid command line argument: ", arg
          stop

       end if

    end do

    call check_arg_validity()

    command_line_string = ""  ! Initialize the string
  
    i = 1

    do while (i <= command_argument_count())
      call get_command_argument(i, value=arg)
      command_line_string = trim(command_line_string) // " " // trim(arg) ! Concatenate with spaces
      i = i + 1
    end do

    ! Create file to save parameters
    open(30,file=trim(directory)//"/parameters.dat", status='unknown', position='append')
    write(30, '(a)') "These are the parameters used in this simulation:"
    write(30, '(a)')
    write(30, '(a)') trim(command_line_string)
    write(30, '(a)')
    close(30)

end subroutine parse_command_line_arguments

subroutine check_arg_validity()

    implicit none

    print*, "Simulation parameters..."

    if (solver == "convective_explicit") then
        print*, "Simulation performed using the convective solver with explicit Coriolis"
    else if (solver == "convective_implicit") then
        print*, "Simulation performed using the convective solver with implicit Coriolis"
    else if (solver == "newton_convective_explicit") then
        print*, "Simulation performed using the Newton convective solver with explicit Coriolis"
    else if (solver == "newton_convective_implicit") then
        print*, "Simulation performed using the Newton convective solver with implicit Coriolis"
    else if (solver == "continuation_convective_explicit") then
        print*, "Simulation performed using the continuation convective solver with explicit Coriolis"
    else if (solver == "continuation_convective_implicit") then
        print*, "Simulation performed using the continuation convective solver with implicit Coriolis"
    else
        print*, "Simulation stopped - solver has an invalid value"
        stop
    end if

    if (((solver == 'convective_explicit') .or. (solver == 'convective_implicit')) .and. &
        & (time_step == 'no')) then
        print*, 'No timestepping scheme selected, switching to predictor-corrector'
        time_step = 'pc'
    else if (((solver == 'newton_convective_explicit') .or. (solver == 'newton_convective_implicit') .or. & 
        (solver == 'continuation_convective_explicit') .or. (solver == 'continuation_convective_implicit')) .and. &
        & (time_step == 'no')) then
        print*, 'No timestepping scheme selected, switching to Forward-Backward Euler'
        time_step = 'fbe'
    else if (((solver == 'newton_convective_explicit') .or. (solver == 'newton_convective_implicit') .or. & 
        (solver == 'continuation_convective_explicit') .or. (solver == 'continuation_convective_implicit')) .and. &
        & ((time_step == 'pc') .or. (time_step == 'bdf2'))) then
        print*, 'Predictor corrector or BDF2 timestep selected but Newton cannot work with them. ' // &
                & 'Switching to Forward-Backward Euler'
        time_step = 'fbe'
    end if

    if (time_step == 'pc') then
        print*, 'Simulation is being executed with predictor corrector time step'
    else if (time_step == 'cn') then
        print*, 'Simulation is being executed with Crank-Nicolson time step'
    else if (time_step == 'fbe') then
        print*, 'Simulation is being executed with Forward-Backward Euler (BDF1) time step'
    else if (time_step == 'bdf2') then
        print*, 'Simulation is being executed with BDF2 time step'
    else
        print*, "Simulation stopped - time_step has an invalid value"
        stop
    end if

    if ((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) then
        if ((restart == "no") .or. (restart == 'n')) then
            print*, "Newton or continuation solver should start from restart file"
            stop
        else if ((init /= "no") .or. (init_amp /= 0.)) then
            print*, "Newton solver should not have initial condition"
            stop
        end if
    end if

    if (((restart == "yes") .or. (restart == 'y')) .and. (init /= "no")) then
        print*, "Initial condition option and restart from file have both been selected, choose just one"
        print*, "Stopping simulation..."
        stop
    else if (((restart == "no") .or. (restart == 'n')) .and. (init == "no")) then
        print*, "Neither initial condition nor restart from file have been selected"
        print*, "Switching to default initial condition"
        init = "christensen"
    end if

    if ((restart == "yes") .or. (restart == "y")) then
        print*, "Beginning simulation from restart file"
    else if ((restart == "no") .or. (restart == "n")) then
        print*, "Beginning simulation from initial condition"
    else
        print*, "Simulation stopped - restart has an invalid value"
        stop
    end if

    if (((restart == "yes") .or. (restart == "y")) .and. (restart_filename /= "no")) then
        print*, "Restaring fields from ", restart_filename
    else if (((restart == "yes") .or. (restart == "y")) .and. (restart_filename == "no")) then
        print*, "Restart option activated but no filename given for field restart."
        print*, "Restarting fields from Restart.b file"
        restart_filename = "Restart.b"
    end if

    if (((restart == "yes") .or. (restart == "y")) .and. (dim_filename /= "no")) then
        print*, "Restaring dimensions from ", dim_filename
    else if (((restart == "yes") .or. (restart == "y")) .and. (dim_filename == "no")) then
        print*, "Restart option activated but no filename given for dimensions restart."
        print*, "Restarting dimension from Dim.b file"
        dim_filename = "Dim.b"
    end if

    if (((restart == "no") .or. (restart == "n")) .and. &
        & ((restart_filename /= "no") .or. (dim_filename /= "no"))) then
        print*, "Restart option not activated, restart_filename or dim_filename will not be used"
    end if

    if (init == "christensen") then
        print*, "Starting simulation from Christensen et al. (2001) initial condition"
    else if (init == "symmetric") then
        print*, "Starting simulation from symmetric initial condition"
    else if (init /= 'no') then
        print*, "Simulation stopped - init has an invalid value"
        stop
    end if

    if (init == "symmetric") then
        if (init_amp <= 0.0) then
            print *, "Simulation stopped - init_amp has an invalid value"
            stop
        else
            print*, "The amplitude for the symmetric initial condition is ", init_amp
        end if
    end if

    if (init == "symmetric") then
        if (sym <= 0) then
            print*, "Simulation stopped - sym has an invalid value"
            stop
        else
            print*, "The azimuthal symmetry of the initial condition is ", sym
        end if
    end if

    if (delta_t <= 0.0) then
        print *, "Simulation stopped - delta_t has an invalid value"
        stop
    else
        print*, "delta_t = ", delta_t
    end if

    if (((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) &
        & .and. (NTS /= 0)) then
        print*, "Newton solver selected, NTS is of no use"
    else if (((solver == "convective_explicit") .or. (solver == "convective_implicit")) &
        & .and. (NTS <= 0)) then
        print *, "Simulation stopped - NTS has an invalid value"
        stop
    else if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
        print*, "Number of iterations = ", NTS
    end if

    if (((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) &
        & .and. (save_every /= 0)) then
        print*, "Newton solver selected, save_every is of no use"
    else if (((solver == "convective_explicit") .or. (solver == "convective_implicit")) &
        & .and. (save_every <= 0)) then
        print *, "Simulation stopped - save_every has an invalid value"
        stop
    else if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
        print*, "Flow files will be saved every ", save_every, " iterations"
    end if

    if (((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) &
        & .and. (save_restart /= 0)) then
        print*, "Newton solver selected, save_restart is of no use"
    else if (((solver == "convective_explicit") .or. (solver == "convective_implicit")) &
        & .and. (save_restart <= 0)) then
        print *, "Simulation stopped - save_restart has an invalid value"
        stop
    else if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
        print*, "Restart files will be saved every ", save_restart, " iterations"
    end if

    if (((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") &
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) &
        & .and. (save_Ur_mgep_t /= 0)) then
        print*, "Newton solver selected, save_Ur_mgep_t is of no use"
    else if (((solver == "convective_explicit") .or. (solver == "convective_implicit")) &
        & .and. (save_Ur_mgep_t < 0)) then
        print *, "Simulation stopped - save_Ur_mgep_t has an invalid value"
        stop
    else if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
        print*, "Time series of Ur in equatorial plane-mid gap will be saved in the last ", save_Ur_mgep_t, " iterations"
    end if

    if (KK <= 0) then
        print *, "Simulation stopped - KK has an invalid value"
        stop
    else
        print*, "The number of Chebyshev modes is ", KK
    end if

    if (LL <= 0) then
        print *, "Simulation stopped - LL has an invalid value"
        stop
    else
        print*, "The number of Legendre modes is ", LL
    end if

    if (MM <= 0) then
        print *, "Simulation stopped - MM has an invalid value"
        stop
    else
        print*, "The number of Fourier modes is ", MM
    end if

    if (mres < 0) then
        print *, "Simulation stopped - mres has an invalid value"
        stop
    else if (mres == 0) then
        print*, "mres not fixed, setting to mres = 1"
        mres = 1
    else
        print*, "The resolution of the Fourier transform is ", mres
    end if

    if (((solver == "convective_explicit") .or. (solver == "convective_implicit")) .and. &
        & ((max_newt /= 0) .or. (max_gmres /= 0) .or. (restart_gmres /= 0) .or. (newt_eps /= 0) & 
        & .or. (newt_delta /= 0) .or. (tol_gmres /= 0) .or. (M_wave /= 0))) then
        print*, "Newton solver not selected, values of max_newt, max_gmres, restart_gmres, newt_eps, " // &
                & "newt_delta, tol_gmres or M_wave will not be used"
    else if ((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") & 
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) then
        if (max_newt <= 0) then
            print*, "Simulation stopped - max_newt has an invalid value"
            stop
        else
            print*, "The number of maximum Newton iterations for the Newton solver is ", max_newt
        end if

        if (newt_eps <= 0) then
            print*, "Simulation stopped - newt_eps has an invalid value"
            stop
        else
            print*, "The Newton epsilon criteria for the Newton solver is ", newt_eps
        end if

        if (newt_delta <= 0) then
            print*, "Simulation stopped - newt_delta has an invalid value"
            stop
        else
            print*, "The Newton delta criteria for the Newton solver is ", newt_delta
        end if

        if (max_gmres <= 0) then
            print*, "Simulation stopped - max_gmres has an invalid value"
            stop
        else
            print*, "The number of maximum GMRESm iterations for the Newton solver is ", max_gmres
        end if

        if (restart_gmres <= 0) then
            print*, "Simulation stopped - restart_gmres has an invalid value"
            stop
        else
            print*, "GMRESm will restart every ", restart_gmres, " iterations"
        end if

        if (tol_gmres <= 0) then
            print*, "Simulation stopped - tol_gmres has an invalid value"
            stop
        else
            print*, "The GMRESm tolarance for the Newton solver is ", tol_gmres
        end if

        if (m_wave <= 0) then
            print*, "Simulation stopped - M_wave has an invalid value"
            stop
        else
            print*, "The solution the Newton solver will look for has an m = ", M_wave, ' azimuthal symmetry'
        end if

    end if

    if ((solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) then

        if ((Ek_final /= 0.) .and. (Ra_final /= 0.)) then
            print*, "Continuation solver selected but final values for both Ekman and Rayleigh were given, choose one."
            print*, "Stopping simulation..."
            stop
        end if

        if (Ek_final < 0.) then
            print*, "Simulation stopped - Ek_final has an invalid value"
            stop
        else if (Ek_final > 0.) then
            print*, "Continuation in Ekman until Ek_final = ", Ek_final
        end if

        if (Ra_final < 0.) then
            print*, "Simulation stopped - Ra_final has an invalid value"
            stop
        else if (Ra_final > 0.) then
            print*, "Continuation in Rayleigh until Ra_final = ", Ra_final
        end if

        if ((adapt_param_string == "yes") .or. (adapt_param_string == "y")) then
            print*, "Continuation will adapt the step of the parameter"
            adapt_param = .true.
        else
            print*, "Continuation will not adapt the step of the parameter (default setting)"
            adapt_param = .false.
        end if

        if (Nopt <= 0.) then
            print*, "Simulation stopped - Nopt has an invalid value"
            stop
        else
            print*, "Continuation using Nopt = ", Nopt
        end if

        if (delta_param == 0.) then
            print*, "Simulation stopped - delta_param has an invalid value"
            stop
        else
            print*, "Continuation using delta_param = ", delta_param
        end if

        if (gamma <= 0.) then
            print*, "Simulation stopped - gamma has an invalid value"
            stop
        else
            print*, "Continuation using gamma = ", gamma
        end if

        if ((grid_refine_string == "yes") .or. (grid_refine_string == "y")) then
            print*, "Continuation will be performed with grid refinement"
            grid_refine = .true.
            if (gr_threshold <= 0.) then
                print*, "Simulation stopped - gr_threshold has an invalid value"
                stop
            else
                print*, "Grid refinement using gr_threshold = ", gr_threshold
            end if
            if (gr_threshold > 1.0e-4) print*, "Carefull, grid refinement threshold might be too high"
        else
            print*, "Continuation will be performed without grid refinement (default setting)"
            grid_refine = .false.
        end if

    end if

    if ((dealiasing /= "yes") .or. (dealiasing /= "y")) then
        print*, "Using SHTns dealiasing"
    else if ((dealiasing /= "no") .or. (dealiasing /= "n")) then
        print*, "Dealiasing desabled"
    else
        print *, "Simulation stopped - dealiasing has an invalid value"
        stop
    end if

    if (Rin <= 0) then
        print *, "Simulation stopped - Rin has an invalid value"
        stop
    else
        print*, "The inner radius is ", Rin
    end if

    if (Rout <= 0) then
        print *, "Simulation stopped - Rout has an invalid value"
        stop
    else
        print*, "The outer radius is ", Rout
    end if

    
    if (Pr <= 0) then
        print *, "Simulation stopped - Pr has an invalid value"
        stop
    else
        print*, "The Prandtl number is ", Pr
    end if

    
    if (Ek <= 0) then
        print *, "Simulation stopped - Ek has an invalid value"
        stop
    else
        print*, "The Ekman number is ", Ek
    end if

    if (((solver == "convective_explicit") .or. (solver == "convective_implicit").or. &
        & (solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit") &
        & .or. (solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) .and. &
        & (Ro /= 0.)) then
        print*, "Convective solver selected, the Rossby number will not be used"
    else if ((solver == "mhd_explicit") .or. (solver == "mhd_implicit")) then
        if (Ro <= 0) then
            print *, "Simulation stopped - Ro has an invalid value"
            stop
        else
            print*, "The Rossby number is ", Ro
        end if
    end if
    
    if (Ra <= 0) then
        print *, "Simulation stopped - Ra has an invalid value"
        stop
    else
        print*, "The Rayleigh number is ", Ra
    end if

    if (time_step == "fbe") then
        IER = 1.
    else if (time_step == "bdf2") then
        IER = 2. / 3.
    else
        if (IER <= 0) then
            print *, "Simulation stopped - IER has an invalid value"
            stop
        else
            print*, "The implicit to explicit ratio is ", IER
        end if
    end if

end subroutine check_arg_validity

end module mod_CommandLineParser