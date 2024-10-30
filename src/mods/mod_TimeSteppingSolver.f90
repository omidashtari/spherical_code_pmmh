module mod_TimeSteppingSolver

    use mod_Globalvars
    use mod_Matrices
    use mod_PrecompSH
    use mod_Tchebyshev
    use mod_precompXY
    use mod_read
    use mod_init
    use mod_TimeStep
    use mod_IterativeSolvers
    use mod_Output
  
    implicit none
  
contains

subroutine convective_solver()

    implicit none

    procedure(), pointer :: TimeStep_ptr

    ! Check for solver type
    if (solver == "convective_explicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xe and Ye matrices..."
        call precompXeYe()
        print*, "Precomputing the Xf and Yf matrices..."
        call precompXfYf()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        ! Check for time stepping scheme
        if (time_step == "pc") then
            TimeStep_ptr => compute_time_step_convective_explicit_PC
        else if (time_step == "iee") then
            TimeStep_ptr => compute_time_step_convective_explicit_IEE
            Explicit_RHS_ptr => comp_ExplicitRHS
        end if

    else if (solver == "convective_implicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        print*, "Precomputing the Xef and Yef matrices..."
        call PrecompimplicitXY()
        print*, "Precomputing the XT and YT matrices..."
        call precompXTYT()

        ! Check for time stepping scheme
        if (time_step == "pc") then
            TimeStep_ptr => compute_time_step_convective_implicit_PC
        else if (time_step == "iee") then
            TimeStep_ptr => compute_time_step_convective_implicit_IEE
            Implicit_RHS_ptr => comp_ImplicitRHS
        end if

    end if

    print*
    print*, "--------------------- Initial values ---------------------------"
    print*

    if ((restart == 'y') .or. (restart == 'yes')) then
    print*, "Reading restart files"
    call readDim()
    call readRestart()
    print*, "Starting from t = ", t0
    print*, "Starting from iteration = ", step_0
    else
    print*, "Starting simulation from initial condition"
    call init_EF()
    call init_TEMP()
    end if

    time = t0
    NTS = NTS + step_0

    ! Open file for KE time series
    open(51,file=trim(directory)//"/KE_timeserie.dat", status='unknown', position='append')
    write(14, "(A12, 4x, A16)") "time", "E_kin"

    ! Open file for saven Ur in mid gap, equatorial plane
    if (save_Ur_mgep_t /= 0) open(52,file=trim(directory)//"/Ur_mgep_timeserie.dat", status='unknown', position='append')
    write(14, "(A12, 4x, A16)") "time", "Ur"


    print*, "Volume of the shell = ", Vol

    print*
    print*, "################## Starting the simulation #####################"
    print*

    do step = step_0 + 1,NTS

    print*, "############## Starting time step n°", step
    print*

    time = time + delta_t

    print*, "Time = ",time
    print*

    ! Now we compute the explicit timestep
    call TimeStep_ptr()

    ! Save restart files
    if ((mod(step, save_restart) == 0) .or. step == NTS) then
        call writeRestart()
        call writeDim()
    end if

    ! Compute, print and save kinetic energy
    call comp_KineticEnergy(Ur, Ut, Up)
    print*, "Kinetic energy:", Ekin
    write(51,"(E16.6,3x,E16.6)") time, Ekin

    ! Save time series of Ur in mid gap, equatorial plane
    if ((save_Ur_mgep_t > 0) .and. (step >= (NTS-save_Ur_mgep_t))) then
        write(52,"(E16.6,3x,E16.6)") time, Ur(kN / 2, lN / 2, 1)
    end if


    if (mod(step, save_every)==0) then
        call Output_files(Ur, Up, T_real, step=step)
    end if

    print*
    print*, "############## End of time step n°  ",step
    print*

    end do

end subroutine convective_solver

end module mod_TimeSteppingSolver