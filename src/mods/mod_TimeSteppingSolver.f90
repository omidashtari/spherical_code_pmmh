! spherical_code_pmmh: solve the Boussinesq equations in a rotating spherical shell
! Copyright (C) 2025 Alan Riquier, Camille Rambert, Juan Cruz Gonzalez Sembla and Laurette Tuckerman

! Contact: juan-cruz.gonzalez-sembla@espci.fr, laurette.tuckerman@espci.fr

! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

! This module was written by J. C. Gonzalez Sembla.
! It contains the convective timestepping solver.

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
    use, intrinsic :: ieee_arithmetic
  
    implicit none
  
contains

subroutine convective_solver()

    implicit none

    integer :: ios

    procedure(), pointer :: TimeStep_ptr

    ! Check for solver type
    if (solver == "convective_explicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        call precompBuildXY() ! Compute building blocks for the matrices
        
        ! Check for time stepping scheme
        select case (time_step)

            case ("pc")

                print*, "Precomputing the Xe and Ye matrices..."
                call precompXeYe()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
                TimeStep_ptr => compute_time_step_convective_explicit_PC

            case ("cn")

                print*, "Precomputing the Xe and Ye matrices..."
                call precompXeYe()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
                TimeStep_ptr => compute_time_step_convective_explicit_CN_FBE
                Explicit_RHS_ptr => comp_ExplicitRHS

            case ("fbe")

                print*, "Precomputing the Xe and Ye matrices..."
                call precompXeYe_BDF()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                TimeStep_ptr => compute_time_step_convective_explicit_CN_FBE
                Explicit_RHS_ptr => comp_ExplicitRHS

            case ("bdf2")

                print*, "Precomputing the Xe and Ye matrices..."
                call precompXeYe_BDF()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                TimeStep_ptr => compute_time_step_convective_explicit_BDF2

        end select

    else if (solver == "convective_implicit") then

        ! We begin by precomputing the matrices for the resolution of the linear system
        call precompBuildXY() ! Compute building blocks for the matrices
        
        ! Check for time stepping scheme
        select case (time_step)

            case ("pc")

                print*, "Precomputing the Xef and Yef matrices..."
                call PrecompimplicitXY()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
                TimeStep_ptr => compute_time_step_convective_implicit_PC

            case ("cn")

                print*, "Precomputing the Xef and Yef matrices..."
                call PrecompimplicitXY()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
                TimeStep_ptr => compute_time_step_convective_implicit_CN
                Implicit_RHS_ptr => comp_ImplicitRHS

            case ("fbe")

                print*, "Precomputing the Xef, Ye and Yf matrices..."
                call PrecompimplicitXY_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                TimeStep_ptr => compute_time_step_convective_implicit_FBE
                Implicit_RHS_ptr => comp_ImplicitRHS

            case ("bdf2")

                print*, "Precomputing the Xef, Ye and Yf matrices..."
                call PrecompimplicitXY_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                TimeStep_ptr => compute_time_step_convective_implicit_BDF2

        end select

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

    ! Open file for KE time series, checking if file already exists
    open(51, file=trim(directory)//"/KE_timeserie.dat", status='unknown', action='read', iostat=ios)
    close(51)
    if (ios /= 0) then 
        ! File does not exist
        open(51,file=trim(directory)//"/KE_timeserie.dat", status='unknown', position='append')
        write(51, "(A15, 4x, A23)") "time", "E_kin"
    else
        ! File exists
        open(51,file=trim(directory)//"/KE_timeserie.dat", status='unknown', position='append')
    end if

    ! Open file to save timeseries of Ur in mid gap, equatorial plane. Checking if file already exists
    open(52, file=trim(directory)//"/Ur_mgep_timeserie.dat", status='unknown', action='read', iostat=ios)
    close(52)
    if (ios /= 0) then 
        ! File does not exist
        open(52,file=trim(directory)//"/Ur_mgep_timeserie.dat", status='unknown', position='append')
        write(52, "(A15, 4x, A21)") "time", "Ur"
    else
        ! File exists
        open(52,file=trim(directory)//"/Ur_mgep_timeserie.dat", status='unknown', position='append')
    end if

    print*, "Volume of the shell = ", Vol

    print*
    print*, "################## Starting the simulation #####################"
    print*

    ! If we are using BDF2 but we were not able to read data on t-1 from restart file
    if ((read_tm1 .eqv. .false.) .and. (time_step == "bdf2")) then
        ! Initialise states so that the first timestep works as an FBE with dt = 2/3*delta_t
        ! Save the states
        E_tm1 = E ; F_tm1 = F ; T_tm1 = T ;
        ! Now we compute and save the RHS
        if (solver == "convective_explicit") then
            call comp_ExplicitRHS(DE, DF, DT)
        else
            call comp_ImplicitRHS(DE, DF, DT)
        end if
        DE_tm1 = DE ; DF_tm1 = DF ; DT_tm1 = DT ;
    end if

    do step = step_0 + 1, NTS

        print*, "############## Starting time step n°", step
        print*

        time = time + delta_t

        print*, "Time = ", time
        print*

        ! Now we compute the explicit timestep
        call TimeStep_ptr()

        ! Compute, print and save kinetic energy
        call comp_KineticEnergy(Ur, Ut, Up)
        ! Check if Ekin is NaN
        if (ieee_is_nan(Ekin)) then
            print*, "The Kinetic Energy has turned into a NaN"
            print*, "Stopping simulation..."
            stop
        end if
        print*, "Kinetic energy:", Ekin
        write(51,"(E24.16,3x,E24.16)") time, Ekin

        ! Save time series of Ur in mid gap, equatorial plane
        write(52,"(E24.16,3x,E24.16)") time, Ur(kN / 2, lN / 2, 1)

        ! Save restart files and output files
        if (mod(step, save_every)==0) then
            call Output_files(Ur, Up, T_real, step=step)
            call writeRestart(step=step)
            call writeDim(step=step)
        end if

        print*
        print*, "############## End of time step n°  ",step
        print*

    end do

end subroutine convective_solver

end module mod_TimeSteppingSolver