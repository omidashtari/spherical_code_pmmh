module mod_tools

    use iso_c_binding

    use mod_GlobalVars
    use mod_PrecompSH
    use mod_Tchebyshev
    use mod_PrecompXY
    use mod_Output
    use mod_read

contains

subroutine grid_refinement()

    implicit none

    integer :: delta_SH ! Increase in modes for SH
    character(len = 20) :: Ek_str, suffix ! To store Ek in scientific notation
    integer :: i

    refinement_count = refinement_count + 1

    ! First we deallocate everything that was previously allocated.
    print*, "Deallocating old fields"

    if (ur_Cheb) then
        deallocate(xN, rN, xK, rK)
        deallocate(Chb, ChbR, ChbR2, ChbR3, ChbD1, ChbD1R, Chbderiv, Chbderiv2, Chbinv, Chb2, Chb4, Chb_mulR_deriv)
        deallocate(Chb0XY, ChbD1XY, ChbD2XY, ChbD4XY, divR)
        deallocate(BC_mat)
    end if

    deallocate(ell, SH_norm, phi, ll1, dphi, wts, costh_mul, sth_dth)
    deallocate(E, F, T)
    deallocate(DE, DF, DT)
    deallocate(DEp, DFp, DTp)
    deallocate(wT, wDT)
    deallocate(E_base, F_base, T_base)
    deallocate(Ur, Ut, Up, T_real)
    deallocate(E_nm1, F_nm1, T_nm1, E_nm2, F_nm2, T_nm2, E_nm3, F_nm3, T_nm3)

    if (solver == "continuation_convective_explicit") then
        deallocate(wE, wDE, wF, wDF)
        deallocate(Xe_inv, Xe_invYe, Xf_inv, Xf_invYf, XT_inv, XT_invYT) 
    end if

    if (solver == "continuation_convective_implicit") then
        deallocate(A, Ap, A1, DA, DA1)
        deallocate(DEr, DFr, DEpr, DFpr)
        deallocate(Xef, Yef, PIVOT, XT_inv, XT_invYT)
        deallocate(idx_even, idx_odd)
    end if

    ! Destroy previous SHTns grid
    call shtns_unset_grid(shtns_c)
    call shtns_destroy(shtns_c)

    print*, "Deallocation completed"
    print*, "Computing new grid"
    
    ! Defining new resolution
    if (ur_Cheb) then
        KK = KK + ceiling(0.2 * KK) ! Increase by 20%
    end if

    if (ur_SH) then
        delta_SH = ceiling((0.2 * LL) / mres)
        LL = LL + delta_SH * mres
        MM = MM + delta_SH
    end if

    ! Output new grid parameters in file
    open(30,file=trim(directory)//"/parameters.dat", status='unknown', position='append')
    write(30, '(a,3x,E24.16)') "Underresolution detected when doing continuation in Ekman at Ek = ", Ek
    write(30, '(a)') "The new grid parameters are:"
    write(30, '(a,3x,I5,a,3x,I5,a,3x,I5)') "KK = ", KK, " LL = ", LL, " MM = ", MM
    write(30, '(a,3x,I5)') "r, theta and phi files will henceforth be identified by the following number: ", refinement_count
    write(30, '(a)')
    close(30)

    if (ur_Cheb) then
        print*, "Precomputing the Chebyshev terms..."
        call PrecompTchebyshev()
        call precompBuildXY() ! Compute building blocks for the matrices
    end if
    print*, "Precomputing the spherical harmonics grid..."
    call PrecompSH()
    call output_coordinates(step=refinement_count)

    ! Now we re-compute the matrices
    if (solver == "continuation_convective_explicit") then

        if (time_step == "cn") then

            print*, "Precomputing the Xe and Ye matrices..."
            call precompXeYe()
            print*, "Precomputing the Xf and Yf matrices..."
            call precompXfYf()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT()

        else if (time_step == "iee") then

            print*, "Precomputing the Xe and Ye matrices..."
            call precompXeYe_BDF()
            print*, "Precomputing the Xf and Yf matrices..."
            call precompXfYf_BDF()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT_BDF()

        end if

    else if (solver == "continuation_convective_implicit") then

        if (time_step == "cn") then

            print*, "Precomputing the Xef and Yef matrices..."
            call PrecompimplicitXY()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT()

        else if (time_step == "iee") then

            print*, "Precomputing the Xef and Yef matrices..."
            call PrecompimplicitXY_BDF()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT_BDF()

        end if

    end if

    ! Now we restart from previous continuation step
    print*, "Reading restart files from previous continuation step using new grid"

    ! First we change the names of the restart files
    write(Ek_str, "(ES10.1E2)") Ek  ! Scientific notation with format ES (e.g., 1e-3)
    do i = 1, len(Ek_str)
        if (Ek_str(i:i) == 'E') then
        Ek_str(i:i) = 'e'
        end if
    end do
    suffix = "_Ek_" // trim(adjustl(Ek_str))
    restart_filename = "Restart"//trim(adjustl(suffix))//".b"
    dim_filename = "Dim"//trim(adjustl(suffix))//".b"
    
    ! And then we call the restarting subroutines
    call readDim()
    call readRestart()

    ! Setting pointers
    ! Set pointers for E, F and T
    call c_f_pointer(c_loc(E), E_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F), F_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T), T_ptr, [2 * KK2 * shtns%nlm])

    ! Set pointers for E_base, F_base and T_base
    call c_f_pointer(c_loc(E_base), E_base_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F_base), F_base_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T_base), T_base_ptr, [2 * KK2 * shtns%nlm])

    print*, "Grid refinement completed"

end subroutine grid_refinement

end module mod_tools