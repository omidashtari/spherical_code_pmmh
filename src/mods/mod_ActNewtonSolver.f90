module mod_ActNewtonSolver

    use iso_c_binding

    use mod_GlobalVars
    use mod_Matrices
    use mod_PrecompSH
    use mod_ExplicitTerms
    use mod_TimeStep

    implicit none

contains

subroutine act_on_RHS(time_stepping_sub, E_RHS, F_RHS, T_RHS)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2 * KK2 * shtns%nlm), intent(out) :: E_RHS
    double precision, dimension(2 * KK4 * shtns%nlm), intent(out) :: F_RHS
    double precision, dimension(2 * KK2 * shtns%nlm), intent(out) :: T_RHS

    E_base = E
    F_base = F
    T_base = T

    call time_stepping_sub()

    ! Compute F(U)
    E_RHS = E_ptr - E_base_ptr
    F_RHS = F_ptr - F_base_ptr
    T_RHS = T_ptr - T_base_ptr

end subroutine act_on_RHS

subroutine subA(time_stepping_sub, EFT_in, EFT_out)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2 * shtns%nlm * (2 * KK2 + KK4)), intent(inout) :: EFT_in
    double precision, dimension(2 * shtns%nlm * (2 * KK2 + KK4)), intent(out) :: EFT_out

    ! Extract c_per and set phase
    c_per = EFT_in(wavespeed_loc) ! We are storing the wavespeed c in the place where we inforce the phase
    EFT_in(wavespeed_loc) = 0.    ! Fixing the phase

    ! Set EFT_out = EFT_in
    EFT_out = EFT_in

    call act_lin_non_lin(time_stepping_sub, EFT_in(1), EFT_in(2 * shtns%nlm * KK2 + 1), &
                        & EFT_in(2 * shtns%nlm * (KK2 + KK4) + 1))

    ! Compute difference between input and output
    EFT_out = EFT_in - EFT_out

end subroutine subA

subroutine act_lin_non_lin(time_stepping_sub, E_per_ptr, F_per_ptr, T_per_ptr)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2 * KK2 * shtns%nlm), intent(inout), target :: E_per_ptr
    double precision, dimension(2 * KK4 * shtns%nlm), intent(inout), target :: F_per_ptr
    double precision, dimension(2 * KK2 * shtns%nlm), intent(inout), target :: T_per_ptr

    ! Set pointers for E_per, F_per and T_per
    call c_f_pointer(c_loc(E_per_ptr), E_per, [KK2, shtns%nlm])
    call c_f_pointer(c_loc(F_per_ptr), F_per, [KK4, shtns%nlm])
    call c_f_pointer(c_loc(T_per_ptr), T_per, [KK2, shtns%nlm])

    ! Call linearized time step
    call time_stepping_sub(E_per, F_per, T_per)

end subroutine act_lin_non_lin

subroutine update_EFT(E_best, F_best, T_best, norm_EFT_new)

    implicit none

    double precision, dimension(2 * KK2 * shtns%nlm), intent(in) :: E_best
    double precision, dimension(2 * KK4 * shtns%nlm), intent(in) :: F_best
    double precision, dimension(2 * KK2 * shtns%nlm), intent(in) :: T_best

    double precision, intent(out) :: norm_EFT_new

    ! Update
    E_ptr = E_base_ptr - E_best
    F_ptr = F_base_ptr - F_best
    T_ptr = T_base_ptr - T_best

    ! Compute the norm
    norm_EFT_new = sqrt(dot_product(E_ptr, E_ptr) + dot_product(F_ptr, F_ptr) + &
                    & dot_product(T_ptr, T_ptr) + C_base ** 2)

end subroutine update_EFT

end module mod_ActNewtonSolver