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

! This module was written by A. Riquier, C. Rambert and J. C. Gonzalez Sembla.
! It contains the subroutines to compute the RHS and transform into physical 
! space and back. Subroutines:
! - comp_ExplicitRHS: compute RHS for timestep using explicit Coriolis.
! - comp_ImplicitRHS: compute RHS for timestep using implicit Coriolis.
! - comp_RHS_with_rot: compute RHS taking into account rotation from drifting 
! frequency for Newton solver.
! - comp_LinNonLin: compute linearized non-linear term.
! - comp_cdphi: compute d/dphi and multiply by frequency c.
! - ToReal: go to real space.
! - BackToSpectral: go to spectral space.
! - comp_radial_derivatives: compute radial derivatives from E, F and T and saves in auxiliary fields.
! - comp_U_cU_gT: compute U, curl(U) and grad(T) based on the auxiliary fields.
! - comp_U: compute U based on the auxiliary fields.
! - comp_curlU: compute curl(U) based on the auxiliary fields.
! - comp_gradT: compute grad(T) based on the auxiliary fields.
! - comp_U_from_EF: compute U from poloidal and toroidal potentials.
! - rCurl: compute rCurlF and rCurlCurlF.


module mod_ExplicitTerms

  use mod_Globalvars
  use mod_Matrices
  use mod_PrecompSH
  use mod_Tchebyshev

  implicit none

contains

subroutine comp_ExplicitRHS(DE, DF, DT)

  !###########################################################################
  !   Computation of RHS for the convective solver with explicit Coriolis
  !###########################################################################

  implicit none

  double complex, dimension(KK2, shtns%nlm), intent(out) :: DE
  double complex, dimension(KK4, shtns%nlm), intent(out) :: DF
  double complex, dimension(KK2, shtns%nlm), intent(out) :: DT

  double precision, dimension(kN, lN, mN) :: cUr, cUt, cUp ! Components of the curl of U in real space

  double precision, dimension(kN, lN, mN) :: gTr, gTt, gTp ! Components of the gradient of T in real space

  double precision, dimension(kN, lN, mN) :: Fr, Gt, Gp    ! Forcing terms in real space
  double precision, dimension(kN, lN, mN) :: UgradT        ! Temp equaqtion

  double complex, dimension(KK, shtns%nlm) :: Fr_spec, Gt_spec, Gp_spec ! Forcing terms in spectral space

  double complex, dimension(KK, shtns%nlm) :: rCF_spec, rCCF_spec ! rcurl(F) and rcurl(curl(F)) in spectral space

  double complex, dimension(KK, shtns%nlm) :: UgradT_spec ! UgradT in spectral space

  integer :: k, l

  DE = 0.
  DF = 0.
  DT = 0.

  call comp_radial_derivatives(E, F, T, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
  call comp_U_cU_gT(Ur, Ut, Up, cUr, cUt, cUp, gTr, gTt, gTp)
  call ToReal(T, T_real, KK2)

  !--- Forcing terms for U

  ! Note 1: we write (U.grad)U = curl(U) x U + 1/2 * grad(U²)
  !       and there is no need to compute grad(U²) since we will take
  !       the curl of F later
  ! Note 2: As in comp_U, comp_curlU and comp_gradT the theta and phi
  !       components are multiplied by sin(theta) we need to take care of
  !       this when we compute the forcing terms

  do k = 1, kN
    do l = 1, lN
      Fr(k, l, :) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) &
                  & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) &
                  & + 2.d0 * Up(k,l,:) * SinTh(l)
      Gt(k, l, :) = (Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) &
                  & + 2.d0 * CosTh(l) * Up(k,l,:)) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
      Gp(k, l, :) = (Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) &
                  & - 2.d0 * (CosTh(l)*Ut(k,l,:) + SinTh(l)*Ur(k,l,:))) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
      UgradT(k, l, :) = Ur(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                      & - Ut(k,l,:) * gTt(k,l,:) - Up(k,l,:) * gTp(k,l,:)
    end do
  end do

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurl(Fr_spec, Gt_spec, Gp_spec, rCF_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  ! When solving for total temperature:
  ! DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
                             ! We fix the l=m=0 mode (cosine part) by fixing it to 1. / norm(0, 0), which is equal to 1.

  ! For E
  DE(KK + 1, :) = 0.
  DE(KK + 2, :) = 0.

  ! For F
  DF(KK + 1, :) = 0.
  DF(KK + 2, :) = 0.
  DF(KK + 3, :) = 0.
  DF(KK + 4, :) = 0.

end subroutine comp_ExplicitRHS

!----------------------------------------------------------------------------

subroutine comp_ImplicitRHS(DE, DF, DT)

  !###########################################################################
  !   Computation of RHS for the convective solver with implicit Coriolis
  !###########################################################################

  implicit none

  double complex, dimension(KK2, shtns%nlm), intent(out) :: DE
  double complex, dimension(KK4, shtns%nlm), intent(out) :: DF
  double complex, dimension(KK2, shtns%nlm), intent(out) :: DT

  double precision, dimension(kN, lN, mN) :: cUr, cUt, cUp ! Components of the curl of U in real space

  double precision, dimension(kN, lN, mN) :: gTr, gTt, gTp ! Components of the gradient of T in real space

  double precision, dimension(kN, lN, mN) :: Fr, Gt, Gp    ! Forcing terms in real space
  double precision, dimension(kN, lN, mN) :: UgradT        ! Temp equaqtion

  double complex, dimension(KK, shtns%nlm) :: Fr_spec, Gt_spec, Gp_spec ! Forcing terms in spectral space

  double complex, dimension(KK, shtns%nlm) :: rCF_spec, rCCF_spec ! rcurl(F) and rcurl(curl(F)) in spectral space

  double complex, dimension(KK, shtns%nlm) :: UgradT_spec ! UgradT in spectral space

  integer :: k, l

  DE = 0.
  DF = 0.
  DT = 0.

  call comp_radial_derivatives(E, F, T, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
  call comp_U_cU_gT(Ur, Ut, Up, cUr, cUt, cUp, gTr, gTt, gTp)
  call ToReal(T, T_real, KK2)
  
  do k = 1, kN
    do l = 1, lN
      Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) &
                  & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1))
      Gt(k,l,:) = Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
      Gp(k,l,:) = Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
      UgradT(k,l,:) = Ur(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                    & - Ut(k,l,:) * gTt(k,l,:) - Up(k,l,:) * gTp(k,l,:)
    end do
  end do

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurl(Fr_spec, Gt_spec, Gp_spec, rCF_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  ! When solving for total temperature:
  ! DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
                             ! We fix the l=m=0 mode (cosine part) by fixing it to 1. / norm(0, 0), which is equal to 1.

  ! For E
  DE(KK + 1, :) = 0.
  DE(KK + 2, :) = 0.

  ! For F
  DF(KK + 1, :) = 0.
  DF(KK + 2, :) = 0.
  DF(KK + 3, :) = 0.
  DF(KK + 4, :) = 0.

end subroutine comp_ImplicitRHS

!----------------------------------------------------------------------------

subroutine comp_RHS_with_rot(DE, DF, DT)

  !###########################################################################
  !   Computation of RHS for the convective solver with explicit Coriolis
  !###########################################################################

  implicit none

  double complex, dimension(KK2, shtns%nlm), intent(out) :: DE
  double complex, dimension(KK4, shtns%nlm), intent(out) :: DF
  double complex, dimension(KK2, shtns%nlm), intent(out) :: DT

  double precision, dimension(kN, lN, mN) :: cUr, cUt, cUp ! Components of the curl of U in real space

  double precision, dimension(kN, lN, mN) :: gTr, gTt, gTp ! Components of the gradient of T in real space

  double precision, dimension(kN, lN, mN) :: Fr, Gt, Gp    ! Forcing terms in real space
  double precision, dimension(kN, lN, mN) :: UgradT        ! Temp equaqtion

  double complex, dimension(KK, shtns%nlm) :: Fr_spec, Gt_spec, Gp_spec ! Forcing terms in spectral space

  double complex, dimension(KK, shtns%nlm) :: rCF_spec, rCCF_spec ! rcurl(F) and rcurl(curl(F)) in spectral space

  double complex, dimension(KK, shtns%nlm) :: UgradT_spec ! UgradT in spectral space

  double complex, dimension(KK2, shtns%nlm) :: E_base_rot ! To store rotated E_base
  double complex, dimension(KK4, shtns%nlm) :: F_base_rot ! To store rotated F_base
  double complex, dimension(KK2, shtns%nlm) :: T_base_rot ! To store rotated T_base

  double precision, dimension(kN, lN, mN) :: Ur_rot, Ut_rot, Up_rot, T_real_rot ! To store rotated velocity components and 
                                                                                ! temperature of rotated base E, F and T

  integer :: k, l

  DE = 0.
  DF = 0.
  DT = 0.

  call comp_radial_derivatives(E, F, T, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
  call comp_U_cU_gT(Ur, Ut, Up, cUr, cUt, cUp, gTr, gTt, gTp)
  call ToReal(T, T_real, KK2)

  !--- Differentiating by phi and multiplying by C_base
  call comp_cdphi(C_base, E_base, F_base, T_base, E_base_rot, F_base_rot, T_base_rot)
  call comp_U_from_EF(E_base_rot, F_base_rot, Ur_rot, Ut_rot, Up_rot)
  call ToReal(T_base_rot, T_real_rot, KK2)

  !--- Forcing terms for U

  ! Note 1: we write (U.grad)U = curl(U) x U + 1/2 * grad(U²)
  !       and there is no need to compute grad(U²) since we will take
  !       the curl of F later
  ! Note 2: As in comp_U, comp_curlU and comp_gradT the theta and phi
  !       components are multiplied by sin(theta) we need to take care of
  !       this when we compute the forcing terms

  if ((solver == 'newton_convective_implicit') .or. (solver == 'continuation_convective_implicit')) then
    do k = 1,kN
      do l = 1,lN
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) &
                & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) + Ur_rot(k,l,:) * Ek
        Gt(k,l,:) = Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:) + Ut_rot(k,l,:)) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:) + Up_rot(k,l,:)) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = Ur(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                  & - Ut(k,l,:) * gTt(k,l,:) - Up(k,l,:) * gTp(k,l,:) + T_real_rot(k,l,:)
      end do
    end do
  else 
    do k = 1, kN
      do l = 1, lN
        Fr(k, l, :) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) &
                    & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) &
                    & + 2.d0 * Up(k,l,:) * SinTh(l) + Ur_rot(k,l,:) * Ek
        Gt(k, l, :) = (Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) &
                    & + 2.d0 * CosTh(l) * Up(k,l,:) + Ut_rot(k,l,:) * Ek) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
        Gp(k, l, :) = (Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) &
                    & - 2.d0 * (CosTh(l)*Ut(k,l,:) +  SinTh(l)*Ur(k,l,:)) + Up_rot(k,l,:) * Ek) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k, l, :) = Ur(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                        & - Ut(k,l,:) * gTt(k,l,:) - Up(k,l,:) * gTp(k,l,:) + T_real_rot(k,l,:)
      end do
    end do
  end if

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurl(Fr_spec, Gt_spec, Gp_spec, rCF_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  ! When solving for total temperature:
  ! DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
                             ! We fix the l=m=0 mode (cosine part) by fixing it to 1. / norm(0, 0), which is equal to 1.

  ! For E
  DE(KK + 1, :) = 0.
  DE(KK + 2, :) = 0.

  ! For F
  DF(KK + 1, :) = 0.
  DF(KK + 2, :) = 0.
  DF(KK + 3, :) = 0.
  DF(KK + 4, :) = 0.

end subroutine comp_RHS_with_rot

!----------------------------------------------------------------------------

subroutine comp_LinNonLin(E_base, F_base, T_base, E_per, F_per, T_per, DE, DF, DT)

  !###########################################################################
  !   Computation of linearized non linear term on perturbationRHS 
  !       for the steady states solver with explicit Coriolis
  !###########################################################################

  implicit none

  ! Input base fields and perturbations
  double complex, dimension(KK2, shtns%nlm), intent(in) :: E_base, E_per
  double complex, dimension(KK4, shtns%nlm), intent(in) :: F_base, F_per
  double complex, dimension(KK2, shtns%nlm), intent(in) :: T_base, T_per

  ! Output RHS for the timestep of the perturbation
  double complex, dimension(KK2, shtns%nlm), intent(out) :: DE
  double complex, dimension(KK4, shtns%nlm), intent(out) :: DF
  double complex, dimension(KK2, shtns%nlm), intent(out) :: DT

  double precision, dimension(kN, lN, mN) :: Ur_per, Ut_per, Up_per, T_real_per ! Velocity components and temperature for the perturbation

  double precision, dimension(kN, lN, mN) :: cUr, cUt, cUp             ! Components of the curl of U in real space
  double precision, dimension(kN, lN, mN) :: cUr_per, cUt_per, cUp_per ! Components of the curl of u in real space

  double precision, dimension(kN, lN, mN) :: gTr, gTt, gTp             ! Components of the gradient of T in real space
  double precision, dimension(kN, lN, mN) :: gTr_per, gTt_per, gTp_per ! Components of the gradient of t in real space

  double precision, dimension(kN, lN, mN) :: Fr, Gt, Gp    ! Forcing terms in real space
  double precision, dimension(kN, lN, mN) :: UgradT        ! Temp equaqtion

  double complex, dimension(KK, shtns%nlm) :: Fr_spec, Gt_spec, Gp_spec ! Forcing terms in spectral space

  double complex, dimension(KK, shtns%nlm) :: rCF_spec, rCCF_spec ! rcurl(F) and rcurl(curl(F)) in spectral space

  double complex, dimension(KK, shtns%nlm) :: UgradT_spec ! UgradT in spectral space

  double complex, dimension(KK2, shtns%nlm) :: E_base_rot, E_per_rot ! To store rotated E_base and E_per
  double complex, dimension(KK4, shtns%nlm) :: F_base_rot, F_per_rot ! To store rotated F_base and F_per
  double complex, dimension(KK2, shtns%nlm) :: T_base_rot, T_per_rot ! To store rotated T_base and T_per
  double precision, dimension(kN, lN, mN) :: Ur_rot, Ut_rot, Up_rot, T_real_rot ! To store rotated velocity components and temperature of rotated base state
  double precision, dimension(kN, lN, mN) :: Ur_per_rot, Ut_per_rot, Up_per_rot, T_real_per_rot  ! To store rotated velocity components and temperature of rotated perturbation state

  integer :: k, l

  DE = 0.
  DF = 0.
  DT = 0.

  call comp_radial_derivatives(E_base, F_base, T_base, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
  call comp_U_cU_gT(Ur, Ut, Up, cUr, cUt, cUp, gTr, gTt, gTp)

  call comp_radial_derivatives(E_per, F_per, T_per, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
  call comp_U_cU_gT(Ur_per, Ut_per, Up_per, cUr_per, cUt_per, cUp_per, gTr_per, gTt_per, gTp_per)

  call ToReal(T_base, T_real, KK2)                         ! Turn T to real
  call ToReal(T_per, T_real_per, KK2)                      ! Turn t to real

  !--- Differentiating by phi and multiplying by C_base
  call comp_cdphi(c_per, E_base, F_base, T_base, E_base_rot, F_base_rot, T_base_rot)
  call comp_cdphi(C_base, E_per, F_per, T_per, E_per_rot, F_per_rot, T_per_rot)

  call comp_U_from_EF(E_base_rot, F_base_rot, Ur_rot, Ut_rot, Up_rot)
  call comp_U_from_EF(E_per_rot, F_per_rot, Ur_per_rot, Ut_per_rot, Up_per_rot)

  call ToReal(T_base_rot, T_real_rot, KK2)
  call ToReal(T_per_rot, T_real_per_rot, KK2)
  !------

  !--- Forcing terms for U

  ! Note 1: we write (U.grad)U = curl(U) x U + 1/2 * grad(U²)
  !       and there is no need to compute grad(U²) since we will take
  !       the curl of F later
  ! Note 2: As in comp_U, comp_curlU and comp_gradT the theta and phi
  !       components are multiplied by sin(theta) we need to take care of
  !       this when we compute the forcing terms

  if ((solver == 'newton_convective_implicit') .or. (solver == 'continuation_convective_implicit')) then
    do k = 1,kN
      do l = 1,lN
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut_per(k,l,:) - CUt(k,l,:)*Up_per(k,l,:)) &
                    & + Ek * (CUp_per(k,l,:)*Ut(k,l,:) - CUt_per(k,l,:)*Up(k,l,:)) &
                    & + Ra * T_real_per(k,l,:) / (Rout*rN(k)**(-1))  + Ek * (Ur_rot(k,l,:) + Ur_per_rot(k,l,:))
        Gt(k,l,:) = (Ek * (CUr(k,l,:)*Up_per(k,l,:) - CUp(k,l,:)*Ur_per(k,l,:)) &
                    & + Ek * (CUr_per(k,l,:)*Up(k,l,:) - CUp_per(k,l,:)*Ur(k,l,:)) & 
                    & + Ek * (Ut_rot(k,l,:) + Ut_per_rot(k,l,:))) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = (Ek * (CUt(k,l,:)*Ur_per(k,l,:) - CUr(k,l,:)*Ut_per(k,l,:)) &
                    & + Ek * (CUt_per(k,l,:)*Ur(k,l,:) - CUr_per(k,l,:)*Ut(k,l,:)) & 
                    & + Ek * (Up_rot(k,l,:) + Up_per_rot(k,l,:))) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = - Ur(k,l,:) * gTr_per(k,l,:) &
                      & - Ut(k,l,:) * gTt_per(k,l,:) - Up(k,l,:) * gTp_per(k,l,:) &
                      & + Ur_per(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                      & - Ut_per(k,l,:) * gTt(k,l,:) - Up_per(k,l,:) * gTp(k,l,:) &
                      & + T_real_rot(k,l,:) + T_real_per_rot(k,l,:)
      end do
    end do
  else 
    do k = 1,kN
      do l = 1,lN
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut_per(k,l,:) - CUt(k,l,:)*Up_per(k,l,:)) &
                  & + Ek * (CUp_per(k,l,:)*Ut(k,l,:) - CUt_per(k,l,:)*Up(k,l,:)) &
                  & + Ra * T_real_per(k,l,:) / (Rout*rN(k)**(-1)) &
                  & + 2.d0 * Up_per(k,l,:) * SinTh(l) + Ek * (Ur_rot(k,l,:) + Ur_per_rot(k,l,:))
        Gt(k,l,:) = (Ek * (CUr(k,l,:)*Up_per(k,l,:) - CUp(k,l,:)*Ur_per(k,l,:)) &
                  & + Ek * (CUr_per(k,l,:)*Up(k,l,:) - CUp_per(k,l,:)*Ur(k,l,:)) &
                  & + 2.d0 * CosTh(l) * Up_per(k,l,:) & 
                  & + Ek * (Ut_rot(k,l,:) + Ut_per_rot(k,l,:))) / SinTh(l) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = (Ek * (CUt(k,l,:)*Ur_per(k,l,:) - CUr(k,l,:)*Ut_per(k,l,:)) &
                  & + Ek * (CUt_per(k,l,:)*Ur(k,l,:) - CUr_per(k,l,:)*Ut(k,l,:)) &
                  & - 2.d0 * (CosTh(l)*Ut_per(k,l,:) +  SinTh(l)*Ur_per(k,l,:)) & 
                  & + Ek * (Up_rot(k,l,:) + Up_per_rot(k,l,:))) / SinTh(l) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = - Ur(k,l,:) * gTr_per(k,l,:) &
                    & - Ut(k,l,:) * gTt_per(k,l,:) - Up(k,l,:) * gTp_per(k,l,:) &
                    & + Ur_per(k,l,:) * (Rout * Rin * rN(k)**(-2) - gTr(k,l,:)) &
                    & - Ut_per(k,l,:) * gTt(k,l,:) - Up_per(k,l,:) * gTp(k,l,:) &
                    & + T_real_rot(k,l,:) + T_real_per_rot(k,l,:)
      end do
    end do
  end if

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurl(Fr_spec, Gt_spec, Gp_spec, rCF_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.

  ! For E_per
  DE(KK + 1, :) = 0.
  DE(KK + 2, :) = 0.

  ! For F_per
  DF(KK + 1, :) = 0.
  DF(KK + 2, :) = 0.
  DF(KK + 3, :) = 0.
  DF(KK + 4, :) = 0.

end subroutine comp_LinNonLin

!----------------------------------------------------------------------------

subroutine comp_cdphi(z, U1, U2, U3, zdU1dphi, zdU2dphi, zdU3dphi)

  !###########################################################################
  !              Computation of z * \partial U / \partial phi
  !###########################################################################

  ! Note: z is the frequency of the rotating wave C (or c for the perturbation)
  ! U1 is E or E_per, U2 is F or F_per and U3 is T or Tt. 

  implicit none

  double precision, intent(in) :: z

  double complex, dimension(KK2, shtns%nlm), intent(in) :: U1 ! E or E_per
  double complex, dimension(KK4, shtns%nlm), intent(in) :: U2 ! F or F_per
  double complex, dimension(KK2, shtns%nlm), intent(in) :: U3 ! T or Tt

  double complex, dimension(KK2, shtns%nlm), intent(out) :: zdU1dphi ! z * dE/dphi or z * dE_per/dphi
  double complex, dimension(KK4, shtns%nlm), intent(out) :: zdU2dphi ! z * dF/dphi or z * dF_per/dphi
  double complex, dimension(KK2, shtns%nlm), intent(out) :: zdU3dphi ! z * dT/dphi or z * dTt/dphi

  integer :: lm

  ! Compute multiplication by z and differentiation by phi
  do lm = 1, shtns%nlm
    zdU1dphi(:, lm)%re = - z * aimag(U1(:, lm)) * dphi(lm)
    zdU1dphi(:, lm)%im = z * real(U1(:, lm)) * dphi(lm)

    zdU2dphi(:, lm)%re = - z * aimag(U2(:, lm)) * dphi(lm)
    zdU2dphi(:, lm)%im = z * real(U2(:, lm)) * dphi(lm)

    zdU3dphi(:, lm)%re = - z * aimag(U3(:, lm)) * dphi(lm)
    zdU3dphi(:, lm)%im = z * real(U3(:, lm)) * dphi(lm)
  end do

end subroutine comp_cdphi

!----------------------------------------------------------------------------

subroutine ToReal(T_spec, T_real, k_max)

  !###########################################################################
  !                Transform into real space from spectral
  !###########################################################################

  implicit none

  integer, intent(in) :: k_max ! Number of Chebyshev modes

  double complex,   dimension(k_max, shtns%nlm), intent(in) :: T_spec       ! Input spectral field
  double precision, dimension(kN, lN, mN),    intent(out) :: T_real      ! Output real field

  double precision, dimension(lN, mN) :: Sh             ! Intermediate array for SH transform
  double complex,   dimension(shtns%nlm) :: Slm         ! Intermediate array for SH transform
  double complex,   dimension(kN, shtns%nlm) :: T_inter ! Intermediate field

  integer :: k, lm

  Sh = 0.
  Slm = 0.
  T_inter = 0.
  T_real = 0.

  !#### 1. Chebyshev -> Real
  do lm = 1, shtns%nlm
    T_inter(:kN, lm)%re = real(T_spec(:k_max, lm)) .dot. Chb(:k_max, :kN)
    T_inter(:kN, lm)%im = aimag(T_spec(:k_max, lm)) .dot. Chb(:k_max, :kN)
  end do

  !#### 2. SH -> Real
  do k = 1, kN
    Slm = T_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    T_real(k, :, :) = Sh
  end do

end subroutine ToReal

!----------------------------------------------------------------------------

subroutine BackToSpectral(T_real, T_spec, k_max)

  !###########################################################################
  !                Transform into spectral space from real
  !###########################################################################

  implicit none

  integer, intent(in) :: k_max ! Number of Chebyshev modes

  double precision, dimension(kN, lN, mN),    intent(in) :: T_real  ! Input real field
  double complex,   dimension(k_max, shtns%nlm), intent(out) :: T_spec ! Output spectral field

  double precision, dimension(lN, mN) :: Sh             ! Intermediate array for SH transform
  double complex,   dimension(shtns%nlm) :: Slm         ! Intermediate array for SH transform
  double complex,   dimension(kN, shtns%nlm) :: T_inter ! Intermediate field

  integer :: k, lm

  Sh = 0.
  Slm = 0.
  T_inter = 0.
  T_spec = 0.

  ! #### 1. Real --> SH
  do k = 1, kN
    Sh = T_real(k, :, :)
    call spat_to_SH(shtns_c, Sh, Slm)
    T_inter(k, :) = Slm
  end do

  ! ### 2. Real --> Chebyshev
  do lm = 1, shtns%nlm
    T_spec(:k_max, lm)%re = real(T_inter(:kN, lm)) .dot. Chbinv(:kN, :k_max)
    T_spec(:k_max, lm)%im = aimag(T_inter(:kN, lm)) .dot. Chbinv(:kN, :k_max)
  end do

end subroutine BackToSpectral

!----------------------------------------------------------------------------

subroutine comp_radial_derivatives(E, F, T, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)

  !###########################################################################
  !         Compute Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St for future 
  !                 computation of U, curl(U) and grad(T)
  !###########################################################################

  implicit none

  double complex, dimension(KK2, shtns%nlm), intent(in) :: E, T ! Input e and T spectral fields
  double complex, dimension(KK4, shtns%nlm), intent(in) :: F    ! Input f spectral field

  double complex, dimension(shtns%nlm, kN), intent(out) :: Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St ! Resulting auxiliary arrays

  integer :: lm

  do lm = 1, shtns%nlm
    ! For U
    Qu(lm, :) = cmplx(real(F(:, lm)) .dot. ChbR2(:KK4, :kN), &
                      aimag(F(:, lm)) .dot. ChbR2(:KK4, :kN)) * ll1(lm)
    Su(lm, :) = cmplx(real(F(:, lm)) .dot. ChbD1R(:KK4, :kN), &
                      aimag(F(:, lm)) .dot. ChbD1R(:KK4, :kN))
    Tu(lm, :) = cmplx(real(E(:, lm)) .dot. ChbR(:KK2, :kN), &
                      aimag(E(:, lm)) .dot. ChbR(:KK2, :kN))

    ! For curl(U)
    Qcu(lm, :) = cmplx(real(E(:, lm)) .dot. ChbR2(:KK2, :kN), &
                      aimag(E(:, lm)) .dot. ChbR2(:KK2, :kN)) * ll1(lm)
    Scu(lm, :) = cmplx(real(E(:, lm)) .dot. ChbD1R(:KK2, :kN), &
                      aimag(E(:, lm)) .dot. ChbD1R(:KK2, :kN))
    Tcu(lm, :) = cmplx(real(F(:, lm)) .dot. (ChbR3(:KK4, :kN) * ll1(lm) - ChbD2R(:KK4, :kN)), &
                      aimag(F(:, lm)) .dot. (ChbR3(:KK4, :kN) * ll1(lm) - ChbD2R(:KK4, :kN)))

    ! For grad(T)
    Tr(lm, :) = cmplx(real(T(:, lm)) .dot. ChbD1(:KK2, :kN), &
                      aimag(T(:, lm)) .dot. ChbD1(:KK2, :kN))
    St(lm, :) = cmplx(real(T(:, lm)) .dot. ChbR(:KK2, :kN), &
                      aimag(T(:, lm)) .dot. ChbR(:KK2, :kN))

  end do  

end subroutine comp_radial_derivatives

!----------------------------------------------------------------------------

subroutine comp_U_cU_gT(Ur, Ut, Up, cUr, cUt, cUp, gTr, gTt, gTp)

  !###########################################################################
  !    Compute U, curl(U) and grat(T) from auxiliary fields Qu, Su and Tu
  !                       Qcu, Scu, Tcu, Tr and St
  !###########################################################################

  implicit none

  double precision, dimension(kN, lN, mN),  intent(out) :: Ur, Ut, Up ! Output real fields
  double precision, dimension(kN, lN, mN),  intent(out) :: cUr, cUt, cUp ! Output real fields
  double precision, dimension(kN, lN, mN),  intent(out) :: gTr, gTt, gTp ! Output real fields

  integer :: k

  Sh1 = 0. ; Sh2 = 0. ; Sh3 = 0. ;
  Ur = 0. ; Ut = 0. ; Up = 0. ;
  cUr = 0. ; cUt = 0. ; cUp = 0. ;
  gTr = 0. ; gTt = 0. ; gTp = 0. ;

  do k = 1, kN

    call SHqst_to_spat(shtns_c, Qu(:, k), Su(:, k), Tu(:, k), Sh1, Sh2, Sh3)
    Ur(k, :, :) = Sh1 ; Ut(k, :, :) = Sh2 ; Up(k, :, :) = Sh3 ;

    call SHqst_to_spat(shtns_c, Qcu(:, k), Scu(:, k), Tcu(:, k), Sh1, Sh2, Sh3)
    cUr(k, :, :) = Sh1 ; cUt(k, :, :) = Sh2 ; cUp(k, :, :) = Sh3 ;

    call SH_to_spat(shtns_c, Tr(:, k), Sh1)
    gTr(k, :, :) = Sh1

    call SHsph_to_spat(shtns_c, St(:, k), Sh1, Sh2)
    gTt(k, :, :) = Sh1 ; gTp(k, :, :) = Sh2 ;

  end do

end subroutine comp_U_cU_gT

!----------------------------------------------------------------------------

subroutine comp_U(Ur, Ut, Up)

  !###########################################################################
  !      Compute velocity components from auxiliary fields Qu, Su and Tu
  !###########################################################################

  implicit none

  double precision, dimension(kN, lN, mN),  intent(out) :: Ur, Ut, Up ! Output real fields

  integer :: k

  Sh1 = 0. ; Sh2 = 0. ; Sh3 = 0. ;
  Ur = 0. ; Ut = 0. ; Up = 0. ;

  do k = 1, kN

    call SHqst_to_spat(shtns_c, Qu(:, k), Su(:, k), Tu(:, k), Sh1, Sh2, Sh3)
    Ur(k, :, :) = Sh1 ; Ut(k, :, :) = Sh2 ; Up(k, :, :) = Sh3 ;

  end do

end subroutine comp_U

!----------------------------------------------------------------------------

subroutine comp_curlU(cUr, cUt, cUp)

  !###########################################################################
  !     Compute components of curl(U) auxiliary fields Qcu, Scu and Tcu
  !###########################################################################

  implicit none

  double precision, dimension(kN, lN, mN),  intent(out) :: cUr, cUt, cUp ! Output real fields

  integer :: k

  Sh1 = 0. ; Sh2 = 0. ; Sh3 = 0. ;
  cUr = 0. ; cUt = 0. ; cUp = 0. ;

  do k = 1, kN

    call SHqst_to_spat(shtns_c, Qcu(:, k), Scu(:, k), Tcu(:, k), Sh1, Sh2, Sh3)
    cUr(k, :, :) = Sh1 ; cUt(k, :, :) = Sh2 ; cUp(k, :, :) = Sh3 ;

  end do

end subroutine comp_curlU

!----------------------------------------------------------------------------

subroutine comp_gradT(gTr, gTt, gTp)

  !###########################################################################
  !         Compute components of grat(T) auxiliary fields Tr and St
  !###########################################################################

  implicit none

  double precision, dimension(kN, lN, mN),  intent(out) :: gTr, gTt, gTp ! Output real fields

  integer :: k

  Sh1 = 0. ; Sh2 = 0. ;
  gTr = 0. ; gTt = 0. ; gTp = 0. ;

  do k = 1, kN

    call SH_to_spat(shtns_c, Tr(:, k), Sh1)
    gTr(k, :, :) = Sh1

    call SHsph_to_spat(shtns_c, St(:, k), Sh1, Sh2)
    gTt(k, :, :) = Sh1 ; gTp(k, :, :) = Sh2 ;

  end do

end subroutine comp_gradT

!----------------------------------------------------------------------------

subroutine comp_U_from_EF(E, F, Ur, Ut, Up)

  !###########################################################################
  !                  Compute velocity components from E and F
  !###########################################################################

  implicit none

  double complex,   dimension(KK2, shtns%nlm), intent(in) :: E          ! Input e spectral field
  double complex,   dimension(KK4, shtns%nlm), intent(in) :: F          ! Input f spectral field

  double precision, dimension(kN, lN, mN),  intent(out) :: Ur, Ut, Up   ! Output real fields

  integer :: k, lm

  Sh1 = 0. ; Sh2 = 0. ; Sh3 = 0. ;
  Ur = 0. ; Ut = 0. ; Up = 0. ;

  do lm = 1, shtns%nlm

    ! For U
    Qu(lm, :) = cmplx(real(F(:, lm)) .dot. ChbR2(:KK4, :kN), &
                      aimag(F(:, lm)) .dot. ChbR2(:KK4, :kN)) * ll1(lm)
    Su(lm, :) = cmplx(real(F(:, lm)) .dot. ChbD1R(:KK4, :kN), &
                      aimag(F(:, lm)) .dot. ChbD1R(:KK4, :kN))
    Tu(lm, :) = cmplx(real(E(:, lm)) .dot. ChbR(:KK2, :kN), &
                      aimag(E(:, lm)) .dot. ChbR(:KK2, :kN))

  end do 

  do k = 1, kN

    call SHqst_to_spat(shtns_c, Qu(:, k), Su(:, k), Tu(:, k), Sh1, Sh2, Sh3)
    Ur(k, :, :) = Sh1 ; Ut(k, :, :) = Sh2 ; Up(k, :, :) = Sh3 ;

  end do

end subroutine comp_U_from_EF

!----------------------------------------------------------------------------

subroutine rCurl(Fr_spec, Gt_spec, Gp_spec, rCF_spec, rCCF_spec)

  !###########################################################################
  !     Compute r * e_r \cdot curl(F) and r^2 * e_r \cdot curl(curl(F))
  !###########################################################################

  implicit none
 
  double complex, dimension(KK, shtns%nlm), intent(in) :: Fr_spec, Gt_spec, Gp_spec ! Input spectral fields

  double complex, dimension(KK, shtns%nlm), intent(out) :: rCF_spec, rCCF_spec      ! Output spectral field

  double complex,   dimension(shtns%nlm) :: Slm, Slm_mul ! Intermediate arrays for SH transform

  double complex, dimension(KK, shtns%nlm) :: Gt_spec_th ! Intermediate array to contain (2 cos(theta) + sin(theta) * d/dtheta) Gt

  integer :: k, lm, m

  Slm = 0. ; Slm_mul = 0. ;
  Gt_spec_th = 0.
  rCF_spec = 0. ; rCCF_spec = 0. ;

  ! #### 1. We compute:
  ! - (2 cos(theta) + sin(theta) d/dtheta) Gp for  rCF
  ! - (2 cos(theta) + sin(theta) d/dtheta) Gt for rCCF
  do k = 1, KK
    ! For rCF
    Slm = Gp_spec(k, :)

    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_mul)
    rCF_spec(k, :) = Slm_mul
    
    call SH_mul_mx(shtns_c, costh_mul, Slm, Slm_mul)
    rCF_spec(k, :) = rCF_spec(k, :) + 2. * Slm_mul

    ! For rCCF
    Slm = Gt_spec(k, :)
    
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_mul)
    Gt_spec_th(k, :) = Slm_mul

    call SH_mul_mx(shtns_c, costh_mul, Slm, Slm_mul)
    Gt_spec_th(k, :) = Gt_spec_th(k, :) + 2. * Slm_mul
  end do

  ! #### 2. Final computation
  do lm = 1, shtns%nlm
    ! For rCF
    rCF_spec(:, lm) = rCF_spec(:, lm) + cmplx(aimag(Gt_spec(:, lm)), - real(Gt_spec(:, lm))) * dphi(lm) ! This is -d/dphi(G_t)
    
    ! For rCCF
    rCCF_spec(:, lm) = cmplx((Chb_mulR_deriv(:KK, :KK) .dot. real(Gt_spec_th(1:KK, lm))) + &
                        & - (Chb_mulR_deriv(:KK, :KK) .dot. aimag(Gp_spec(1:KK, lm)) * dphi(lm)) + &
                        & + (real(Fr_spec(:, lm)) * ll1(lm)), &
                        & (Chb_mulR_deriv(:KK, :KK) .dot. aimag(Gt_spec_th(1:KK, lm))) + &
                        & + (Chb_mulR_deriv(:KK, :KK) .dot. real(Gp_spec(1:KK, lm)) * dphi(lm)) + &
                        & + (aimag(Fr_spec(:, lm)) * ll1(lm))) 
  end do

  ! The following loop fixes a numerical issue and allows for correct computation of rCF up to the l mode LL (excluded)
  do m = 0, MM*mres, mres
    lm = shtns_lmidx(shtns_c, LL+1, m)
    rCF_spec(:, lm) = 0.
    rCCF_spec(:, lm) = 0.
  end do

end subroutine rCurl

!----------------------------------------------------------------------------

end module mod_ExplicitTerms
