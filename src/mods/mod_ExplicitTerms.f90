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

  call comp_U(E, F, Ur, Ut, Up) ! theta and phi components are multiplied by sin(theta)
  call comp_curlU(E, F, cUr, cUt, cUp) ! theta and phi components are multiplied by sin(theta)
  call ToReal(T, T_real, KK2)
  call comp_gradT(T, gTr, gTt, gTp) ! theta and phi components are multiplied by sin(theta)

  !--- Forcing terms for U

  ! Note 1: we write (U.grad)U = curl(U) x U + 1/2 * grad(U²)
  !       and there is no need to compute grad(U²) since we will take
  !       the curl of F later
  ! Note 2: As in comp_U, comp_curlU and comp_gradT the theta and phi
  !       components are multiplied by sin(theta) we need to take care of
  !       this when we compute the forcing terms

  do k = 1, kN
    do l = 1, lN
      Fr(k, l, :) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                  & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) &
                  & + 2.d0 * Up(k,l,:)
      Gt(k, l, :) = (Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) &
                  & + 2.d0 * CosTh(l) * Up(k,l,:)) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
      Gp(k, l, :) = (Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) &
                  & - 2.d0 * (CosTh(l)*Ut(k,l,:) +  SinTh(l)**2*Ur(k,l,:))) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
      UgradT(k, l, :) = - Ur(k,l,:) * gTr(k,l,:) &
                      & - (Ut(k,l,:) * gTt(k,l,:) + Up(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2)
    end do
  end do

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurlF(Gt_spec, Gp_spec, rCF_spec)
  call rCurlCurlF(Fr_spec, Gt_spec, Gp_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
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

  call comp_U(E, F, Ur, Ut, Up) ! theta and phi components are multiplied by sin(theta)
  call comp_curlU(E, F, cUr, cUt, cUp) ! theta and phi components are multiplied by sin(theta)
  call ToReal(T, T_real, KK2)
  call comp_gradT(T, gTr, gTt, gTp) ! theta and phi components are multiplied by sin(theta)
  
  do k = 1, kN
    do l = 1, lN
      Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                  & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1))
      Gt(k,l,:) = Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
      Gp(k,l,:) = Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
      UgradT(k,l,:) = - Ur(k,l,:) * gTr(k,l,:) &
                    & - (Ut(k,l,:) * gTt(k,l,:) + Up(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2)
    end do
  end do

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurlF(Gt_spec, Gp_spec, rCF_spec)
  call rCurlCurlF(Fr_spec, Gt_spec, Gp_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
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

  call comp_U(E, F, Ur, Ut, Up) ! theta and phi components are multiplied by sin(theta)
  call comp_curlU(E, F, cUr, cUt, cUp) ! theta and phi components are multiplied by sin(theta)
  call ToReal(T, T_real, KK2)
  call comp_gradT(T, gTr, gTt, gTp) ! theta and phi components are multiplied by sin(theta)

  !--- Differentiating by phi and multiplying by C_base
  call comp_cdphi(C_base, E_base, F_base, T_base, E_base_rot, F_base_rot, T_base_rot)
  ! E_base_rot = E_base; F_base_rot = F_base; T_base_rot = T_base;
  ! call comp_Rotation(- C_base, E_base_rot, KK2)
  ! call comp_Rotation(- C_base, F_base_rot, KK4)
  ! call comp_Rotation(- C_base, T_base_rot, KK2)
  call comp_U(E_base_rot, F_base_rot, Ur_rot, Ut_rot, Up_rot) ! theta and phi components are multiplied by sin(theta)
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
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) + Ur_rot(k,l,:) * Ek
        Gt(k,l,:) = Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:) + Ut_rot(k,l,:)) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:) + Up_rot(k,l,:)) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = - Ur(k,l,:) * gTr(k,l,:) &
                  & - (Ut(k,l,:) * gTt(k,l,:) + Up(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2) + T_real_rot(k,l,:)
      end do
    end do
  else 
    do k = 1, kN
      do l = 1, lN
        Fr(k, l, :) = Ek * (CUp(k,l,:)*Ut(k,l,:) - CUt(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                    & + Ra * T_real(k,l,:) / (Rout*rN(k)**(-1)) &
                    & + 2.d0 * Up(k,l,:) + Ur_rot(k,l,:) * Ek
        Gt(k, l, :) = (Ek * (CUr(k,l,:)*Up(k,l,:) - CUp(k,l,:)*Ur(k,l,:)) &
                    & + 2.d0 * CosTh(l) * Up(k,l,:) + Ut_rot(k,l,:) * Ek) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
        Gp(k, l, :) = (Ek * (CUt(k,l,:)*Ur(k,l,:) - CUr(k,l,:)*Ut(k,l,:)) &
                    & - 2.d0 * (CosTh(l)*Ut(k,l,:) +  SinTh(l)**2*Ur(k,l,:)) + Up_rot(k,l,:) * Ek) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k, l, :) = - Ur(k,l,:) * gTr(k,l,:) &
                        & - (Ut(k,l,:) * gTt(k,l,:) + Up(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2) + T_real_rot(k,l,:)
      end do
    end do
  end if

  !--- Going back to spectral space
  call BackToSpectral(Fr, Fr_spec, KK)    
  call BackToSpectral(Gt, Gt_spec, KK)
  call BackToSpectral(Gp, Gp_spec, KK)
  call BackToSpectral(UgradT, UgradT_spec, KK)

  !--- Curl(F) and Curl(Curl(F))
  call rCurlF(Gt_spec, Gp_spec, rCF_spec)
  call rCurlCurlF(Fr_spec, Gt_spec, Gp_spec, rCCF_spec)

  ! Now we allocate the RHS
  DE(:KK, :) = rCF_spec
  DF(:KK, :) = rCCF_spec
  DT(:KK, :) = UgradT_spec

  ! Then we fix the boundary conditions
  ! For the temperature
  DT(KK + 1, :) = 0.
  DT(KK + 2, :) = 0.
  DT(KK + 1, 1) = (1., 0.) ! To get T_inner = 1. and T_outer = 0. 
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

  call comp_U(E_base, F_base, Ur, Ut, Up)                  ! For U - theta and phi components are multiplied by sin(theta)
  call comp_curlU(E_base, F_base, cUr, cUt, cUp)           ! For U - theta and phi components are multiplied by sin(theta)
  call comp_U(E_per, F_per, Ur_per, Ut_per, Up_per)        ! For u - theta and phi components are multiplied by sin(theta)
  call comp_curlU(E_per, F_per, cUr_per, cUt_per, cUp_per) ! For u - theta and phi components are multiplied by sin(theta)
  call ToReal(T_base, T_real, KK2)                         ! Turn T to real
  call comp_gradT(T_base, gTr, gTt, gTp)                   ! For T - theta and phi components are multiplied by sin(theta)
  call ToReal(T_per, T_real_per, KK2)                      ! Turn t to real
  call comp_gradT(T_per, gTr_per, gTt_per, gTp_per)                 ! For t - theta and phi components are multiplied by sin(theta)

  !--- Differentiating by phi and multiplying by C_base
  call comp_cdphi(c_per, E_base, F_base, T_base, E_base_rot, F_base_rot, T_base_rot)
  call comp_cdphi(C_base, E_per, F_per, T_per, E_per_rot, F_per_rot, T_per_rot)
  ! E_base_rot = E_base; F_base_rot = F_base; T_base_rot = T_base;
  ! call comp_Rotation(- c_per, E_base_rot, KK2)
  ! call comp_Rotation(- c_per, F_base_rot, KK4)
  ! call comp_Rotation(- c_per, T_base_rot, KK2)
  ! E_per_rot = E_per; F_per_rot = F_per; T_per_rot = T_per;
  ! call comp_Rotation(- C_base, E_per_rot, KK2)
  ! call comp_Rotation(- C_base, F_per_rot, KK4)
  ! call comp_Rotation(- C_base, T_per_rot, KK2)
  call comp_U(E_base_rot, F_base_rot, Ur_rot, Ut_rot, Up_rot)
  call comp_U(E_per_rot, F_per_rot, Ur_per_rot, Ut_per_rot, Up_per_rot)
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
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut_per(k,l,:) - CUt(k,l,:)*Up_per(k,l,:)) / (SinTh(l)**2) &
                    & + Ek * (CUp_per(k,l,:)*Ut(k,l,:) - CUt_per(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                    & + Ra * T_real_per(k,l,:) / (Rout*rN(k)**(-1))  + Ek * (Ur_rot(k,l,:) + Ur_per_rot(k,l,:))
        Gt(k,l,:) = (Ek * (CUr(k,l,:)*Up_per(k,l,:) - CUp(k,l,:)*Ur_per(k,l,:)) &
                    & + Ek * (CUr_per(k,l,:)*Up(k,l,:) - CUp_per(k,l,:)*Ur(k,l,:)) & 
                    & + Ek * (Ut_rot(k,l,:) + Ut_per_rot(k,l,:))) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = (Ek * (CUt(k,l,:)*Ur_per(k,l,:) - CUr(k,l,:)*Ut_per(k,l,:)) &
                    & + Ek * (CUt_per(k,l,:)*Ur(k,l,:) - CUr_per(k,l,:)*Ut(k,l,:)) & 
                    & + Ek * (Up_rot(k,l,:) + Up_per_rot(k,l,:))) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = - Ur(k,l,:) * gTr_per(k,l,:) &
                      & - (Ut(k,l,:) * gTt_per(k,l,:) + Up(k,l,:) * gTp_per(k,l,:)) / (SinTh(l)**2) &
                      & - Ur_per(k,l,:) * gTr(k,l,:) &
                      & - (Ut_per(k,l,:) * gTt(k,l,:) + Up_per(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2) &
                      & + T_real_rot(k,l,:) + T_real_per_rot(k,l,:)
      end do
    end do
  else 
    do k = 1,kN
      do l = 1,lN
        Fr(k,l,:) = Ek * (CUp(k,l,:)*Ut_per(k,l,:) - CUt(k,l,:)*Up_per(k,l,:)) / (SinTh(l)**2) &
                  & + Ek * (CUp_per(k,l,:)*Ut(k,l,:) - CUt_per(k,l,:)*Up(k,l,:)) / (SinTh(l)**2) &
                  & + Ra * T_real_per(k,l,:) / (Rout*rN(k)**(-1)) &
                  & + 2.d0 * Up_per(k,l,:) + Ek * (Ur_rot(k,l,:) + Ur_per_rot(k,l,:))
        Gt(k,l,:) = (Ek * (CUr(k,l,:)*Up_per(k,l,:) - CUp(k,l,:)*Ur_per(k,l,:)) &
                  & + Ek * (CUr_per(k,l,:)*Up(k,l,:) - CUp_per(k,l,:)*Ur(k,l,:)) &
                  & + 2.d0 * CosTh(l) * Up_per(k,l,:) & 
                  & + Ek * (Ut_rot(k,l,:) + Ut_per_rot(k,l,:))) / (SinTh(l)**2) ! Here we compute Gt = Ft / sin(theta)
        Gp(k,l,:) = (Ek * (CUt(k,l,:)*Ur_per(k,l,:) - CUr(k,l,:)*Ut_per(k,l,:)) &
                  & + Ek * (CUt_per(k,l,:)*Ur(k,l,:) - CUr_per(k,l,:)*Ut(k,l,:)) &
                  & - 2.d0 * (CosTh(l)*Ut_per(k,l,:) +  SinTh(l)**2*Ur_per(k,l,:)) & 
                  & + Ek * (Up_rot(k,l,:) + Up_per_rot(k,l,:))) / (SinTh(l)**2) ! Here we compute Gp = Fp / sin(theta)
        UgradT(k,l,:) = - Ur(k,l,:) * gTr_per(k,l,:) &
                    & - (Ut(k,l,:) * gTt_per(k,l,:) + Up(k,l,:) * gTp_per(k,l,:)) / (SinTh(l)**2) &
                    & - Ur_per(k,l,:) * gTr(k,l,:) &
                    & - (Ut_per(k,l,:) * gTt(k,l,:) + Up_per(k,l,:) * gTp(k,l,:)) / (SinTh(l)**2) &
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
  call rCurlF(Gt_spec, Gp_spec, rCF_spec)
  call rCurlCurlF(Fr_spec, Gt_spec, Gp_spec, rCCF_spec)

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

subroutine comp_Rotation(alpha, U, KK)

  implicit none

  double precision, intent(in) :: alpha ! Angle

  integer, intent(in) :: KK  ! Number of Chebyshev modes

  double complex, dimension(KK, shtns%nlm), intent(inout) :: U ! E, E_per, F, F_per, T or Tt

  double complex, dimension(shtns%nlm) :: Slm, Slm_rot ! Intermediate array for SH rotation

  integer :: k

  ! Perform rotation by angle alpha
  do k = 1, KK
    Slm = U(k, :)
    call SH_Zrotate(shtns_c, Slm, alpha, Slm_rot)
    U(k, :) = Slm_rot
  end do

end subroutine comp_Rotation

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

subroutine comp_U(E, F, Ur, Ut, Up)

  !###########################################################################
  !         Compute velocity components from spectral fields e and f
  !###########################################################################

  implicit none

  double complex,   dimension(KK2, shtns%nlm), intent(in) :: E           ! Input e spectral field
  double complex,   dimension(KK4, shtns%nlm), intent(in) :: F           ! Input f spectral field

  double precision, dimension(kN, lN, mN),     intent(out) :: Ur, Ut, Up      ! Output real fields

  double precision, dimension(lN, mN) :: Sh                                   ! Intermediate array for SH transform
  double complex,   dimension(shtns%nlm) :: Slm, Slm_sth_dth                  ! Intermediate arrays for SH transform
  double complex,   dimension(kN, shtns%nlm) :: Ur_inter, Ut_inter, Up_inter  ! Intermediate fields

  double complex, dimension(KK4, shtns%nlm) :: Ut_spec, Up_spec  ! theta and phi velocity components in spectral space

  integer :: k, lm

  Sh = 0. ; Slm = 0. ; Slm_sth_dth = 0. ;
  Ur_inter = 0. ; Ut_inter = 0. ;  Up_inter = 0. ;
  Ut_spec = 0. ; Up_spec = 0. ;
  Ur = 0. ; Ut = 0. ; Up = 0. ;

  ! #### 1. We compute the sin(theta).d/dtheta constributions to Up_spec and Ut_spec
  do k = 1, KK2
    Slm = F(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    Ut_spec(k, :) = Slm_sth_dth

    Slm = E(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    Up_spec(k, :) = - Slm_sth_dth
  end do

  Slm = F(KK+3, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  Ut_spec(KK+3, :) = Slm_sth_dth

  Slm = F(KK+4, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  Ut_spec(KK+4, :) = Slm_sth_dth

  ! #### 2. We go to real in Chebyshev
  do lm = 1, shtns%nlm
    Ur_inter(:kN, lm)%re = (real(F(:, lm)) .dot. ChbR2(:KK4, :kN)) * ll1(lm)
    Ur_inter(:kN, lm)%im = (aimag(F(:, lm)) .dot. ChbR2(:KK4, :kN)) * ll1(lm)

    Ut_inter(:kN, lm)%re = (real(Ut_spec(:, lm)) .dot. ChbD1R(:KK4, :kN)) + &
                         & - ((aimag(E(:, lm)) * dphi(lm)) .dot. ChbR(:KK2, :kN))
    Ut_inter(:kN, lm)%im = (aimag(Ut_spec(:, lm)) .dot. ChbD1R(:KK4, :kN)) + &
                         & + ((real(E(:, lm)) * dphi(lm)) .dot. ChbR(:KK2, :kN))

    Up_inter(:kN, lm)%re = (real(Up_spec(:, lm)) .dot. ChbR(:KK4, :kN)) + &
                         & - ((aimag(F(:, lm)) * dphi(lm)) .dot. ChbD1R(:KK4, :kN))
    Up_inter(:kN, lm)%im = (aimag(Up_spec(:, lm)) .dot. ChbR(:KK4, :kN)) + &
                         & + ((real(F(:, lm)) * dphi(lm)) .dot. ChbD1R(:KK4, :kN))
  end do

  ! #### 3. We go to real in SH
  do k = 1, kN
    Slm = Ur_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    Ur(k, :, :) = Sh

    Slm = Ut_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    Ut(k, :, :) = Sh

    Slm = Up_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    Up(k, :, :) = Sh
  end do

end subroutine comp_U

!----------------------------------------------------------------------------

subroutine comp_curlU(E, F, cUr, cUt, cUp)

  !###########################################################################
  !         Compute velocity components from spectral fields e and f
  !###########################################################################

  implicit none

  double complex,   dimension(KK2, shtns%nlm), intent(in) :: E              ! Input e spectral field
  double complex,   dimension(KK4, shtns%nlm), intent(in) :: F              ! Input f spectral field

  double precision, dimension(kN, lN, mN),     intent(out) :: cUr, cUt, cUp ! Output real fields

  double precision, dimension(lN, mN) :: Sh                                     ! Intermediate array for SH transform
  double complex,   dimension(shtns%nlm) :: Slm, Slm_sth_dth                    ! Intermediate arrays for SH transform
  double complex,   dimension(kN, shtns%nlm) :: cUr_inter, cUt_inter, cUp_inter ! Intermediate fields

  double complex, dimension(KK2, shtns%nlm) :: E_r  ! Auxiliary field to contain de/dr in spectral space
  double complex, dimension(KK4, shtns%nlm) :: F_r2 ! Auxiliary field to contain d2f/dr2 in spectral space

  double complex, dimension(KK4, shtns%nlm) :: cUt_spec1, cUp_spec1, cUt_spec2, cUp_spec2 ! theta and phi components in spectral space

  integer :: k, lm

  Sh = 0. ; Slm = 0. ; Slm_sth_dth = 0. ;
  cUr_inter = 0. ; cUt_inter = 0. ; cUp_inter = 0. ;
  E_r = 0. ; F_r2 = 0.
  cUt_spec1 = 0. ; cUp_spec1 = 0. ; cUt_spec2 = 0. ; cUp_spec2 = 0. ;
  cUr = 0. ; cUt = 0. ; cUp = 0. ;

  ! #### 1. We compute de/dr and df/dr in spectral space and the laplacian and d/dphi contributions to cUp_spec and cUt_spec
  do lm = 1, shtns%nlm
    E_r(:, lm)%re = Chbderiv(0:KK2-1, 0:KK2-1) .dot. real(E(1:KK2, lm))
    E_r(:, lm)%im = Chbderiv(0:KK2-1, 0:KK2-1) .dot. aimag(E(1:KK2, lm))
    F_r2(:, lm)%re = Chbderiv2(0:KK4-1, 0:KK4-1) .dot. real(F(1:KK4, lm))
    F_r2(:, lm)%im = Chbderiv2(0:KK4-1, 0:KK4-1) .dot. aimag(F(1:KK4, lm))
    cUt_spec1(:, lm)%re = aimag(F_r2(:, lm)) * dphi(lm) ! This is -d/dphi(d2f/dr2)
    cUt_spec1(:, lm)%im = - real(F_r2(:, lm)) * dphi(lm)
    cUt_spec2(:, lm)%re = - aimag(F(:, lm)) * dphi(lm) * ll1(lm)
    cUt_spec2(:, lm)%im = real(F(:, lm)) * dphi(lm) * ll1(lm)
    cUp_spec1(:KK2, lm)%re = - aimag(E_r(:, lm)) * dphi(lm)
    cUp_spec1(:KK2, lm)%im = real(E_r(:, lm)) * dphi(lm)
    cUp_spec2(:, lm)%re = real(F(:, lm)) * ll1(lm)
    cUp_spec2(:, lm)%im = aimag(F(:, lm)) * ll1(lm)
  end do

  ! #### 2. We compute the sin(theta).d/dtheta constributions to cUt_spec and cUp_spec
  do k = 1, KK2
    Slm = E_r(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    cUt_spec1(k, :) = cUt_spec1(k, :) + Slm_sth_dth
    
    Slm = F_r2(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    cUp_spec1(k, :) = cUp_spec1(k, :) + Slm_sth_dth

    Slm = cUp_spec2(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    cUp_spec2(k, :) = - Slm_sth_dth
  end do

  Slm = F_r2(KK+3, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  cUp_spec1(KK+3, :) = cUp_spec1(KK+3, :) + Slm_sth_dth

  Slm = F_r2(KK+4, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  cUp_spec1(KK+4, :) = cUp_spec1(KK+4, :) + Slm_sth_dth

  Slm = cUp_spec2(KK+3, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  cUp_spec2(KK+3, :) = - Slm_sth_dth

  Slm = cUp_spec2(KK+4, :)
  call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
  cUp_spec2(KK+4, :) = - Slm_sth_dth

  ! #### 3. We go to real in Chebyshev
  do lm = 1, shtns%nlm
    cUr_inter(:kN, lm)%re = (real(E(:, lm)) .dot. ChbR2(:KK2, :kN)) * ll1(lm)
    cUr_inter(:kN, lm)%im = (aimag(E(:, lm)) .dot. ChbR2(:KK2, :kN)) * ll1(lm)

    cUt_inter(1:kN, lm)%re = (real(cUt_spec1(:, lm)) .dot. ChbR(:KK4, 1:kN)) + &
                           & (real(cUt_spec2(:, lm)) .dot. ChbR3(:KK4, 1:kN))

    cUt_inter(:kN, lm)%im = (aimag(cUt_spec1(:, lm)) .dot. ChbR(:KK4, :kN)) + &
                          & (aimag(cUt_spec2(:, lm)) .dot. ChbR3(:KK4, :kN))

    cUp_inter(:kN, lm)%re = (real(cUp_spec1(:, lm)) .dot. ChbR(:KK4, :kN)) + &
                          & (real(cUp_spec2(:, lm)) .dot. ChbR3(:KK4, :kN))

    cUp_inter(:kN, lm)%im = (aimag(cUp_spec1(:, lm)) .dot. ChbR(:KK4, :kN)) + &
                          & (aimag(cUp_spec2(:, lm)) .dot. ChbR3(:KK4, :kN))
  end do

  ! #### 4. We go to real in SH
  do k = 1, kN
    Slm = cUr_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    cUr(k, :, :) = Sh

    Slm = cUt_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    cUt(k, :, :) = Sh

    Slm = cUp_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    cUp(k, :, :) = Sh
  end do

end subroutine comp_curlU

!----------------------------------------------------------------------------

subroutine comp_gradT(T, gTr, gTt, gTp)

  !###########################################################################
  !         Compute velocity components from spectral fields e and f
  !###########################################################################

  implicit none

  double complex,   dimension(KK2, shtns%nlm), intent(in) :: T              ! Input T spectral field

  double precision, dimension(kN, lN, mN),     intent(out) :: gTr, gTt, gTp ! Output real fields

  double precision, dimension(lN, mN) :: Sh                                     ! Intermediate array for SH transform
  double complex,   dimension(shtns%nlm) :: Slm, Slm_sth_dth                    ! Intermediate arrays for SH transform
  double complex,   dimension(kN, shtns%nlm) :: gTr_inter, gTt_inter, gTp_inter ! Intermediate fields

  double complex, dimension(KK2, shtns%nlm) :: T_spec_r  ! Auxiliary field to contain dT/dr in spectral space

  double complex, dimension(KK2, shtns%nlm) :: gTt_spec, gTp_spec ! theta and phi components in spectral space

  integer :: k, lm

  Sh = 0. ; Slm = 0. ; Slm_sth_dth = 0. ;
  gTr_inter = 0. ; gTt_inter = 0. ; gTp_inter = 0. ;
  T_spec_r = 0.
  gTt_spec = 0. ; gTp_spec = 0. ;
  gTr = 0. ; gTt = 0. ; gTp = 0. ;

  ! #### 1. We compute dT/dr in spectral space and the d/dphi contributions to gTp_spec
  do lm = 1, shtns%nlm
    T_spec_r(:, lm)%re = Chbderiv(0:KK2-1, 0:KK2-1) .dot. real(T(1:KK2, lm))
    T_spec_r(:, lm)%im = Chbderiv(0:KK2-1, 0:KK2-1) .dot. aimag(T(1:KK2, lm))
    gTp_spec(:KK2, lm)%re = - aimag(T(:, lm)) * dphi(lm)
    gTp_spec(:KK2, lm)%im = real(T(:, lm)) * dphi(lm)
  end do

  ! #### 2. We compute sin(theta).dT/dtheta for the theta component
  do k = 1, KK2
    Slm = T(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_sth_dth)
    gTt_spec(k, :) = Slm_sth_dth
  end do

  ! #### 3. We go to real in Chebyshev
  do lm = 1, shtns%nlm
    gTr_inter(:kN, lm)%re = real(T_spec_r(:, lm)) .dot. Chb(:KK2, :kN)
    gTr_inter(:kN, lm)%im = aimag(T_spec_r(:, lm)) .dot. Chb(:KK2, :kN)

    gTt_inter(1:kN, lm)%re = real(gTt_spec(:, lm)) .dot. ChbR(:KK2, 1:kN)
    gTt_inter(:kN, lm)%im = aimag(gTt_spec(:, lm)) .dot. ChbR(:KK2, :kN)

    gTp_inter(:kN, lm)%re = real(gTp_spec(:, lm)) .dot. ChbR(:KK2, :kN)
    gTp_inter(:kN, lm)%im = aimag(gTp_spec(:, lm)) .dot. ChbR(:KK2, :kN)
  end do

  ! #### 4. We go to real in SH
  do k = 1, kN
    Slm = gTr_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    gTr(k, :, :) = Sh

    Slm = gTt_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    gTt(k, :, :) = Sh

    Slm = gTp_inter(k, :)
    call SH_to_spat(shtns_c, Slm, Sh)
    gTp(k, :, :) = Sh
  end do

end subroutine comp_gradT

!----------------------------------------------------------------------------

subroutine rCurlF(Gt_spec, Gp_spec, rCF_spec)

  !###########################################################################
  !                       Compute e_r \cdot curl(F)
  !###########################################################################

  implicit none
 
  double complex, dimension(KK, shtns%nlm), intent(in) :: Gt_spec, Gp_spec     ! Input spectral fields

  double complex, dimension(KK, shtns%nlm), intent(out) :: rCF_spec            ! Output spectral field

  double complex,   dimension(shtns%nlm) :: Slm, Slm_mul ! Intermediate arrays for SH transform

  integer :: k, lm, m

  Slm = 0. ; Slm_mul = 0. ;
  rCF_spec = 0.

  ! #### 1. We compute the d/dphi contributions to rCF
  do lm = 1, shtns%nlm
    rCF_spec(:, lm)%re = aimag(Gt_spec(:, lm)) * dphi(lm) ! This is -d/dphi(G_t)
    rCF_spec(:, lm)%im = - real(Gt_spec(:, lm)) * dphi(lm)
  end do

  ! #### 2. We compute the sin(theta).d/dtheta and multiplication by cos(theta) contributions to rCF
  do k = 1, KK
    Slm = Gp_spec(k, :)
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_mul)
    rCF_spec(k, :) = rCF_spec(k, :) + Slm_mul
    
    Slm = Gp_spec(k, :)
    call SH_mul_mx(shtns_c, costh_mul, Slm, Slm_mul)
    rCF_spec(k, :) = rCF_spec(k, :) + 2. * Slm_mul
  end do

  do m = 1, MM
    lm = shtns_lmidx(shtns_c, LL+1, 1)
    rCF_spec(:, lm) = 0.
  end do

end subroutine rCurlF

subroutine rCurlCurlF(Fr_spec, Gt_spec, Gp_spec, rCCF_spec)

  !###########################################################################
  !                       Compute e_r \cdot curl(F)
  !###########################################################################

  implicit none
 
  double complex, dimension(KK, shtns%nlm), intent(in) :: Fr_spec, Gt_spec, Gp_spec ! Input spectral fields

  double complex, dimension(KK, shtns%nlm), intent(out) :: rCCF_spec                ! Output spectral field

  double complex, dimension(shtns%nlm) :: Slm, Slm_mul ! Intermediate arrays for SH transform

  double complex, dimension(KK, shtns%nlm) :: Gt_spec_th ! Intermediate array to contain (2 cos(theta) + sin(theta) * d/dtheta) Gt

  integer :: k, lm, m

  Slm = 0. ; Slm_mul = 0. ;
  Gt_spec_th = 0.
  rCCF_spec = 0.

  ! #### 1. We compute the (2 cos(theta) + sin(theta) * d/dtheta) Gt
  do k = 1, KK
    Slm = Gt_spec(k, :)
    
    call SH_mul_mx(shtns_c, sth_dth, Slm, Slm_mul)
    Gt_spec_th(k, :) = Slm_mul

    call SH_mul_mx(shtns_c, costh_mul, Slm, Slm_mul)
    Gt_spec_th(k, :) = Gt_spec_th(k, :) + 2. * Slm_mul
  end do

  ! #### 2. We compute the d/dphi and d/dr and laplacian
  do lm = 1, shtns%nlm
    rCCF_spec(:, lm)%re = (Chb_mulR_deriv(:KK, :KK) .dot. real(Gt_spec_th(1:KK, lm))) + &
                        & - (Chb_mulR_deriv(:KK, :KK) .dot. aimag(Gp_spec(1:KK, lm)) * dphi(lm)) + &
                        & + (real(Fr_spec(:, lm)) * ll1(lm))
    rCCF_spec(:, lm)%im = (Chb_mulR_deriv(:KK, :KK) .dot. aimag(Gt_spec_th(1:KK, lm))) + &
                        & + (Chb_mulR_deriv(:KK, :KK) .dot. real(Gp_spec(1:KK, lm)) * dphi(lm)) + &
                        & + (aimag(Fr_spec(:, lm)) * ll1(lm))
  end do

  do m = 1, MM
    lm = shtns_lmidx(shtns_c, LL+1, 1)
    rCCF_spec(:, lm) = 0.
  end do

end subroutine rCurlCurlF

end module mod_ExplicitTerms
