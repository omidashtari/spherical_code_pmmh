module mod_precompXY
 !$ USE OMP_LIB
  use mod_Globalvars
  use mod_Tchebyshev
  use mod_Matrices
  use mod_PrecompSH

  implicit none

contains

subroutine precompBuildXY()

  implicit none

  integer :: k

  ! This subroutine will construct the building blocs of the different matrices of the solve.

  allocate(Chb0XY(KK, KK4), ChbD1XY(KK, KK4), ChbD2XY(KK, KK4), ChbD4XY(KK, KK4), divR(KK, KK4))
  allocate(BC_mat(KK4, 4))

  do k = 1, KK4
    Chb0XY(:, k) = ChbPoly(k - 1, xK)
    ChbD1XY(:, k) = ChbPolDeriv(k - 1, xK, 1) * dxdr
    ChbD2XY(:, k) = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
    ChbD4XY(:, k) = ChbPolDeriv(k - 1, xK, 4) * dxdr ** 4
    divR(:, k) = 1 / rK
    ! Boundary conditions
    BC_mat(k, 1) = (-1.) ** (k - 1)
    BC_mat(k, 2) =   1.
    BC_mat(k, 3) = (-1.) ** k * (k - 1.) ** 2
    BC_mat(k, 4) =   1. * (k - 1.) ** 2
  end do

end subroutine precompBuildXY

! ----------------------------------------------------------------------------------------
! -----------------------------  Creation of blocks --------------------------------------

! Remark: The E equaiton is multiplied by r on both sides and the F one by r^2

function create_Xe(l, param) result(Xe)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param

  double precision, dimension(KK2, KK2) :: Xe

  Xe = 0.

  Xe(:KK, :) = dble(l * (l + 1)) * divR(:, :KK2) * (param * Chb0XY(:, :KK2) - IER * delta_t * Ek * &
               & (ChbD2XY(:, :KK2) - dble(l * (l + 1)) * divR(:, :KK2) ** 2 * Chb0XY(:, :KK2)))

  !Boundary conditions
  Xe(KK + 1, :) = BC_mat(:KK2, 1)
  Xe(KK + 2, :) = BC_mat(:KK2, 2) 
  
end function create_Xe

! ---------------------------------------------------------------------------------------

function create_XTh(l) result(XTh)

  implicit none

  integer, intent(in) :: l

  double precision, dimension(KK2, KK2) :: XTh

  XTh = 0.

  XTh(:KK, :) = Chb0XY(:, :KK2) - IER * delta_t * (ChbD2XY(:, :KK2) + 2.d0 * divR(:, :KK2) * ChbD1XY(:, :KK2) &
              & - l * (l + 1.d0) * divR(:, :KK2) ** 2 * Chb0XY(:, :KK2)) / Pr

  ! Boundary conditions
  XTh(KK + 1, :) = BC_mat(:KK2, 1)
  XTh(KK + 2, :) = BC_mat(:KK2, 2)
  
end function create_XTh

! ---------------------------------------------------------------------------------------

function create_Xf(l, param) result(Xf)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param

  double precision, dimension(KK4, KK4) :: Xf
  double precision, dimension(KK, KK4)  :: fac

  Xf = 0.
  fac = dble(l * (l + 1)) * divR ** 2

  Xf(:KK, :) = - dble(l * (l + 1)) * (param * (ChbD2XY - fac * Chb0XY) - IER * delta_t * Ek * (ChbD4XY - 2. * fac * ChbD2XY &
               & + 4. * fac * divR * ChbD1XY + fac ** 2 * Chb0XY - 6. * fac * divR ** 2 * Chb0XY))

  !Boundary conditions
  Xf(KK + 1, :) = BC_mat(:, 1)
  Xf(KK + 2, :) = BC_mat(:, 2)
  Xf(KK + 3, :) = BC_mat(:, 3)
  Xf(KK + 4, :) = BC_mat(:, 4)

end function create_Xf

! -----------------------------------------------------------------------------------

function create_Be(m, l, COEF) result(Be)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK2, KK2) :: Be
  
  Be = 0.

  Be(:KK, :) = - delta_t * COEF * 2 * dble(m) * divR(:, :KK2) * Chb0XY(:, :KK2)

end function create_Be

! ----------------------------------------------------------------------------------------

function create_Bf(m, l, COEF) result(Bf)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK4, KK4) :: Bf

  Bf = 0.

  Bf(:KK, :) =  delta_t * COEF * 2. * dble(m) * (ChbD2XY - dble(l * (l + 1)) * divR ** 2 * Chb0XY)   

end function create_Bf

! ----------------------------------------------------------------------------------------

function create_Xue(m, l, COEF) result(Xue)
    
  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK4, KK2) :: Xue
   
  Xue = 0.

  if ((l>=m) .and. (l + 1 > 1)) then

    Xue(:KK, :) = - delta_t * 2. * COEF * dble(l * (l + abs(m) + 1) * (l + 2)) / dble(2 * l + 3) * &
                & (dble(l + 1) * Chb0XY(:, :KK2) * divR(:, :KK2) + ChbD1XY(:, :KK2)) * &
                & sqrt(dble(2*l+3)/dble(2*l+1)*dble(l+1-m)/dble(l+1+m))

  end if
  
end function create_Xue

! ----------------------------------------------------------------------------------------

function create_Xuf(m, l, COEF) result(Xuf)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK2, KK4) :: Xuf
  
  Xuf = 0.

  if ((l>=m) .and. (l + 1 > 1)) then

    Xuf(:KK, :) = - delta_t * 2. * COEF * divR * dble(l * (l + abs(m) + 1) * (l + 2)) / dble(2 * l + 3) * &
                & (dble(l + 1) * Chb0XY * divR + ChbD1XY) * sqrt(dble(2*l+3)/dble(2*l+1)*dble(l+1-m)/dble(l+1+m))

  end if
   
end function create_Xuf

! ----------------------------------------------------------------------------------------

function create_Xle(m, l, COEF) result(Xle)
        
  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK4, KK2) :: Xle

  Xle = 0.

  if (l - 1 < LL) then

    Xle(:KK, :) = delta_t * 2. * COEF * dble((l + 1) * (l - 1) * (l - abs(m))) / dble(2 * l - 1) * &
                & (dble(l) * Chb0XY(:, :KK2) * divR(:, :KK2) - ChbD1XY(:, :KK2)) * &
                & sqrt(dble(2*l-1)/dble(2*l+1)*dble(l+m)/dble(l-m))

  end if
   
end function create_Xle

! ----------------------------------------------------------------------------------------

function create_Xlf(m, l, COEF) result(Xlf)
     
  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: COEF

  double precision, dimension(KK2,KK4) :: Xlf

  Xlf = 0.

  if (l - 1 < LL) then
    
    Xlf(:KK, :) = delta_t * 2. * COEF * divR * dble((l + 1) * (l - 1) * (l - abs(m))) / dble(2 * l - 1) * &
                  & (dble(l) * Chb0XY * divR - ChbD1XY) * sqrt(dble(2*l-1)/dble(2*l+1)*dble(l+m)/dble(l-m))

  end if
 
end function create_Xlf

! ----------------------------------------------------------------------------------------

function create_Ye(l, param) result(Ye)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param

  double precision, dimension(KK2, KK2) :: Ye

  Ye = 0.

  Ye(:KK, :) = dble(l * (l + 1)) * divR(:, :KK2) * (param * Chb0XY(:, :KK2) + (1- IER) * delta_t * Ek * ( &
              & ChbD2XY(:, :KK2) - dble(l * (l + 1)) * divR(:, :KK2) ** 2 * Chb0XY(:, :KK2)))

end function create_Ye

! ----------------------------------------------------------------------------------------

function create_Ye_BDF(l, param) result(Ye)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param

  double precision, dimension(KK2, KK2) :: Ye

  Ye = 0.

  Ye(:KK, :) = dble(l * (l + 1)) * divR(:, :KK2) * param * Chb0XY(:, :KK2)

end function create_Ye_BDF

! ----------------------------------------------------------------------------------------

function create_YTh(l) result(YTh)

  implicit none

  integer, intent(in) :: l

  double precision, dimension(KK2, KK2) :: YTh

  YTh = 0.

  YTh(:KK, :) = Chb0XY(:, :KK2) + (1.d0 - IER) * delta_t * (ChbD2XY(:, :KK2) + 2.d0 * divR(:, :KK2) * ChbD1XY(:, :KK2) &
              & - l * (l + 1.d0) * divR(:, :KK2) ** 2 * Chb0XY(:, :KK2)) / Pr

end function create_YTh

! ----------------------------------------------------------------------------------------

function create_YTh_BDF(l) result(YTh)

  implicit none

  integer, intent(in) :: l

  double precision, dimension(KK2, KK2) :: YTh

  YTh = 0.

  YTh(:KK, :) = Chb0XY(:, :KK2)

end function create_YTh_BDF

! ----------------------------------------------------------------------------------------

function create_Yf(l, param) result(Yf)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param
  
  double precision, dimension(KK4, KK4) :: Yf
  double precision, dimension(KK, KK4) :: fac

  Yf = 0.

  fac = dble(l * (l + 1)) * divR ** 2

  Yf(:KK, :) = - dble(l * (l + 1)) * (param * (ChbD2XY - fac * Chb0XY) + & 
            & (1 - IER) * delta_t * Ek * (ChbD4XY - 2. * fac * ChbD2XY + &
            & 4. * fac * divR * ChbD1XY + fac ** 2 * Chb0XY - 6. * fac * divR ** 2 * Chb0XY))
   
end function create_Yf

! ----------------------------------------------------------------------------------------

function create_Yf_BDF(l, param) result(Yf)

  implicit none

  integer, intent(in) :: l
  double precision, intent(in) :: param
  
  double precision, dimension(KK4, KK4) :: Yf
  double precision, dimension(KK, KK4) :: fac

  Yf = 0.

  fac = dble(l * (l + 1)) * divR ** 2

  Yf(:KK, :) = - dble(l * (l + 1)) * param * (ChbD2XY - fac * Chb0XY)
   
end function create_Yf_BDF

! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! -----------------------------  Explicit Coriolis ---------------------------------------

subroutine precompXeYe()

  implicit none

  double precision, dimension(KK2, KK2) :: Xe, Ye
  double precision :: param
  integer :: l

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! The l=m=0 mode is of no use
    
    Xe = create_Xe(l, param)
    Ye = create_Ye(l, param)

    Xe_inv  (:, :, l) = invL(Xe)
    Xe_invYe(:, :, l) = transpose(Xe_inv(:, :, l) .dot. Ye)
    Xe_inv  (:, :, l) = transpose(Xe_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXeYe
   
! -------------------------------------------------------------------------

subroutine precompXeYe_BDF()

  implicit none

  double precision, dimension(KK2, KK2) :: Xe, Ye
  double precision :: param
  integer :: l

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! The l=m=0 mode is of no use
    
    Xe = create_Xe(l, param)
    Ye = create_Ye_BDF(l, param)

    Xe_inv  (:, :, l) = invL(Xe)
    Xe_invYe(:, :, l) = transpose(Xe_inv(:, :, l) .dot. Ye)
    Xe_inv  (:, :, l) = transpose(Xe_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXeYe_BDF
   
! -------------------------------------------------------------------------

subroutine precompXfYf()

  implicit none

  double precision, dimension(KK4, KK4) :: Xf, Yf
  double precision :: param
  integer :: l

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! The l=m=0 mode is of no use

    Xf = create_Xf(l, param)
    Yf = create_Yf(l, param)
    
    Xf_inv  (:, :, l) = invL(Xf)
    Xf_invYf(:, :, l) = transpose(Xf_inv(:, :, l) .dot. Yf)
    Xf_inv  (:, :, l) = transpose(Xf_inv(:, :, l) .dot. Chb4)

  end do

end subroutine precompXfYf

!-----------------------------------------------------------------------------

subroutine precompXfYf_BDF()

  implicit none

  double precision, dimension(KK4, KK4) :: Xf, Yf
  double precision :: param
  integer :: l

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! The l=m=0 mode is of no use

    Xf = create_Xf(l, param)
    Yf = create_Yf_BDF(l, param)
    
    Xf_inv  (:, :, l) = invL(Xf)
    Xf_invYf(:, :, l) = transpose(Xf_inv(:, :, l) .dot. Yf)
    Xf_inv  (:, :, l) = transpose(Xf_inv(:, :, l) .dot. Chb4)

  end do

end subroutine precompXfYf_BDF

!-----------------------------------------------------------------------------

subroutine precompXTYT()

  implicit none

  double precision, dimension(KK2, KK2) :: XTh, YTh
  integer :: l

  do l = 0, LL + 1 ! We must include the l=m=0 mode - cf Hollerback's 2000 paper.

    XTh = create_XTh(l)
    YTh = create_YTh(l)

    XT_inv  (:, :, l) = invL(XTh)
    XT_invYT(:, :, l) = transpose(XT_inv(:, :, l) .dot. YTh)
    XT_inv  (:, :, l) = transpose(XT_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXTYT

!-----------------------------------------------------------------------------

subroutine precompXTYT_BDF()

  implicit none

  double precision, dimension(KK2, KK2) :: XTh, YTh
  integer :: l

  do l = 0, LL + 1 ! We must include the l=m=0 mode - cf Hollerback's 2000 paper.

    XTh = create_XTh(l)
    YTh = create_YTh_BDF(l)

    XT_inv  (:, :, l) = invL(XTh)
    XT_invYT(:, :, l) = transpose(XT_inv(:, :, l) .dot. YTh)
    XT_inv  (:, :, l) = transpose(XT_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXTYT_BDF

! ----------------------------------------------------------------------------------------
! -----------------------------  Implicit Coriolis ---------------------------------------


! Creation of 2lU(KK2+KK4)x2lU(KK2+KK4) banded matrices X and Y using blocks :
!
!         |Xe    Be    Xuf   0                                                                     |
!         |Be    Xe    0     Xuf                                                                   |
!         |Xle   0     Xf    Bf    Xue   0                                                         |
!         |0     Xle   Bf    Xf    0     Xue                                                       |
!         |            Xlf   0     Xe    Be   Xuf  0                                               |
!         |            0     Xlf   Be    Xe   0    Xuf                                             |
!         |                        Xle   0    Xf   Bf                                              |
!         |                        0     Xle  Bf   Xf                                              | 
!         |                                            Xf    Bf    Xue   0                         |
!         |                                            Bf    Xf    0     Xue                       |
!         |                                            Xlf   0     Xe    Be    Xuf   0             | 
!         |                                            0     Xlf   Be    Xe    0     Xuf           |
!         |                                                        Xle   0     Xf    Bf   Xue  0   |
!         |                                                        0     Xle   Bf    Xf   0    Xue |
!         |                                                                    Xlf   0    Xe   Be  |
!         |                                                                    0     Xlf  Be   Xe  |
!         |
          

! ----------------------------------------------------------------------------------------
! -----------------------------  Creation of columns -------------------------------------

function CXF(l, m, param)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: param

  double precision, dimension(4 * KK2 + 2 * KK4, 2 * KK4)::CXF
  integer :: step

  step = 2 * (KK2 + KK4)

  CXF = 0.

  CXF(1 : KK2, 1 : KK4) = create_Xuf(m, l - 1, IER)
  CXF(KK2 + 1 : 2 * KK2, KK4 + 1 : 2 * KK4) = create_Xuf(m, l - 1, IER)
  CXF(2 * KK2 + 1 : 2 * KK2 + KK4, 1 : KK4) = create_Xf(l, param)
  CXF(2 * KK2 + 1 : 2 * KK2 + KK4, KK4 + 1 : 2 * KK4) = - create_Bf(m, l, IER)
  CXF(2 * KK2 + KK4 + 1 : step, 1 : KK4) = create_Bf(m, l, IER)
  CXF(2 * KK2 + KK4 + 1 : step, KK4 + 1 : 2 * KK4) = create_Xf(l, param)
  CXF(1 + step : step + KK2, 1 : KK4) = create_Xlf(m, l + 1, IER)
  CXF(step + KK2 + 1 : step + 2 * KK2, KK4 + 1 : 2 * KK4) = create_Xlf(m, l + 1, IER)

end function CXF

! ----------------------------------------------------------------------------------------
   
function CYF(l, m, param)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: param

  double precision, dimension(4 * KK2 + 2 * KK4, 2 * KK4) :: CYF
  integer :: step

  step = 2 * (KK2 + KK4)

  CYF = 0.

  CYF(1 : KK2, 1 : KK4) = - create_Xuf(m, l - 1, (1 - IER))
  CYF(KK2 + 1 : 2 * KK2, KK4 + 1 : 2 * KK4) = - create_Xuf(m, l - 1, (1 - IER))
  CYF(2 * KK2 + 1 : 2 * KK2 + KK4, 1 : KK4) = create_Yf(l, param)
  CYF(2 * KK2 + 1 : 2 * KK2 + KK4, KK4 + 1 : 2 * KK4) = create_Bf(m, l, (1 - IER))
  CYF(2 * KK2 + KK4 + 1 : step, 1 : KK4) = - create_Bf(m, l, (1 - IER))
  CYF(2 * KK2 + KK4 + 1 : step, KK4 + 1 : 2 * KK4) = create_Yf(l, param)
  CYF(1 + step : step + KK2, 1 : KK4) = - create_Xlf(m, l + 1, (1 - IER))
  CYF(step + KK2 + 1 : step + 2 * KK2, KK4 + 1 : 2 * KK4) = - create_Xlf(m, l + 1, (1 - IER))

end function CYF

! -----------------------------------------------------------------------------------------

function CXE(l, m, param)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: param

  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: CXE
  integer :: step

  step = 2 * (KK2 + KK4)

  CXE = 0.

  CXE(1 : KK4, 1 : KK2) = create_Xue(m, l - 1, IER)                                           
  CXE(KK4 + 1 : 2 * KK4, KK2 + 1 : 2 * KK2) = create_Xue(m, l - 1, IER)                               
  CXE(2 * KK4 + 1 : 2 * KK4 + KK2, 1 : KK2) = create_Xe(l, param)                                   
  CXE(2 * KK4 + 1 : 2 * KK4 + KK2, KK2 + 1 : 2 * KK2) = - create_Be(m, l, IER)                          
  CXE(2 * KK4 + KK2 + 1 : step, 1 : KK2)= create_Be(m, l, IER)                                   
  CXE(2 * KK4 + KK2 + 1 : step, KK2 + 1 : 2 * KK2) = create_Xe(l, param)                                 
  CXE(1 + step : step + KK4, 1 : KK2) = create_Xle(m, l + 1, IER)                                       
  CXE(step + KK4 + 1 : step + 2 * KK4, KK2 + 1 : 2 * KK2) = create_Xle(m, l + 1, IER)

end function CXE

! ----------------------------------------------------------------------------------------

function CYE(l, m, param)

  implicit none

  integer, intent(in) :: l, m
  double precision, intent(in) :: param

  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: CYE
  integer :: step

  step = 2 * (KK2 + KK4)

  CYE = 0.

  CYE(1 : KK4, 1 : KK2) = - create_Xue(m, l - 1, (1 - IER))
  CYE(KK4 + 1 : 2 * KK4, KK2 + 1 : 2 * KK2) = - create_Xue(m, l - 1, (1 - IER))
  CYE(2 * KK4 + 1 : 2 * KK4 + KK2, 1 : KK2) = create_Ye(l, param)
  CYE(2 * KK4 + 1 : 2 * KK4 + KK2, KK2 + 1 : 2 * KK2) = create_Be(m, l, (1 - IER))
  CYE(2 * KK4 + KK2 + 1 : step , 1 : KK2) = - create_Be(m, l, (1 - IER))
  CYE(2 * KK4 + KK2 + 1 : step, KK2 + 1 : 2 * KK2) = create_Ye(l, param)
  CYE(1 + step : step + KK4, 1 : KK2) = - create_Xle(m, l + 1, (1 - IER))
  CYE(step + KK4 + 1 : step + 2 * KK4, KK2 + 1 : 2 * KK2) = - create_Xle(m, l + 1, (1 - IER))

end function CYE

! ----------------------------------------------------------------------------------------
! ---Creation of matrices X and Y in the form allowing LU-banded-matrice-decomposition----
! -------------------------and banded-matrice-multiplication------------------------------
   
subroutine PrecompimplicitXY()

  implicit none

  integer :: step, mid, top, i, j, l, m, p
  integer :: l_odd, l_even, n_even, n_odd, m_idx
  double precision, dimension(2 * KK4 + 4 * KK2, 2 * KK4) :: ColumnFX, ColumnFY
  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: ColumnEX, ColumnEY
  double precision :: param

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  Xef = 0.0d0
  Yef = 0.0d0
  PIVOT = 0

  step = 2 * (KK2 + KK4)

  m_idx = 0
    
  do m = 0, MM*mres, mres

    ColumnEX = 0.
    ColumnFX = 0.
    ColumnEY = 0.
    ColumnFY = 0.

    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    top = 2 * (KK2 + KK4) * (n_even + n_odd)
    p = 1 - mod(max(m, 1), 2) ! Indicates if we start with the even or the odd

    l_odd = 0
    l_even = 0

    do l = (m + mod(m + 1, 2)), LL, 2 ! Loop over odds. The (m + mod(m + 1, 2)) is because if m=4 the loop begins in l=5
      ColumnEX = CXE(l, m, param)
      ColumnFX = CXF(l, m, param)
      ColumnEY = CYE(l, m, param)
      ColumnFY = CYF(l, m, param)

      do i = 1, 2 * KK2 + 4 * KK4 ! Odds of e are in the first symmetry group
        do j = 1, 2 * KK2
          Xef(step + i - j + 2 * KK2 - 1, 2 * KK4 * p + l_odd * step + j, m_idx) = ColumnEX(i, j)
          Yef(i - j + 2 * KK2 - 1 + 1, 2 * KK4 * p + l_odd * step + j, m_idx) = ColumnEY(i, j)
        end do
      end do

      do i = 1, 2 * KK4 + 4 * KK2 ! Odds of f are in the second symmetry group
        do j = 1, 2 * KK4
          Xef(step + i - j + 2 * KK4 - 1, 2 * KK2 * p + mid + l_odd * step + j, m_idx) = ColumnFX(i, j) 
          Yef(i - j + 2 * KK4 - 1 + 1, 2 * KK2 * p + mid + l_odd * step + j, m_idx) = ColumnFY(i, j)
        end do
      end do

      l_odd = l_odd + 1
    end do

    do l = (max(m, 1) + mod(max(m, 1), 2)), LL, 2 ! Loop over evens. Same thing goes for the first l as in the odd case.
      ColumnFX = CXF(l, m, param)
      ColumnEX = CXE(l, m, param)
      ColumnFY = CYF(l, m, param)
      ColumnEY = CYE(l, m, param)

      do i = 1, 2 * KK4 + 4 * KK2 ! Evens of f are in the first symmetry group
        do j = 1, 2 * KK4
          Xef(step + i - j + 2 * KK4 - 1, 2 * KK2 * (1 - p) + l_even * step + j, m_idx) = ColumnFX(i, j)
          Yef(i - j + 2 * KK4 - 1 + 1, 2 * KK2 * (1 - p) + l_even * step + j, m_idx) = ColumnFY(i, j)
        end do
      end do

      do i = 1, 2 * KK2 + 4 * KK4 ! Evens of e are in the second symmetry group
        do j = 1, 2 * KK2
          Xef(step + i - j + 2 * KK2 - 1, 2 * KK4 * (1 - p) + mid + l_even * step + j, m_idx) = ColumnEX(i, j) 
          Yef(i - j + 2 * KK2 - 1 + 1, 2 * KK4 * (1 - p) + mid + l_even * step + j, m_idx) = ColumnEY(i, j)
        end do
      end do

      l_even = l_even + 1
    end do

    call Banded_LU_decomp(top, 2 * (KK2 + KK4) - 1, 6 * (KK2 + KK4) - 2, Xef(:, :top, m_idx), PIVOT(:top, m_idx))

    m_idx = m_idx + 1

  end do

end subroutine PrecompimplicitXY

subroutine PrecompimplicitXY_BDF()

  implicit none

  integer :: step, mid, top, i, j, l, m, p
  integer :: l_odd, l_even, n_even, n_odd, m_idx
  double precision, dimension(2 * KK4 + 4 * KK2, 2 * KK4) :: ColumnFX
  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: ColumnEX
  double precision :: param

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  Xef = 0.0d0
  PIVOT = 0

  step = 2 * (KK2 + KK4)

  m_idx = 0

  ! First we compute Xef
    
  do m = 0, MM*mres, mres

    ColumnEX = 0.
    ColumnFX = 0.

    n_even = ((LL + mod(LL, 2)) - max(m, 1)) / 2 + mod(LL + 1, 2)
    n_odd = ((LL + mod(LL, 2)) + 1 - m) / 2
    mid = 2 * (KK2 * n_odd + KK4 * n_even)
    top = 2 * (KK2 + KK4) * (n_even + n_odd)
    p = 1 - mod(max(m, 1), 2) ! Indicates if we start with the even or the odd

    l_odd = 0
    l_even = 0

    do l = (m + mod(m + 1, 2)), LL, 2 ! Loop over odds. The (m + mod(m + 1, 2)) is because if m=4 the loop begins in l=5
      ColumnEX = CXE(l, m, param)
      ColumnFX = CXF(l, m, param)

      do i = 1, 2 * KK2 + 4 * KK4 ! Odds of e are in the first symmetry group
        do j = 1, 2 * KK2
          Xef(step + i - j + 2 * KK2 - 1, 2 * KK4 * p + l_odd * step + j, m_idx) = ColumnEX(i, j)
        end do
      end do

      do i = 1, 2 * KK4 + 4 * KK2 ! Odds of f are in the second symmetry group
        do j = 1, 2 * KK4
          Xef(step + i - j + 2 * KK4 - 1, 2 * KK2 * p + mid + l_odd * step + j, m_idx) = ColumnFX(i, j)
        end do
      end do

      l_odd = l_odd + 1
    end do

    do l = (max(m, 1) + mod(max(m, 1), 2)), LL, 2 ! Loop over evens. Same thing goes for the first l as in the odd case.
      ColumnFX = CXF(l, m, param)
      ColumnEX = CXE(l, m, param)

      do i = 1, 2 * KK4 + 4 * KK2 ! Evens of f are in the first symmetry group
        do j = 1, 2 * KK4
          Xef(step + i - j + 2 * KK4 - 1, 2 * KK2 * (1 - p) + l_even * step + j, m_idx) = ColumnFX(i, j)
        end do
      end do

      do i = 1, 2 * KK2 + 4 * KK4 ! Evens of e are in the second symmetry group
        do j = 1, 2 * KK2
          Xef(step + i - j + 2 * KK2 - 1, 2 * KK4 * (1 - p) + mid + l_even * step + j, m_idx) = ColumnEX(i, j)
        end do
      end do

      l_even = l_even + 1
    end do

    call Banded_LU_decomp(top, 2 * (KK2 + KK4) - 1, 6 * (KK2 + KK4) - 2, Xef(:, :top, m_idx), PIVOT(:top, m_idx))

    m_idx = m_idx + 1

  end do

  ! And then we do Ye_mat and Yf_mat

  Ye_mat = 0.
  Yf_mat = 0.

  do l = 1, LL + 1 ! The l=m=0 mode is of no use
    Ye_mat(:, :, l) = create_Ye_BDF(l, param)
    Yf_mat(:, :, l) = create_Yf_BDF(l, param)
  end do

end subroutine PrecompimplicitXY_BDF

end module mod_precompXY
