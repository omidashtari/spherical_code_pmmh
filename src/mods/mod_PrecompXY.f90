module mod_precompXY
 !$ USE OMP_LIB
  use mod_Globalvars
  use mod_Tchebyshev
  use mod_Matrices
  use mod_PrecompSH

  implicit none

contains

subroutine precompXeYe()

  implicit none

  double precision, dimension(KK) :: Chb0v, ChbD2v

  double precision, dimension(KK2, KK2) :: Xe, Ye ! eq (15)
  double precision :: param
  integer :: k, l, m, j

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! There is one Xe and one Ye matrix for each l. The Xe and Ye matrices are the ones in Hollerbach's paper but multiplied by r
    do k = 1, KK2
      Chb0v = ChbPoly    (k - 1, xK)
      ChbD2v = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
      Xe(:KK, k) = dble(l * (l + 1)) / rK * (param*Chb0v - IER * delta_t * Ek * ( &
                 & ChbD2v - dble(l * (l + 1)) / (rK ** 2) * Chb0v))       
      Ye(:KK, k) = dble(l * (l + 1)) / rK * (param * Chb0v + (1. - IER) * delta_t * Ek * ( &
                 & ChbD2v - dble(l * (l + 1)) / (rK ** 2) * Chb0v))
      !Boundary conditions
      Xe(KK + 1, k) = (-1.) ** (k - 1)
      Xe(KK + 2, k) =   1.
      Ye(KK + 1, k) =   0.
      Ye(KK + 2, k) =   0.
    end do

    Xe_inv  (:, :, l) = invL(Xe)
    Xe_invYe(:, :, l) = transpose(Xe_inv(:, :, l) .dot. Ye)
    Xe_inv  (:, :, l) = transpose(Xe_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXeYe
   
! -------------------------------------------------------------------------

subroutine precompXfYf()

  implicit none

  double precision, dimension(KK) :: Chb0v, ChbD1v, ChbD2v, ChbD4v, fac

  double precision, dimension(KK4, KK4) :: Xf, Yf ! eq (20)
  double precision :: param
  integer :: k, l, m, j, i

  if ((solver == "convective_explicit") .or. (solver == "newton_convective_explicit") & 
    & .or. (solver == "continuation_convective_explicit")) then
    param = Ek
  end if

  do l = 1, LL + 1 ! There is one Xe and one Ye matrix for each l. The Xe and Ye matrices are the ones in Hollerback's paper but multiplied by r**2
    fac = dble(l * (l + 1)) / rK ** 2
    do k = 1, KK4
      Chb0v = ChbPoly    (k - 1, xK)
      ChbD1v = ChbPolDeriv(k - 1, xK, 1) * dxdr
      ChbD2v = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
      ChbD4v = ChbPolDeriv(k - 1, xK, 4) * dxdr ** 4
      Xf(:KK, k) = - dble(l * (l + 1)) * (param * (ChbD2v - fac * Chb0v) &
                 & - IER * delta_t * Ek * (ChbD4v - 2.d0 * fac * ChbD2v &
                 & + 4.d0 * fac / rK * ChbD1v + fac ** 2 * Chb0v &
                 & - 6.d0 * fac / rK ** 2 * Chb0v))
      Yf(:KK, k) = - dble(l * (l + 1)) * (param * (ChbD2v - fac * Chb0v) &
                 & + (1 - IER) * delta_t * Ek * (ChbD4v - 2.d0 * fac * ChbD2v &
                 & + 4.d0 * fac / rK * ChbD1v + fac ** 2 * Chb0v &
                 & - 6.d0 * fac / rK ** 2 * Chb0v))
      !Boundary conditions
      Xf(KK + 1, k) = (-1.) ** (k - 1)
      Xf(KK + 2, k) =   1.
      Xf(KK + 3, k) = (-1.) ** k * (k - 1.) ** 2 ! I think this should be Xf(KK+3,k) = (-1.d0)**(k-1) * (k-1.d0)**2. (JCGS) 26-03-24
      Xf(KK + 4, k) =   1. * (k - 1.) ** 2
      Yf(KK + 1, k) =   0.
      Yf(KK + 2, k) =   0.
      Yf(KK + 3, k) =   0.
      Yf(KK + 4, k) =   0.

    end do
    
    Xf_inv  (:, :, l) = invL(Xf)
    Xf_invYf(:, :, l) = transpose(Xf_inv(:, :, l) .dot. Yf)
    Xf_inv  (:, :, l) = transpose(Xf_inv(:, :, l) .dot. Chb4)

  end do

end subroutine precompXfYf

!-----------------------------------------------------------------------------

subroutine precompXTYT()

  implicit none

  double precision, dimension(KK) :: Chb0v, ChbD1v, ChbD2v

  double precision, dimension(KK2, KK2) :: XTh, YTh ! eq (41)
  integer :: k, l, m, j

  do l = 0, LL + 1 ! We must include the l=m=0 mode - Hollerback's paper. (JCGS) 26-03-24
    do k = 1, KK2
      Chb0v = ChbPoly    (k - 1, xK)
      ChbD1v = ChbPolDeriv(k - 1, xK, 1) * dxdr
      ChbD2v = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
      ! Note that we're using the non-magnetic equations as we use Pr and not the Roberts number as Hollerback does. (JCGS) 26-03-24
      XTh(:KK, k) = Chb0v - IER * delta_t * (ChbD2v + 2.d0 / rK * ChbD1v         &
                  & - l * (l + 1.d0) / rK ** 2 * Chb0v) / Pr
      YTh(:KK, k) = Chb0v + (1.d0 - IER) * delta_t * (ChbD2v + 2.d0 / rK * ChbD1v         &
                  & - l * (l + 1.d0) / rK ** 2 * Chb0v) / Pr
      ! Boundary conditions
      XTh(KK + 1, k) = (-1.) ** (k - 1)
      XTh(KK + 2, k) =   1.
      YTh(KK + 1, k) =   0.
      YTh(KK + 2, k) =   0.
    end do

    XT_inv  (:, :, l) = invL(XTh)
    XT_invYT(:, :, l) = transpose(XT_inv(:, :, l) .dot. YTh)
    XT_inv  (:, :, l) = transpose(XT_inv(:, :, l) .dot. Chb2)

  end do

end subroutine precompXTYT

! ------------------------------------------------------------------------------------------------------
! --------------  Implicit Coriolis --------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------ 


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
! -----------------------------  Creation of blocks --------------------------------------

function Xe(m, l) 

  implicit none
  double precision, dimension(KK)       :: Chb0, Chb2
  integer, intent(in)                   :: l,m
  double precision, dimension(KK2, KK2) :: Xe
  double precision :: param
  integer                               :: k

  Xe = 0.

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  do k = 1, KK2
    Chb0 = ChbPoly    (k - 1, xK)
    Chb2 = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2

    Xe(:KK, k) = dble(l * (l + 1)) / rK * (param * Chb0 - IER * delta_t * Ek * &
               & (Chb2 - dble(l * (l + 1)) / (rK ** 2) * Chb0))

    !Boundary conditions
    Xe(KK + 1, k) = (-1.) ** (k - 1)
    Xe(KK + 2, k) = 1.
  end do  
  
end function Xe

! ---------------------------------------------------------------------------------------

function Xf(m, l) 

  implicit none

  integer, intent(in)                   :: l,m
  double precision, dimension(KK4, KK4) :: Xf
  double precision :: param
  double precision, dimension(KK)       :: Chb0, Chb1, Chb2, Chb4 
  integer                               :: k
  double precision, dimension(KK)       :: fac

  Xf = 0.
  fac = dble(l * (l + 1)) / rK ** 2

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  do k = 1, KK4
    Chb0 = ChbPoly    (k - 1, xK)
    Chb1 = ChbPolDeriv(k - 1, xK, 1) * dxdr
    Chb2 = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
    Chb4 = ChbPolDeriv(k - 1, xK, 4) * dxdr ** 4

    Xf(:KK, k) = - dble(l * (l + 1)) * (param * (Chb2 - fac * Chb0) - IER * delta_t * Ek * (Chb4 - 2. * fac * Chb2 &
               & + 4. * fac / rK * Chb1 + fac ** 2 * Chb0 - 6. * fac / rK ** 2 * Chb0))

    !Boundary conditions
    Xf(KK + 1, k) = (-1.) ** (k - 1)
    Xf(KK + 2, k) = 1.
    Xf(KK + 3, k) = (-1.) ** k * (k - 1.) ** 2
    Xf(KK + 4, k) = 1. * (k - 1.) ** 2
  end do  

end function Xf

! -----------------------------------------------------------------------------------

function Be(m, l, COEF)

  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK2, KK2) :: Be
  double precision, dimension(KK)       :: Chb0
  integer                               :: k
  double precision :: COEF

  ! Testing
  integer :: p

  Be = 0.

  do k = 1, KK2
    Chb0 = ChbPoly(k - 1, xK)

    Be(:KK, k) = - delta_t * COEF * 2 * dble(m) / rK * Chb0

    !Boundary conditions
    Be(KK + 1, k) =  0.
    Be(KK + 2, k) =  0.
  end do

end function Be

! ----------------------------------------------------------------------------------------

function Bf(m, l, COEF)

  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK4, KK4) :: Bf
  double precision, dimension(KK)       :: Chb0, Chb2
  integer                               :: k
  double precision :: COEF

  Bf = 0.
  
  do k = 1, KK4
    Chb0 = ChbPoly    (k - 1, xK)
    Chb2 = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2

    Bf(:KK, k) =  delta_t * COEF * 2. * dble(m) * (Chb2 - dble(l * (l + 1)) / (rK ** 2) * Chb0)
    
    !Boundary conditions
    Bf(KK + 1, k) =  0.d0
    Bf(KK + 2, k) =  0.d0
    Bf(KK + 3, k) =  0.d0
    Bf(KK + 4, k) =  0.d0
  end do    

end function Bf

! ----------------------------------------------------------------------------------------

function Xue(m, l, COEF)
    
  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK4, KK2) :: Xue
  double precision, dimension(KK)       :: Chb0, Chb1, Chb0R
  integer                               :: k
  double precision :: COEF
   
  Xue = 0.

  if ((l>=m) .and. (l + 1 > 1)) then
    do k = 1, KK2
      Chb0  = ChbPoly (k - 1, xK)
      Chb1  = ChbPolDeriv(k - 1, xK, 1) * dxdr
      Chb0R = Chb0 / rK

      Xue(:KK, k) = - delta_t * 2. * COEF * dble(l * (l + abs(m) + 1) * (l + 2)) / dble(2 * l + 3) * &
                  & (dble(l + 1) * Chb0R + Chb1) * SH_norm(m, l+1) / SH_norm(m, l)
      
      !Boundary conditions
      Xue(KK + 1, k) = 0.
      Xue(KK + 2, k) = 0.
      Xue(KK + 3, k) = 0.
      Xue(KK + 4, k) = 0.
    end do
  else
    Xue = 0.
  end if
  
end function Xue

! ----------------------------------------------------------------------------------------

function Xuf(m, l, COEF) 

  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK2, KK4) :: Xuf
  double precision, dimension(KK)       :: Chb0, Chb1, Chb0R
  integer                               :: k
  double precision :: COEF
  
  Xuf = 0.

  if ((l>=m) .and. (l + 1 > 1)) then
    do k = 1, KK4
      Chb0  = ChbPoly    (k - 1, xK)
      Chb1  = ChbPolDeriv(k - 1, xK, 1) * dxdr
      Chb0R = Chb0 / rK

      Xuf(:KK, k) = - delta_t * 2. * COEF / rK * dble(l * (l + abs(m) + 1) * (l + 2)) / dble(2 * l + 3) * &
                  & (dble(l + 1) * Chb0R + Chb1) * SH_norm(m, l+1) / SH_norm(m, l)

      !Boundary conditions
      Xuf(KK + 1, k) = 0.
      Xuf(KK + 2, k) = 0.
    end do    
  else
    Xuf = 0.
  end if
   
end function Xuf

! ----------------------------------------------------------------------------------------

function Xle(m, l, COEF)
        
  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK4, KK2) :: Xle
  double precision, dimension(KK)       :: Chb0, Chb1, Chb0R
  integer                               :: k
  double precision :: COEF

  Xle = 0.

  if (l - 1 < LL) then
    do k = 1, KK2
      Chb0  = ChbPoly (k - 1, xK)
      Chb1  = ChbPolDeriv(k - 1, xK, 1) * dxdr
      Chb0R = Chb0 / rK

      Xle(:KK, k) = delta_t * 2. * COEF * dble((l + 1) * (l - 1) * (l - abs(m))) / dble(2 * l - 1) * &
                  & (dble(l) * Chb0R - Chb1) * SH_norm(m, l - 1) / SH_norm(m, l)

      !Boundary conditions
      Xle(KK + 1, k) = 0.
      Xle(KK + 2, k) = 0.
      Xle(KK + 3, k) = 0.
      Xle(KK + 4, k) = 0.
    end do
  else
    Xle = 0.
  end if
   
end function Xle

! ----------------------------------------------------------------------------------------

function Xlf(m, l, COEF)
     
  implicit none

  integer, intent(in)                  :: l, m
  double precision, dimension(KK2,KK4) :: Xlf
  double precision, dimension(KK)      :: Chb0, Chb1, Chb0R
  integer                              :: k
  double precision :: COEF

  Xlf = 0.

  if (l - 1 < LL) then
    do k = 1, KK4
      Chb0  = ChbPoly    (k - 1, xK)
      Chb1  = ChbPolDeriv(k - 1, xK, 1) * dxdr
      Chb0R = Chb0 / rK

      Xlf(:KK, k) = delta_t * 2. * COEF / rK * dble((l + 1) * (l - 1) * (l - abs(m))) / dble(2 * l - 1) * &
                  & (dble(l) * Chb0R - Chb1) * SH_norm(m, l-1) / SH_norm(m, l)

      !Boundary conditions
      Xlf(KK + 1, k) = 0.
      Xlf(KK + 2, k) = 0.
    end do
  else
    Xlf = 0.
  end if
 
end function Xlf

! ----------------------------------------------------------------------------------------

function Ye(m, l) 

  implicit none

  double precision, dimension(KK)       :: Chb0, Chb2
  integer, intent(in)                   :: l, m
  double precision, dimension(KK2, KK2) :: Ye
  double precision :: param
  integer                               :: k

  Ye = 0.

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  do k = 1, KK2
    Chb0 = ChbPoly    (k - 1, xK)
    Chb2 = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2

    Ye(:KK, k) = dble(l * (l + 1)) / rK * (param * Chb0 + (1- IER) * delta_t * Ek * ( &
               & Chb2 - dble(l * (l + 1)) / (rK ** 2) * Chb0))

    !Boundary conditions
    Ye(KK + 1, k) = 0.
    Ye(KK + 2, k) = 0.
  end do

end function Ye

! ----------------------------------------------------------------------------------------

function Yf(m, l) 

  implicit none

  integer, intent(in)                   :: l, m
  double precision, dimension(KK4, KK4) :: Yf
  double precision :: param
  double precision, dimension(KK)       :: Chb0, Chb1, Chb2, Chb4 
  integer                               :: k
  double precision, dimension(KK)       :: fac

  Yf = 0.

  fac = dble(l * (l + 1)) / rK ** 2

  if ((solver == "convective_implicit") .or. (solver == "newton_convective_implicit") & 
    & .or. (solver == "continuation_convective_implicit")) then
    param = Ek
  end if

  do k = 1, KK4
    Chb0 = ChbPoly    (k - 1, xK)
    Chb1 = ChbPolDeriv(k - 1, xK, 1) * dxdr
    Chb2 = ChbPolDeriv(k - 1, xK, 2) * dxdr ** 2
    Chb4 = ChbPolDeriv(k - 1, xK, 4) * dxdr ** 4

    Yf(:KK, k) = - dble(l * (l + 1)) * (param * (Chb2 - fac * Chb0) + & 
               & (1 - IER) * delta_t * Ek * (Chb4 - 2. * fac * Chb2 + &
               & 4. * fac / rK * Chb1 + fac ** 2 * Chb0 - 6. * fac / rK ** 2 * Chb0))

    !Boundary conditions
    Yf(KK + 1, k) = 0.
    Yf(KK + 2, k) = 0.
    Yf(KK + 3, k) = 0.
    Yf(KK + 4, k) = 0.
  end do
   
end function Yf

! ----------------------------------------------------------------------------------------
! -----------------------------  Creation of columns -------------------------------------

function CXF(l, m)

  implicit none
  double precision, dimension(4 * KK2 + 2 * KK4, 2 * KK4)::CXF
  integer :: step, l, m

  step = 2 * (KK2 + KK4)

  CXF=0.

  CXF(1 : KK2, 1 : KK4) = Xuf(m, l - 1, IER)                                           
  CXF(KK2 + 1 : 2 * KK2, KK4 + 1 : 2 * KK4) = Xuf(m, l - 1, IER)                               
  CXF(2 * KK2 + 1 : 2 * KK2 + KK4, 1 : KK4) = Xf(m, l)                                    
  CXF(2 * KK2 + 1 : 2 * KK2 + KK4, KK4 + 1 : 2 * KK4) = - Bf(m, l, IER)                          
  CXF(2 * KK2 + KK4 + 1 : step, 1 : KK4) = Bf(m, l, IER)                                   
  CXF(2 * KK2 + KK4 + 1 : step, KK4 + 1 : 2 * KK4) = Xf(m, l)                                  
  CXF(1 + step : step + KK2, 1 : KK4) = Xlf(m, l + 1, IER)                                       
  CXF(step + KK2 + 1 : step + 2 * KK2, KK4 + 1 : 2 * KK4) = Xlf(m, l + 1, IER)  

end function CXF

! ----------------------------------------------------------------------------------------
   
function CYF(l,m)

  implicit none

  double precision, dimension(4 * KK2 + 2 * KK4, 2 * KK4) :: CYF
  integer :: step, l, m

  step = 2 * (KK2 + KK4)

  CYF = 0.

  CYF(1 : KK2, 1 : KK4) = - Xuf(m, l - 1, (1 - IER))
  CYF(KK2 + 1 : 2 * KK2, KK4 + 1 : 2 * KK4) = - Xuf(m, l - 1, (1 - IER))
  CYF(2 * KK2 + 1 : 2 * KK2 + KK4, 1 : KK4) = Yf(m, l)
  CYF(2 * KK2 + 1 : 2 * KK2 + KK4, KK4 + 1 : 2 * KK4) = Bf(m, l, (1 - IER))
  CYF(2 * KK2 + KK4 + 1 : step, 1 : KK4) = - Bf(m, l, (1 - IER))
  CYF(2 * KK2 + KK4 + 1 : step, KK4 + 1 : 2 * KK4) = Yf(m, l)
  CYF(1 + step : step + KK2, 1 : KK4) = - Xlf(m, l + 1, (1 - IER))
  CYF(step + KK2 + 1 : step + 2 * KK2, KK4 + 1 : 2 * KK4) = - Xlf(m, l + 1, (1 - IER))

end function CYF

! -----------------------------------------------------------------------------------------

function CXE(l, m)

  implicit none
  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: CXE
  integer :: step, l, m, i, j

  step = 2 * (KK2 + KK4)

  CXE = 0.

  CXE(1 : KK4, 1 : KK2) = Xue(m, l - 1, IER)                                           
  CXE(KK4 + 1 : 2 * KK4, KK2 + 1 : 2 * KK2) = Xue(m, l - 1, IER)                               
  CXE(2 * KK4 + 1 : 2 * KK4 + KK2, 1 : KK2) = Xe(m, l)                                    
  CXE(2 * KK4 + 1 : 2 * KK4 + KK2, KK2 + 1 : 2 * KK2) = - Be(m, l, IER)                          
  CXE(2 * KK4 + KK2 + 1 : step, 1 : KK2)= Be(m, l, IER)                                   
  CXE(2 * KK4 + KK2 + 1 : step, KK2 + 1 : 2 * KK2) = Xe(m, l)                                  
  CXE(1 + step : step + KK4, 1 : KK2) = Xle(m, l + 1, IER)                                       
  CXE(step + KK4 + 1 : step + 2 * KK4, KK2 + 1 : 2 * KK2) = Xle(m, l + 1, IER)

end function CXE

! ----------------------------------------------------------------------------------------

function CYE(l,m)

  implicit none

  double precision, dimension(2 * KK2 + 4 * KK4, 2 * KK2) :: CYE
  integer :: step, l, m

  step = 2 * (KK2 + KK4)

  CYE = 0.

  CYE(1 : KK4, 1 : KK2) = - Xue(m, l - 1, (1 - IER))
  CYE(KK4 + 1 : 2 * KK4, KK2 + 1 : 2 * KK2) = - Xue(m, l - 1, (1 - IER))
  CYE(2 * KK4 + 1 : 2 * KK4 + KK2, 1 : KK2) = Ye(m, l)
  CYE(2 * KK4 + 1 : 2 * KK4 + KK2, KK2 + 1 : 2 * KK2) = Be(m, l, (1 - IER))
  CYE(2 * KK4 + KK2 + 1 : step , 1 : KK2) = - Be(m, l, (1 - IER))
  CYE(2 * KK4 + KK2 + 1 : step, KK2 + 1 : 2 * KK2) = Ye(m, l)
  CYE(1 + step : step + KK4, 1 : KK2) = - Xle(m, l + 1, (1 - IER))
  CYE(step + KK4 + 1 : step + 2 * KK4, KK2 + 1 : 2 * KK2) = - Xle(m, l + 1, (1 - IER))

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
      ColumnEX = CXE(l, m)
      ColumnFX = CXF(l, m)
      ColumnEY = CYE(l, m)
      ColumnFY = CYF(l, m)

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
      ColumnFX = CXF(l, m)
      ColumnEX = CXE(l, m)
      ColumnFY = CYF(l, m)
      ColumnEY = CYE(l, m)

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

end module mod_precompXY
