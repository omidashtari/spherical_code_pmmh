module mod_Tchebyshev

  use mod_Globalvars
  use mod_Matrices

  implicit none

  double precision :: Vol
  double precision :: dxdr

  double precision, dimension(:), allocatable :: xN, xK  !--- Roots of the kN-th Tchebyshev polynomial
  double precision, dimension(:), allocatable :: rN, rK

  double precision, dimension(:, :), allocatable :: Chb      !--- Tchebyshev polynomials, first coodinate is the order and the second the values
  double precision, dimension(:, :), allocatable :: ChbR     !--- Chebyshev/R
  double precision, dimension(:, :), allocatable :: ChbR2    !--- Chebyshev/R2
  double precision, dimension(:, :), allocatable :: ChbR3    !--- Chebyshev/R3
  double precision, dimension(:, :), allocatable :: ChbD1    !--- dChebyshev/dr
  double precision, dimension(:, :), allocatable :: ChbD1R   !--- dChebyshev/dr * 1/R

  double precision, dimension(:, :), allocatable :: Chbderiv, Chbderiv2 ! To compute first derivatives in Chebyshev space

  double precision, dimension(:, :), allocatable :: Chb_mulR_deriv ! To compute d(r())/dr

  double precision, dimension(:,:), allocatable :: Chbinv !--- inverse Chebyshev transform

  double precision, dimension(:, :), allocatable :: Chb2, Chb4 ! To compute transform with BC for the solve

contains

  subroutine PrecompTchebyshev()

    implicit none

    integer :: k, p

    ! Initialize variables
    kN = max(int(3 * KK / 2) + 1, KK + 4)
    KK2 = KK+2
    KK4 = KK+4

    ! Define dxdr
    dxdr = 2.d0 / (Rout-Rin) ! This is because r = (Rin + Rout)/2.d0 + (Rout - Rin)/2.d0 * x

    !Allocate arrays
    allocate(xN(kN), rN(kN))
    allocate(xK(KK), rK(KK))

    allocate(Chb(kN,kN))
    allocate(ChbR(kN,kN))
    allocate(ChbR2(kN,kN))
    allocate(ChbR3(kN,kN))
    allocate(ChbD1(kN,kN))
    allocate(ChbD1R(kN,kN))

    allocate(Chbderiv(0:KK4, 0:KK4), Chbderiv2(0:KK4, 0:KK4))

    allocate(Chbinv(kN,kN))

    allocate(Chb2(KK2, KK2), Chb4(KK4, KK4))

    ! First of all we compute the volume
    Vol = 4./3. * pi * (Rout**3 - Rin**3)

    !--- Collocation points for Chebyshev polynomials
    ! The roots of the T_n Chebyshev polynomial are defined as:
    ! x_k = cos( pi (k + 1/2 ) / n ), k = 0, ..., n-1
    ! We add a - in front of the expression for the x_k to be ordered from -1 to 1
    xN = -cos( (/(2.d0*k-1.d0, k = 1,kN )/)/kN * pi/2.d0 )
    xK = -cos( (/(2.d0*k-1.d0, k = 1,KK )/)/KK * pi/2.d0 )

    !--- Corresponding points in R
    ! We compute the r discretization as:
    ! r = mid point between Rin and Rout + difference between Rin and Rout over 2 * x_k 
    ! We recall that x_k belongs to the interval [-1, 1]
    rN = (Rin + Rout)/2.d0 + (Rout - Rin)/2.d0 * xN
    rK = (Rin + Rout)/2.d0 + (Rout - Rin)/2.d0 * xK
    
    do k = 1, kN
      Chb   (k, :) = ChbPoly(k - 1, xN) ! In ChbPoly, the size of C is its first parameter + 1, so we get k-1+1 = k.
      ChbR  (k, :) = Chb  (k, :) / rN
      ChbR2 (k, :) = ChbR (k, :) / rN
      ChbR3 (k, :) = ChbR2(k, :) / rN
      ChbD1 (k, :) = ChbPolDeriv(k - 1, xN, 1) * dxdr
      ChbD1R(k, :) = ChbPolDeriv(k - 1, xN, 1) / rN * dxdr
    end do

    ! Now for the derivatives in Chebyshev space
    Chbderiv = 0.
    do k = 0, KK4-1
      do p = k+1, KK4, 2
        Chbderiv(k, p) = 2.0 * p * dxdr
      end do
    end do
    Chbderiv(0, :) =  Chbderiv(0, :) / 2.

    Chbderiv2 = 0.
    do k = 0, KK4-2
      do p = k+2, KK4, 2
          Chbderiv2(k, p) = p * (p**2 - k**2) * dxdr**2
      end do
    end do
    Chbderiv2(0, :) =  Chbderiv2(0, :) / 2.

    ! Create matrix for multiplication by d(r())/dr
    allocate(Chb_mulR_deriv(KK4, KK4))
    Chb_mulR_deriv = 0.
    do k = 2, KK4-1
      Chb_mulR_deriv(k, k) = (Rout + Rin) / 2.
      Chb_mulR_deriv(k, k + 1) = (Rout - Rin) / 4.
      Chb_mulR_deriv(k, k - 1) = (Rout - Rin) / 4.
    end do

    Chb_mulR_deriv(1, 1) = (Rout + Rin) / 2.
    Chb_mulR_deriv(KK4, KK4) = (Rout + Rin) / 2.
    Chb_mulR_deriv(1, 2) = (Rout - Rin) / 4.
    Chb_mulR_deriv(2, 1) = Chb_mulR_deriv(2, 1) + (Rout - Rin) / 4.
    Chb_mulR_deriv(KK4, KK4-1) = (Rout - Rin) / 4.

    Chb_mulR_deriv = Chbderiv(0:KK4-1, 0:KK4-1) .dot. Chb_mulR_deriv(:KK4, :KK4)

    ! And lastly we get the inverse
    Chbinv = invL(Chb)

    ! The only thing missing are the matrices for the transform with BC
    Chb2 = 0.
    Chb4 = 0.

    do k = 1, KK
      Chb2(:KK, k) = ChbPoly(k - 1, xK) * delta_t
      Chb4(:KK, k) = ChbPoly(k - 1, xK) * delta_t
    end do
    Chb2(KK + 1, KK + 1) = 1.
    Chb2(KK + 2, KK + 2) = 1.
    Chb4(KK + 1, KK + 1) = 1.
    Chb4(KK + 2, KK + 2) = 1.
    Chb4(KK + 3, KK + 3) = 1.
    Chb4(KK + 4, KK + 4) = 1.

  end subroutine PrecompTchebyshev

  !############################################################################
  !
  !                           Chebyshev Polynomials
  !
  !############################################################################

  function ChbPolDeriv(k,x,m) result(C)

    !**************************************************************************
    ! Compute the m-th derivative of the k-th Chebyshev polynomial at x.
    ! Need:
    !   * m (INTEGER), the order of the derivative (can be 0)
    !   * k (INTEGER), the order of the Chebyshev polynomial
    !   * x (DOUBLE PRECISION, DIMENSION) contains the points at which the
    !                   derivative must be evaluated
    !**************************************************************************

    implicit none

    integer, intent(in)                          :: m,k
    double precision, dimension(:), intent(in)   :: x
    double precision, dimension(size(x))         :: C ! Here we store the point values of the derivative of the kth polynomial. (JCGS) 25-03-24
    double precision, dimension(k+1)             :: A
    !double precision, dimension(0:10000,size(x)) :: B
    integer                                      :: i,j,N

    A(:) = 0.d0
    A(k+1) = 1.d0

    call DiffChbArray(A,m)

    C = ChbPolyTab(A,x)

  end function ChbPolDeriv

  !----------------------------------------------------------------------------

  subroutine DiffChbArray(A,m)

    !**************************************************************************
    ! Compute the derivative of a function f(x) of the form
    !          f(x) = [A0 A1 A2 ... An]^t * [T_0 T_1 T_2 ... T_n](x)
    ! i.e. returns the coefficients of f'(x) in the Chebyshev basis
    !**************************************************************************

    implicit none

    double precision, dimension(:), intent(inout) :: A
    integer, intent(in)                           :: m
    double precision, dimension(size(A,1))        :: der
    integer                                       :: i,j,n

    n = size(A,1)

    if (m >= n) then ! i.e. if the degree of the derivative is bigger or equal than the order of the polynomial, A = 0. (JCGS) 25-03-24
      A = 0
    else if (m > 0) then
      do i = 1,m
        n = n-1 ! n goes down as i increases. (JCGS) 25-03-24
        do j = n,3,-1 ! j takes on integer values starting from n and decrements by 1 until it reaches 3. (JCGS) 25-03-24
           der(j) = 2.d0 * j * A(j+1)
           A(j-1) = A(j-1) + j * A(j+1) / (j-2.d0)
        end do
        if (n > 1) der(2) = 4*A(3)
        der(1) = A(2)
        A(:) = der(:)
        A(n+1:) = 0.d0
      end do
    end if

  end subroutine DiffChbArray

  !----------------------------------------------------------------------------

  function ChbPoly(k,x) result(Chb)

    !**************************************************************************
    ! Compute the k-th Chebyshev Polynomial at points contained in x
    !   * k (INTEGER) is the order of the polynomial
    !   * x (DOUBLE PRECISION, DIMENSION) contains the points
    ! Note that we only use Chebyshev polynomials of the first kind
    !**************************************************************************

    implicit none

    integer, intent(in)                        :: k
    double precision, dimension(:), intent(in) :: x
    double precision, dimension(k+1)           :: C 
    double precision, dimension(size(x))       :: Chb

    C = 0.d0 ! We initialize the C vector as zeros.
    C(k+1) = 1.d0 ! We set the last element to 1

    Chb = ChbPolyTab(C,x)

  end function ChbPoly

  !----------------------------------------------------------------------------

  function ChbPolyTab(C,x) result(Chb)

    !**************************************************************************
    ! Compute the k-th Chebyshev Polynomial at points contained in x
    !   * k (INTEGER) is the order of the polynomial
    !   * x (DOUBLE PRECISION, DIMENSION) contains the points
    ! Note that we only use Chebyshev polynomials of the first kind
    !**************************************************************************

    implicit none

    double precision, dimension(:), intent(in) :: C
    double precision, dimension(:), intent(in) :: x
    double precision, dimension(size(x))       :: Chb
    double precision, dimension(size(x))       :: buffer1, buffer2
    integer                                    :: i

    if (size(C,1) == 1) then ! Here doing size(C,1) is the same as size(C) given that C is a one dimensional array. (JCGS) 25-03-24
      buffer1 = C(1)
      buffer2 = 0.d0
    else if (size(C) == 2) then
      buffer1 = C(1)
      buffer2 = C(2)
    else
      buffer1 = C(size(C,1)-1)
      buffer2 = C(size(C,1))
      do i = 3, size(C,1)
        Chb = buffer1
        buffer1 = C(size(C,1)-i+1) - buffer2
        buffer2 = Chb + buffer2 * 2.d0 * x
      end do
    end if
    Chb = buffer2*x + buffer1

  end function ChbPolyTab

end module mod_Tchebyshev
