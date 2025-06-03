program iterative_solver_test
  
  use mod_Globalvars
  use mod_CommandLineParser
  use mod_PrecompSH
  use mod_Matrices
  use mod_Tchebyshev
  use mod_precompXY
  use mod_read
  use mod_init
  use mod_TimeSteppingSolver
  use mod_NewtonSolver
  use mod_TimeStep
  use mod_IterativeSolvers_test
  use mod_Output

  implicit none


  ! definition of program inputs and parameters
  integer, parameter :: m = 110                    ! the size of the linear system
  character(len = 100) :: coeff_dir = './A10.txt'  ! coefficient matrix directory (column-major order)
  character(len = 100) :: rhs_dir = './b10.txt'    ! RHS vector directory
  double precision :: res_tol = 1.0d-6             ! solver convergence threshold
  double precision, dimension(m) :: x0             ! initial condition
  
  ! declaring program outputs
  double precision, dimension(m) :: ubest          ! iterative solver best approximation
  integer :: iters                                 ! number of solver itertations
  
  ! internal variables
  double precision, dimension(m, m) :: A_          ! the coefficient matrix
  double precision, dimension(m) :: b_             ! the RHS vector
  integer :: unit = 0                              ! unit number for file input and output
  integer :: ios                                   ! input/output status
  integer :: k                                     ! loop counter

  
  print*, "################################################################"
  print*, "#                                                              #"
  print*, "#             Testing different iterative solvers              #"
  print*, "#                                                              #"
  print*, "# Solver 1: GMRES                                              #"
  print*, "# Solver 2: IDR(s)                                             #"
  print*, "# Solver 3: BiCGSTAB                                           #"
  print*, "#                                                              #"
  print*, "################################################################"


  ! ---------- read coefficient matrix from file
  open(unit, file = coeff_dir, status = 'old', action = 'read', iostat = ios)
    
  if (ios /= 0) then
    print *, "Error opening file."
    stop
  end if

  read(unit, *) A_
  close(unit)
  unit = unit + 1
  print*, "...Coefficient matrix successfully read from file."

  
  ! ---------- read RHS vector from file
  open(unit, file = rhs_dir, status = 'old', action = 'read', iostat = ios)
    
  if (ios /= 0) then
    print *, "Error opening file."
    stop
  end if

  read(unit, *) b_
  close(unit)
  unit = unit + 1
  print*, "...RHS vector successfully read from file."
  
  
  ! ---------- test GMRES
  print*, "                                                                "
  print*, "############################ GMRES #############################"
  
  call GMRES(A_, m, m, m, 1.0d-6, b_, ubest, iters)
  ! call print_vector(ubest, m, 10)
  
  ! ---------- test IDR(s)
  print*, "                                                                "
  print*, "############################ IDR(s) ############################"
  
  x0 = 0
  call IDRs(A_, b_, 8, m, 1.0d-6, m, x0, ubest, iters)
  ! call print_vector(ubest, m, 10)
  
  ! ---------- test BiCGSTAB
  print*, "                                                                "
  print*, "########################### BiCGSTAB ###########################"
  
  x0 = 0
  call BiCGSTAB(A_, b_, m, 1.0d-6, m, x0, ubest, iters)
  ! call print_vector(ubest, m, 10)
  
end program iterative_solver_test
