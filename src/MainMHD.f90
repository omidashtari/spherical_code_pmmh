program MainMHD
  !$ use OMP_LIB
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
  use mod_IterativeSolvers
  use mod_Output

  implicit none


  print*, "################################################################"
  print*, "#                                                              #"
  print*, "#                   MHD in Spherical Geometry                  #"
  print*, "#            This version uses SHTns for SH transforms         #"
  print*, "#                           v2.1, 2024                         #"
  print*, "#                                                              #"
  print*, "################################################################"
  print*
  print*, "--------------- Reading the Simulation parameters --------------"
  print*, "#                                                              #"
  call parse_command_line_arguments()
  print*, "#                                                              #"
  print*, "------------------- Precomputing the matrices ------------------"
  print*
  print*, "Precomputing the Chebyshev terms..."
  call PrecompTchebyshev()
  print*, "Precomputing the spherical harmonics grid..."
  call PrecompSH()
  call output_coordinates()

  if ((solver == "convective_explicit") .or. (solver == "convective_implicit")) then
    ! We call the covective solver with explicit Coriolis
    call convective_solver()
  else if ((solver == "newton_convective_explicit") .or. (solver == "newton_convective_implicit")) then
    ! We call the Newton covective solver
    call convective_newton_solver()
  else if ((solver == "continuation_convective_explicit") .or. (solver == "continuation_convective_implicit")) then
    ! We call the continuation covective solver
    call continuation_convective_solver()
  end if

end program MainMHD
