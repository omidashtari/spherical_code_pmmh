# parameters
DEBUG= -O2 -fdefault-real-8 -fdefault-double-8 -g -fcheck=all # -fopenmp -openmp-report2 -Wall
#DEBUG=  -fopenmp
COMPILE=gfortran $(DEBUG)

# Path to modules and datas
SRC = src
OBJ = obj
MOD = src/mods
SHTNS_DIR = /home/juan/shtns

# path to LAPACK and BLAS libraries
LAPACK_PATH = /usr/local/lib
# name of the BLAS and LAPACK libraries
BLAS_NAME = blas
LAPACK_NAME = lapack
SHTNS_NAME = shtns
FFTW3_NAME = fftw3
# Flags for including LAPACK and BLAS
LAPACK_FLAGS = -L$(LAPACK_PATH) -l$(LAPACK_NAME) -l$(BLAS_NAME) -l$(SHTNS_NAME) -l$(FFTW3_NAME) -lm -lc

# final compilation
MainMHD : MainMHD.o mod_GlobalVars.o mod_CommandLineParser.o mod_PrecompSH.o mod_Matrices.o mod_Tchebyshev.o mod_PrecompXY.o mod_tools.o mod_read.o mod_ExplicitTerms.o mod_init.o mod_TimeStep.o mod_Output.o mod_TimeSteppingSolver.o mod_NewtonSolver.o mod_ActNewtonSolver.o mod_IterativeSolvers.o
	$(COMPILE) -o MainMHD $(OBJ)/*.o $(LAPACK_FLAGS)
	@echo compilation done

# object compilation
MainMHD.o : $(SRC)/MainMHD.f90 mod_GlobalVars.o mod_CommandLineParser.o mod_PrecompSH.o mod_Matrices.o mod_Tchebyshev.o mod_PrecompXY.o mod_tools.o mod_read.o mod_ExplicitTerms.o mod_init.o mod_TimeStep.o mod_Output.o mod_TimeSteppingSolver.o mod_NewtonSolver.o mod_ActNewtonSolver.o mod_IterativeSolvers.o
	$(COMPILE) -I$(OBJ) -c $(SRC)/MainMHD.f90 -o$(OBJ)/MainMHD.o
mod_Output.o : $(MOD)/mod_Output.f90 mod_GlobalVars.o mod_PrecompSH.o mod_Matrices.o  mod_Tchebyshev.o mod_ExplicitTerms.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_Output.f90 -o$(OBJ)/mod_Output.o
mod_TimeStep.o : $(MOD)/mod_TimeStep.f90 mod_GlobalVars.o mod_Matrices.o mod_PrecompXY.o mod_ExplicitTerms.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_TimeStep.f90 -o$(OBJ)/mod_TimeStep.o
mod_TimeSteppingSolver.o : $(MOD)/mod_TimeSteppingSolver.f90 mod_GlobalVars.o mod_Matrices.o mod_PrecompSH.o mod_Tchebyshev.o mod_PrecompXY.o mod_read.o mod_init.o mod_Output.o mod_TimeStep.o mod_IterativeSolvers.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_TimeSteppingSolver.f90 -o$(OBJ)/mod_TimeSteppingSolver.o
mod_IterativeSolvers.o : $(MOD)/mod_IterativeSolvers.f90 mod_GlobalVars.o mod_Matrices.o mod_ActNewtonSolver.o
	$(COMPILE) -I$(OBJ) -I$(SHTNS_DIR) -J$(OBJ) -c $(MOD)/mod_IterativeSolvers.f90 -o$(OBJ)/mod_IterativeSolvers.o $(LAPACK_FLAGS)
mod_NewtonSolver.o : $(MOD)/mod_NewtonSolver.f90 mod_GlobalVars.o mod_PrecompSH.o mod_read.o mod_TimeStep.o mod_ActNewtonSolver.o mod_IterativeSolvers.o mod_Output.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_NewtonSolver.f90 -o$(OBJ)/mod_NewtonSolver.o
mod_ActNewtonSolver.o : $(MOD)/mod_ActNewtonSolver.f90 mod_GlobalVars.o mod_Matrices.o mod_PrecompSH.o mod_TimeStep.o mod_ExplicitTerms.o
	$(COMPILE) -I$(OBJ) -I$(SHTNS_DIR) -J$(OBJ) -c $(MOD)/mod_ActNewtonSolver.f90 -o$(OBJ)/mod_ActNewtonSolver.o $(LAPACK_FLAGS)
mod_init.o : $(MOD)/mod_init.f90 mod_GlobalVars.o mod_PrecompSH.o mod_Matrices.o  mod_Tchebyshev.o mod_ExplicitTerms.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_init.f90 -o$(OBJ)/mod_init.o
mod_ExplicitTerms.o : $(MOD)/mod_ExplicitTerms.f90 mod_GlobalVars.o mod_Matrices.o mod_PrecompSH.o mod_Tchebyshev.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_ExplicitTerms.f90 -o$(OBJ)/mod_ExplicitTerms.o
mod_read.o : $(MOD)/mod_read.f90 mod_GlobalVars.o mod_PrecompSH.o mod_Matrices.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_read.f90 -o$(OBJ)/mod_read.o
mod_tools.o : $(MOD)/mod_tools.f90 mod_GlobalVars.o mod_Tchebyshev.o mod_PrecompSH.o mod_PrecompXY.o mod_Output.o mod_read.o 
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_tools.f90 -o$(OBJ)/mod_tools.o
mod_PrecompXY.o : $(MOD)/mod_PrecompXY.f90 mod_GlobalVars.o mod_Matrices.o mod_Tchebyshev.o mod_PrecompSH.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_PrecompXY.f90 -o$(OBJ)/mod_PrecompXY.o
mod_Tchebyshev.o : $(MOD)/mod_Tchebyshev.f90 mod_GlobalVars.o mod_Matrices.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_Tchebyshev.f90 -o$(OBJ)/mod_Tchebyshev.o
mod_Matrices.o : $(MOD)/mod_Matrices.f90 mod_GlobalVars.o mod_PrecompSH.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_Matrices.f90 -o$(OBJ)/mod_Matrices.o $(LAPACK_FLAGS)
mod_PrecompSH.o : $(MOD)/mod_PrecompSH.f90 mod_GlobalVars.o
	$(COMPILE) -I$(OBJ) -I$(SHTNS_DIR) -J$(OBJ) -c $(MOD)/mod_PrecompSH.f90 -o$(OBJ)/mod_PrecompSH.o $(LAPACK_FLAGS)
mod_CommandLineParser.o : $(MOD)/mod_CommandLineParser.f90 mod_GlobalVars.o
	$(COMPILE) -I$(OBJ) -J$(OBJ) -c $(MOD)/mod_CommandLineParser.f90 -o$(OBJ)/mod_CommandLineParser.o $(LAPACK_FLAGS)
mod_GlobalVars.o : $(MOD)/mod_GlobalVars.f90
	$(COMPILE) -J$(OBJ) -c $(MOD)/mod_GlobalVars.f90 -o$(OBJ)/mod_GlobalVars.o

# cleaning after compilation
clean :
	rm -f obj/*.o obj/*.mod
