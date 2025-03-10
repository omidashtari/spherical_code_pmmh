---
    geometry: margin=1in
---

# Code options

\hspace{1cm} The objective of this document is to outline the different options that we can use when running the code. We can perform time-stepping simulations or use the Newton or the Continuation solvers.

* Time-stepping simulations:

    * delta_t (float): time step for the simulation.
    * NTS (integer): total number of iterations.
    * save_every (integer): Files (output and restart files) will be saved every *save_every* iterations.
    * solver (string): for these simulations it should be either *convective_explicit* or *convective_implicit* for explicit or implicit treatment of the Coriolis term.
    * time_step (string): time stepping scheme for the simulation. It may be *pc, cn, fbe* or *bdf2* for Predictor-Corrector, Crank-Nicolson, Forward-Backward Euler or BDF2.
    * restart (string): it can be *yes, y, no, n*. With this we will indicate the code to begin the simulation by reading the fields from binary files.
    * restart_filename (string): file name for restart file containing fields. If no name is specified the code defaults to *Restart.b* as filename.
    * dim_filename (string): file name for restart file containing dimensions. If no name is specified the code defaults to *Dim.b* as filename.
    * dealiasing (string): it can be *yes, y, no, n*. With this we will indicate the code to use dealiasing techniques for the transformations to real space and back.
    * directory (string): indicate the directory where the simulation will be executed. If doing a restart from binary files, these have to be in that directory.
    * KK (integer): number of Chebyshev modes.
    * LL (integer): number of Legendre modes.
    * MM (integer): number of Fourier modes.
    * mres (integer): $\phi$ will span an angle range between 0 and $2 \pi/$mres.
    * Rin and Rout (floats): inner and outer radii. Default values are set to 7/13 and 20/13, respectively.
    * Pr (float): Prandtl number.
    * Ek (float): Ekman number.
    * Ra (float): Rayleigh number.
    * IER (float): implicit to explicit ratio. Only available in Predictor-Corrector and Crank-Nicolson time-stepping schemes. For more details see R. Hollerbach *A spectral solution of the magneto-convection equations in spherical geometry* in Int. J. Numer. Meth. Fluids 2000.
    * init (string): initial condition. Can be set to *Christensen* for the U.R. Christensen *et al.* *A numerical dynamo benchmark*, in Physics of the Earth and Planetary Interiors 128 (2001) initial condition. Otherwise, it can be set to *symmetric*, for an M-fold symmetry initial condition. Should this last option be chosen, two other parameters are to be set:
        * sym (integer): M-fold symmetry of choice.
        * init_amp (float): amplitude of the initial condition.
    The resulting initial condition will be ${\rm init\_ amp} \cdot \cos ({\rm sym} \cdot \phi)$
    
    An example of a command line to begin a simulation from a symmetric initial condition would be:

    ./MainMHD -delta_t 1.0e-4 -NTS 5000 -save_every 5000 -KK 30 -LL 40 -MM 10 -mres 4 -Pr 1. -Ek 1.0e-3 -Ra 65. -time_step bdf2 -directory 'my_directory' -restart no -init symmetric -sym 4 -init_amp 4 -dealiasing yes -solver convective_explicit

* Newton Solver:

    All the options stated for the time stepping simulations are valid except for: NTS, save_every, init, init_amp and sym. The restart option must be set to *yes* and filenames should be specified unless they go by the default names indicated above.

    * delta_t (float): time step for the simulation. Should be set delta_t $\approx 10^2$.
    * solver (string): for these simulations it should be either *newton_convective_explicit* or *newton_convective_implicit* for explicit or implicit treatment of the Coriolis term.
    * time_step (string): time stepping scheme for these simulations may be *cn* or *fbe* for Crank-Nicolson or Forward-Backward Euler. The latter is the default one.
    * max_newt (integer): maximum amount of Newton iterations allowed. Usual values: 10-15.
    * max_gmres (integer): maximum amount of GMRES iterations allowed in each Newton step. Usual values: 500-1000.
    * restart_gmres (integer): in case GMRES restart is desired by the user. Otherwise we recommend to set it equal to max_gmres for no restart.
    * newt_eps (float): tolerance for the Newton method. Usual values: $1\times 10^{-7}$
    * newt_delta (float): stagnation criterion to evaluate how much the Newton solution is changing. Usual values: $1\times 10^{-16}$.
    * tol_gmres (float): GMRES tolerance. Usual values: $1\times 10^{-10}$.
    * M_wave (integer): symmetry of the rotating wave to be found.

    An example of a command line would be:

    ./MainMHD -delta_t 200. -KK 30 -LL 40 -MM 10 -mres 4 -Pr 1. -Ek 1.0e-3 -Ra 145. -IER 0.8 -directory 'my_directory' -restart yes -restart_filename Restart_Ra_140.b -dim_filename Dim_Ra_140.b -dealiasing yes -max_newt 15 -max_gmres 1000 -restart_gmres 1000 -newt_eps 1.0e-7 -newt_delta 1.0e-16 -tol_gmres 1.0e-10 -M_wave 4 -time_step fbe -solver newton_convective_implicit

* Continuation Solver

    All the options stated for the Newton Solver are valid and mandatory. The following are more options that need to be specified for the continuation solver to work:

    * Ek_final (float): final for continuation in Ekman.
    * Ra_final (float): final for continuation in Rayleigh. 
    * delta_param (float): to indicate how much the parameter will change in each continuation step. Rayleigh changes in linear scale and Ekman in logarithmic scale. Usual values for Rayleigh are 1 to 5 and in Ekman around $1\times 10^{-2}$.
    * adapt_param (string): it can be *yes, y, no, n*. The variation in the parameter will be adapted according to the optimal amount of Newton steps set by the user (Nopt). Default is *no*.
    * Nopt (integer): if adapt_param is set to *yes*, the variation in the parameter will be corrected according to equation (16) from K. Borońska & L. S. Tuckeman *Extreme multiplicity in cylindrical Rayleigh-Bénard convection. II. Bifurcation diagram and symmetry classification*, in Physical Review E (2010). Usual values for Nopt are 3 to 5.
    * grid_refine (string): it can be *yes, y, no, n*. It will indicate the code to perform grid refinement when the spectral resolution is below a threshold. Default is *no*.
    * gr_threshold (float): in case grid_refine is set to *yes* the grid will be refined if the ratio between the last mode and the mode with maximum absolute value is below this threshold. Usual values are $1\times 10^{-7}$.

    Note: one can only set Ek_final *or* Ra_final, not both. 
    
    An example of a command line would be:

    ./MainMHD -delta_t 200. -KK 30 -LL 40 -MM 10 -mres 4 -Pr 1. -Ek 1.0e-2 -Ra 65 -Ek_final 1.0e-3 -directory 'my_directory' -restart yes -restart_filename Restart_20000.b -dim_filename Dim_20000.b -dealiasing yes -max_newt 15 -max_gmres 1000 -restart_gmres 1000 -newt_eps 1.0e-7 -newt_delta 1.0e-16 -tol_gmres 1.0e-10 -M_wave 4 -adapt_param y -Nopt 5 -gamma 100. -delta_param 1. -time_step fbe -solver continuation_convective_implicit