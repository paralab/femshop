# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
# Try making an optional log
@useLog("logfile")

# Set up the configuration (order doesn't matter)
@domain(2, IRREGULAR, UNSTRUCTURED) # dimension, geometry, decomposition
@mesh("circle_2d.msh")              # .msh file or generate our own
@solver(DG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

#============== not yet implemented ======================
# Specify the problem
T = 3;
a = 1;
@timeInterval(0, T)                 # (start, end) using this sets problem to time dependent
@initial("sin(pi*x)")               # initial condition needed if time dependent

@boundary(DIRICHLET, "sin(pi*t)")   # Specify like this to apply b.c. everywhere
@boundary(DIRICHLET, 1.3, bdry)     # Specify like this to apply b.c. where bdry(x) = true

u = trial_function();               # Define the trial and test functions
v = test_function();

# Make an expression with trial and test functions using
# Dt, Dx, grad, etc. and surface() representing the surface integral
@weakForm(Dt(u*v) - a*u*grad(v) + surface(a*v))

@flux(LAXWENDROFF, 0.2)             # Specifies the numerical flux for hyperbolic problems
                                    # Second value is the relevant parameter.

solve(u);                           # solves everything
=======================================================#

# check
c = Femshop.config
log_dump_config(c);
@assert(c.geometry == IRREGULAR)
@assert(Femshop.mesh_data.nx > 0)

# Try making some files
@language(CPP, "tryitcpp", "Just an empty file for testing");
@finalize