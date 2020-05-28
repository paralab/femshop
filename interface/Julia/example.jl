# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
# Try making an optional log
@useLog("logfile")

# Set up the configuration (order doesn't matter)
#@domain(2, IRREGULAR, UNSTRUCTURED) # dimension, geometry, decomposition
@domain(3, IRREGULAR, UNSTRUCTURED) # dimension, geometry, decomposition
@solver(DG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Problem specification
@mesh("circle_2d.msh")              # .msh file or generate our own

T = 3;
@timeInterval(0, T)                 # (start, end) using this sets problem to time dependent
@initial("sin(pi*x)")               # initial condition needed if time dependent

@variable(u)
@variable(q)

#songzhe: maybe start from a simple square case, we also need to tell boundary ids.
@boundary(1, DIRICHLET, "sin(pi*t)")   # Specify like this to apply b.c. with bid
#============== not yet implemented ======================

# Make an expression with trial and test functions using
# Dt, Dx, grad, etc. and surface() representing the surface integral
@weakForm(Dt(u*v) - a*u*grad(v) + surface(a*v))

@flux(LAXWENDROFF, 0.2)             # Specifies the numerical flux for hyperbolic problems
                                    # Second value is the relevant parameter.

solve(u);                           # solves everything
=======================================================#

# check
log_dump_config(Femshop.config);

# Try making some files
#@language(CPP, "tryitcpp", "Just an empty file for testing");
@finalize
