# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
# Try making an optional log
@useLog("logfile")

# Set up the configuration (order doesn't matter)
@domain(1, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(DG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh("line.msh")                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)
@variable(q)

T = 0.2;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial("sin(pi*x)")               # initial condition needed if time dependent

#============== not yet implemented ======================

#@boundary(DIRICHLET, "sin(pi*t)")           # Specify like this to apply b.c. everywhere
#@boundary(DIRICHLET, "sin(pi*t)", bdry)    # Specify like this to apply b.c. where bdry(x) = true
#@boundary(1, DIRICHLET, "sin(pi*t)")       # Specify like this to apply b.c. with bid

# @trialFunction(u)                   # This may not be necessary. Let's think about it.
# @trialFunction(q)

@testFunction(v)                    # define a new symbol that will be used when writing the weak form

# Make an expression with trial and test functions using
# Dt, Dx, grad, etc. and surface() representing the surface integral
# This example is a heat equation Dt(u) = Dx(a*Dx(u)) , u(bdry) = 0, a = 1
# Split into Dt(u) = Dx(q) , q = Dx(u) , u(bdry) = 0, q(bdry) = 0
@weakForm(u, Dt(u*v) + q*Dx(v) + surface(q,v))
@weakForm(q, q*v + u*Dx(v) + surface(u,v))
# -OR-
@LHS(u, DT(u*v))
@RHS(u, -q*Dx(v) - surface(q,v))

@LHS(q, q*v)
@RHS(q, -u*Dx(v) - surface(u,v))

=======================================================#

solve();                            # solves everything

# solution is stored in the variable's "values"
#using Plots
#plot(Femshop.DGSolver.allnodes, u.values)

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
