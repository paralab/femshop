# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
init_femshop("exampleProject");

# Try making an optional log
@useLog("exampleLog")

# Set up the configuration (order doesn't matter)
@domain(1, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(DG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh("line.msh")                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)
@variable(q)

@testFunction(v)                    # sets the symbol for a test function

T = 0.01;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "abs(x-0.5)<0.2 ? 1 : 0")               # initial condition needed if time dependent

@boundary(u, 1, DIRICHLET, "0")       # Specify like this to apply b.c. with bid
@boundary(q, 1, NEUMANN, "0") 

# Make an expression with trial and test functions using
# Dt, grad, etc. and surface() representing the surface integral
# This example is a 1D heat equation Dt(u) = Dx(a*Dx(u)) , u(bdry) = 0, a = 1
# Split into Dt(u) = grad(q) , q = grad(u) , u(bdry) = 0, q(bdry) = 0
@weakForm(u, "Dt(u*v) + grad(q)*v - surface(diff_flux(q)*v)")  # These are changed into LHS and RHS expressions with 
@weakForm(q, "q*v + grad(u)*v - surface(diff_flux(u)*v)")      # symbolic operators that need to be implemented by the solvers

###### Not ready #####
# solve();                            # solves everything

# solution is stored in the variable's "values"
#using Plots
#display(plot(Femshop.DGSolver.allnodes, u.values))
######################

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
