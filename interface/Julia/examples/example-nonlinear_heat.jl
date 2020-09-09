#=
# 1D nonlinear heat eq. Dirichlet bc, CG
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("nonlinear_heat");

# Try making an optional log
@useLog("nonlinear_heatlog")

# Set up the configuration (order doesn't matter)
@domain(1, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(EULER_IMPLICIT)            # time stepper (optional second arg is CFL#)

# Specify the problem
@mesh(LINEMESH, 15, 2)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)
@variable(du)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

T = 1;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "(1.31-0.866)*x-1.31+0.866*2")  # initial condition needed if time dependent
@initial(du, "0")  # initial condition needed if time dependent

@boundary(du, 1, DIRICHLET, 0)
@boundary(du, 2, DIRICHLET, 0)
@boundary(u, 1, DIRICHLET, 0.866)
@boundary(u, 2, DIRICHLET, 1.31)

# Write the weak form
@coefficient(f, "x^(-4)")
@weakForm(du, "(grad(v)*grad(du) - v*5*u*u*u*u*f*du) + (grad(v)*grad(u) - v*u*u*u*u*u*f)")

solve(du);

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes[:,1], Femshop.grid_data.allnodes[:,2], u.values, st = :surface))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
