#=
# 2D heat eq. Dirichlet bc, CG
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("heat2d");

# Try making an optional log
@useLog("heat2dlog")

# Set up the configuration (order doesn't matter)
@domain(2, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(CRANK_NICHOLSON)            # time stepper (optional second arg is CFL#)

# Specify the problem
@mesh(QUADMESH, 10)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

T = 1;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "abs(x-0.5)+abs(y-0.5) < 0.2 ? 1 : 0")  # initial condition needed if time dependent

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form
@coefficient(f, "0.5*sin(6*pi*x)*sin(6*pi*y)")
@weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")

solve(u);

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Femshop.grid_data.allnodes[1,:], Femshop.grid_data.allnodes[2,:], u.values[:], st = :surface))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
