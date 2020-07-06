#=
# 2D Poisson, Dirichlet bc, CG
=#
if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
init_femshop("poisson2d");

# Try making an optional log
@useLog("poisson2dlog")

# Set up the configuration (order doesn't matter)
@domain(2, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, 20)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testFunction(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# solution is stored in the variable's "values"
using Plots
display(plot(Femshop.grid_data.allnodes[:,1], Femshop.grid_data.allnodes[:,2], u.values, st = :surface))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
