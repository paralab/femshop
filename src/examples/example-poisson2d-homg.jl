#=
# 2D Poisson, Dirichlet bc
# Uses HOMG(MATLAB)
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("poisson2dmatlab");

# Try making an optional log
@useLog("poisson2dmatlablog")

# Generate Matlab
@generateFor(HOMG, "poisson2d", "This is an example for 2D Poisson, Dirichlet bc.")

# Set up the configuration (order doesn't matter)
@domain(2, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, 30)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, "sin(3*pi*x)")

# Write the weak form 
@coefficient(f, "-2*pi*pi*sin(pi*x)*sin(pi*y)")
@weakForm(u, "-dot(grad(u),grad(v)) - f*v")

solve(u);

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
