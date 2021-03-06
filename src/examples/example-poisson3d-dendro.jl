#=
# 3D Poisson, Dirichlet bc
# uses dendro
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("poisson3ddendro");

# Try making an optional log
@useLog("poisson3ddendrolog")

# Generate C++ for dendro
@generateFor(DENDRO, "poisson3d", "This is an example for Poisson, Dirichlet bc.")
dendro(max_depth=6, wavelet_tol = 0.1, partition_tol = 0.3, solve_tol = 1e-6, max_iters = 100);

# Set up the configuration (order doesn't matter)
@domain(3)
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, "0")

# Write the weak form 
@coefficient(f, "-14*pi*pi*sin(3*pi*x)*sin(2*pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u),grad(v)) - f*v")

solve(u);

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
