#=
# 1D Poisson, Dirichlet bc
# CG, Linear element
# Simplest test possible
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("poisson1d");

# Try making an optional log
@useLog("poisson1dlog")

# Set up the configuration (order doesn't matter)
@domain(1, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(LINEMESH, 5)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testFunction(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-pi*pi*sin(pi*x)")
@weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# exact solution is sin(pi*x)
# check error
maxerr = 0;
exact(x) = sin(pi*x);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    err = abs(u.values[i] - exact(x));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes, u.values, markershape=:circle))

# check
log_dump_config();
log_dump_prob();

@finalize()
