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
@domain(1)                      # dimension
@functionSpace(LEGENDRE, 4)     # basis function, order

# Specify the problem (mesh comes first)
@mesh(LINEMESH, 20)             # build uniform LINEMESH. 2nd arg=# of elements, (optional)3rd arg=# of BIDs

@variable(u)                    # same as @variable(u, SCALAR)
@testSymbol(v)                  # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, "0")   # boundary condition for BID 1 is Dirichlet with value 0

# Write the weak form 
@coefficient(f, "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
@weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# exact solution is sin(10*pi*x)*sin(pi*x)
# check error
maxerr = 0;
exact(x) = sin(10*pi*x)*sin(pi*x);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    err = abs(u.values[i] - exact(x));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Femshop.grid_data.allnodes, u.values, markershape=:circle, legend=false))

# Dump things to the log if desired
log_dump_config();
log_dump_prob();

@finalize()
