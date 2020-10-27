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

cachesim(true);

# Set up the configuration (order doesn't matter)
@domain(1)                          # dimension
@solver(CG)                         # Use CG solver
@functionSpace(LEGENDRE, 4)         # basis function, order
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(LINEMESH, 5,2)               # build uniform LINEMESH. 2nd arg=# of elements, (optional)3rd arg=# of BIDs

@variable(u)                        # same as @variable(u, SCALAR)
@testSymbol(v)                      # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)       # boundary condition for BID 1
@boundary(u, 2, DIRICHLET, -1)      # and BID 2

# Write the weak form 
@coefficient(f, "-2.25*pi*pi*sin(1.5*pi*x)")
@weakForm(u, "-grad(u)*grad(v) - f*v")

solve(u);

# exact solution is sin(1.5*pi*x)
# check error
maxerr = 0;
exact(x) = sin(1.5*pi*x);

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
