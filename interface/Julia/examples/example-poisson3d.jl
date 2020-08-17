#=
# 1D Poisson, Dirichlet bc
# CG, Linear element
# Simplest test possible
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("poisson3d");

# Try making an optional log
@useLog("poisson3dlog")

n = 10;
ord = 3;

# Set up the configuration (order doesn't matter)
@domain(3, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, ord)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(HEXMESH, n)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testFunction(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

solve(u);

# exact solution is sin(pi*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
exact(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    y = Femshop.grid_data.allnodes[i,2];
    z = Femshop.grid_data.allnodes[i,3];
    err = abs(u.values[i] - exact(x,y,z));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#N = n*ord+1;
#half = Int(round(N/2));
#range = (N*N*half+1):(N*N*(half+1));
#display(plot(Femshop.grid_data.allnodes[range,1], Femshop.grid_data.allnodes[range,2], u.values[range], st=:surface))

# check
log_dump_config();
log_dump_prob();

@finalize()
