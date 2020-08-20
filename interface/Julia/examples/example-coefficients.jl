#=
# Tests out using a variety of coefficients in different places
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("coef");

# Try making an optional log
@useLog("coeflog")

# Set up the configuration
@domain(2, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, 20)                 # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, 0)

# Write the weak form 
@coefficient(a, "x")
@coefficient(b, "3")
@coefficient(c, "y")
@coefficient(d, 0.1)
@coefficient(f, "(-5*pi*pi*x*sin(pi*x)*sin(2*pi*y) + pi*cos(pi*x)*sin(2*pi*y)) + y*sin(pi*x)*sin(2*pi*y)")

@weakForm(u, "-a*b*dot(grad(u), grad(v))*d + c*u*d*v*b - d*f*b*v")

solve(u);

# exact solution is sin(pi*x)*sin(2*pi*y)
# check error
maxerr = 0;
exact(x,y) = sin(pi*x)*sin(2*pi*y);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    y = Femshop.grid_data.allnodes[i,2];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes[:,1], Femshop.grid_data.allnodes[:,2], u.values, st=:surface))

# check
log_dump_config();
log_dump_prob();

@finalize()
