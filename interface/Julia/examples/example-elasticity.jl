#=
# Linear elasticity
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("elasticity");

# Try making an optional log
@useLog("elasticitylog")

n = 10;
ord = 2;

# Set up the configuration (order doesn't matter)
@domain(3, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # Use CG solver
@functionSpace(LEGENDRE, ord)       # basis function, order
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(HEXMESH, n, 2)               # build uniform mesh. 2nd arg=# of elements, (optional)3rd arg=# of BIDs

@variable(u, VECTOR)                        # same as @variable(u, SCALAR)

@testSymbol(v, VECTOR)                      # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, [0,0,0])       # x=0
@boundary(u, 2, NEUMANN, [0,0,0])

# Write the weak form
#@coefficient(mu, "x>0.5 ? 1 : 10") # discontinuous mu
@coefficient(mu, 10) # constant mu
@coefficient(lambda, 1.25)
@coefficient(f, VECTOR, ["0","0","-0.1"])
@weakForm(u, "inner( (lambda * div(u) .* [1 0 0; 0 1 0; 0 0 1] + mu .* (grad(u) + transpose(grad(u)))), grad(v)) - dot(f,v)")

println("solving")
solve(u);
println("solved")

# using Plots
# pyplot();
# N = n*ord+1;
# half = Int(round(N/2));
# range = (N*N*half+1):(N*N*(half+1));
# hair = (N*N*half+1):(N*N*half+N);
# #display(plot(Femshop.grid_data.allnodes[hair,1], u.values[hair,3], reuse=false))
# #display(plot(Femshop.grid_data.allnodes[range,1], Femshop.grid_data.allnodes[range,2], u.values[range,3], st=:surface, reuse=false))
# bent = u.values[hair,3];
# tmp = copy(bent);
# xpos = Femshop.grid_data.allnodes[hair,1];
# slope = 0;
# for i=2:length(bent)
#     change = tmp[i]-tmp[i-1];
#     bent[i] = bent[i-1] + slope + change;
#     global slope = (bent[i] - bent[i-1]);
# end
# display(plot(Femshop.grid_data.allnodes[hair,1], bent, reuse=false, marker=:circle))

# check
log_dump_config();
log_dump_prob();

@finalize()
