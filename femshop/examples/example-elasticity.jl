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

n = [10,4,4];
interval = [0,1,0,0.2,0,0.2];
ord = 2;

# Set up the configuration (order doesn't matter)
@domain(3, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # Use CG solver
@functionSpace(LEGENDRE, ord)       # basis function, order
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(HEXMESH, n, 2, interval)               # build uniform mesh. 2nd arg=# of elements, (optional)3rd arg=# of BIDs

@variable(u, VECTOR)                        # same as @variable(u, SCALAR)

@testSymbol(v, VECTOR)                      # sets the symbol for a test function

@boundary(u, 1, DIRICHLET, [0,0,0])       # x=0
@boundary(u, 2, NEUMANN, [0,0,0])

# Write the weak form
@coefficient(mu, "x>0.5 ? 0.2 : 10") # discontinuous mu
#@coefficient(mu, 1) # constant mu
@coefficient(lambda, 1.25)
@coefficient(f, VECTOR, ["0","0","-0.1"])
@weakForm(u, "inner( (lambda * div(u) .* [1 0 0; 0 1 0; 0 0 1] + mu .* (grad(u) + transpose(grad(u)))), grad(v)) - dot(f,v)")

println("solving")
solve(u);
println("solved")

# using Plots
# pyplot();
# Nx = n[1]*ord+1;
# Ny = n[2]*ord+1;
# Nz = n[3]*ord+1;
# half = Int(round(Nz/2));
# range = (Nx*Ny*half+1):(Nx*Ny*(half+1));
# hair = (Nx*Ny*half+1):(Nx*Ny*half+Nx);
# #display(plot(Femshop.grid_data.allnodes[1,hair], u.values[3,hair], reuse=false))
# #display(plot(Femshop.grid_data.allnodes[1,range], Femshop.grid_data.allnodes[2,range], u.values[3,range], st=:surface, reuse=false))
# bent = u.values[3,hair];
# tmp = copy(bent);
# xpos = Femshop.grid_data.allnodes[1,hair];
# slope = 0;
# for i=2:length(bent)
#     change = tmp[i]-tmp[i-1];
#     bent[i] = bent[i-1] + slope + change;
#     global slope = (bent[i] - bent[i-1]);
# end
# display(plot(Femshop.grid_data.allnodes[1,hair], bent, reuse=false, marker=:circle))

# check
log_dump_config();
log_dump_prob();

@finalize()
