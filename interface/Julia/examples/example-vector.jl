#=
# Tests vector unknown capability with a simple Poisson-like problem.
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("vector");

# Try making an optional log
@useLog("vectorlog")

# Set up the configuration
@domain(2, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 3)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# Specify the problem
@mesh(QUADMESH, 30)                 # .msh file or generate our own

@variable(u, VECTOR)

@testSymbol(v, VECTOR)

@boundary(u, 1, DIRICHLET, [0, 0])

# Write the weak form
@coefficient(f, VECTOR, ["-5*pi*pi*sin(pi*x)*sin(2*pi*y)", "-25*pi*pi*sin(3*pi*x)*sin(4*pi*y)"])

@weakForm(u, "-inner(grad(u), grad(v)) - f*v")

solve(u);

# # exact solution is [sin(pi*x)*sin(2*pi*y), sin(3*pi*x)*sin(4*pi*y)]
# # check error
erroru = zeros(size(u.values));
maxerru = 0
exactu1(x,y) = sin(pi*x)*sin(2*pi*y);
exactu2(x,y) = sin(3*pi*x)*sin(4*pi*y);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    y = Femshop.grid_data.allnodes[i,2];
    exac = [exactu1(x,y), exactu2(x,y)];
    for j=1:size(Femshop.grid_data.allnodes,2)
        erroru[i,j] = u.values[i,j] - exac[j];
        global maxerru;
        maxerru = max(abs(erroru[i,j]),maxerru);
    end
    
end
println("u max error = "*string(maxerru));

# using Plots
# pyplot();
# display(plot(Femshop.grid_data.allnodes[:,1], Femshop.grid_data.allnodes[:,2], u.values[:,1], st=:surface))

# check
log_dump_config();
log_dump_prob();

@finalize()
