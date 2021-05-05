#=
# Import a simple tet or hex mesh from a .msh file.
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("unstruct3d");

@useLog("unstruct3dlog")

@domain(3, SQUARE, UNSTRUCTURED)
@functionSpace(LEGENDRE, 1)

# Unit cube
#@mesh("utet.msh")  # Using tets
@mesh("uhex.msh")   # Using hexes

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, "sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

@coefficient(K, "x*x + 4")
@coefficient(C, 10)
@coefficient(f, "x*pi*cos(pi*0.5*x)*sin(pi*y)*sin(pi*z) + (-10-(x*x+4)*pi*pi*2.25)*sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)")

# PDE: div(K*grad(u)) -C*u = f
@weakForm(u, "K*dot(grad(u), grad(v)) + C*u*v + f*v")

solve(u);

# exact solution is sin(pi*0.5*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
exact(x,y,z) = sin(pi*0.5*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Femshop.grid_data.allnodes,2)
    x = Femshop.grid_data.allnodes[1,i];
    y = Femshop.grid_data.allnodes[2,i];
    z = Femshop.grid_data.allnodes[3,i];
    err = abs(u.values[i] - exact(x,y,z));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));

# using Plots
# pyplot();
# display(plot(Femshop.grid_data.allnodes[1,:], Femshop.grid_data.allnodes[2,:], u.values[:], st=:surface))

@finalize()
