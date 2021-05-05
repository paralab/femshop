#=
# Import a simple triangle or quad mesh from a .msh file.
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("unstruct2dtest");

@useLog("unstruct2dlog")

@domain(2, SQUARE, UNSTRUCTURED)
@functionSpace(LEGENDRE, 2)

# This rectangle covers [0, 0.1]x[0, 0.3]
# Uncomment the desired mesh.
@mesh("utriangle.msh")  # Using triangles
#@mesh("uquad.msh")     # Using quads

@variable(u)
@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-(x+1)*200*pi*pi*sin(10*pi*x)*sin(10*pi*y) + 10*pi*cos(10*pi*x)*sin(10*pi*y)")
@coefficient(k, "x+1")
@weakForm(u, "-k*dot(grad(u), grad(v)) - f*v")

solve(u);

@finalize()

# exact solution is sin(10*pi*x)*sin(10*pi*y)
# check error
maxerr = 0;
exact(x,y) = sin(10*pi*x)*sin(10*pi*y);

for i=1:size(Femshop.grid_data.allnodes,2)
    x = Femshop.grid_data.allnodes[1,i];
    y = Femshop.grid_data.allnodes[2,i];
    err = abs(u.values[i] - exact(x,y));
    global maxerr;
    maxerr = max(err,maxerr);
end
println("max error = "*string(maxerr));
