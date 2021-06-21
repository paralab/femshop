#=
# Import a simple triangle or quad mesh from a .msh file.
=#

### If the Femshop package has already been added, use this line #########
using Femshop # Note: to add the package, first do: ]add "https://github.com/paralab/femshop.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Femshop)
#     include("../Femshop.jl");
#     using .Femshop
# end
##########################################################################

init_femshop("unstruct2dtest");

useLog("unstruct2dlog")

domain(2, grid=UNSTRUCTURED)
functionSpace(order=2)

# This rectangle covers [0, 0.1]x[0, 0.3]
# Uncomment the desired mesh.
mesh("utriangle.msh")  # Using triangles
#mesh("uquad.msh")     # Using quads

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)

coefficient("f", "(-10-(x+1)*200*pi*pi)*sin(10*pi*x)*sin(10*pi*y) + 10*pi*cos(10*pi*x)*sin(10*pi*y)")
coefficient("k", "x+1")
coefficient("C", 10)
weakForm(u, "k*dot(grad(u), grad(v)) + C*u*v+ f*v")

solve(u);

finalize_femshop();

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
