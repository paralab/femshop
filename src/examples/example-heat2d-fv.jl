if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVheat2d");

useLog("FVheat2dlog")

# Configuration setup
domain(2)
solverType(FV)
timeStepper(EULER_EXPLICIT)

# Mesh
n = 20 # number of elements
mesh(QUADMESH, elsperdim=n)

# Variables and BCs
u = variable("u", SCALAR, CELL)
boundary(u, 1, FLUX, "0")

# Time interval and initial condition
T = 0.1;
timeInterval(T)
initial(u, "(sin(pi*x)*sin(2*pi*y))^4")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("D", 0.1) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
fluxAndSource(u, "-D*dot(grad(u),normal())", "0") 

#exportCode("fvheat2dcode") # uncomment to export generated code to a file
#importCode("fvheat2dcode") # uncomment to import code from a file

solve(u)

@finalize()

##### Uncomment below to plot
# xy = Femshop.fv_info.cellCenters

# using Plots
# pyplot();
# display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
