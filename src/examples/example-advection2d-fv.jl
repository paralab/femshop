if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVadvection2d");

@useLog("FVadvection2dlog")

# Configuration setup
@domain(2)
@solver(FV)
@stepper(EULER_EXPLICIT)

# Mesh
n = 15 # number of elements in each direction
@mesh(QUADMESH, n, 4)

# Variables and BCs
@variable(u)
@boundary(u, 1, FLUX, "(abs(y-0.2) < 0.11) ? sin(2*pi*t)^2 : 0") # x=0
@boundary(u, 2, NO_BC) # x=1
@boundary(u, 3, NO_BC) # y=0
@boundary(u, 4, NO_BC) # y=1

# Time interval and initial condition
T = 1;
@timeInterval(T)
@initial(u, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
@coefficient(a, VECTOR, ["cos(pi*x/2)","sin(pi*x/2)"]) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
@fluxAndSource(u, "upwind(a,u)", "0") 

#@exportCode("fvad2dcode") # uncomment to export generated code to a file
#@importCode("fvad2dcode") # uncomment to import code from a file

solve(u)

@finalize()

##### Uncomment below to plot

# xy = Femshop.fv_info.cellCenters

# using Plots
# pyplot();
# display(plot(xy[1,:], xy[2,:], u.values[:], st=:surface))
