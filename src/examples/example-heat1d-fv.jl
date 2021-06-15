if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVheat1d");

useLog("FVheat1dlog")

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4)

# Mesh
n = 40 # number of elements
mesh(LINEMESH, elsperdim=n)

# Variables and BCs
u = variable("u", SCALAR, CELL)
boundary(u, 1, FLUX, "0")

# Time interval and initial condition
T = 0.1;
timeInterval(T)
initial(u, "x<0.5 ? sin(pi*x)^4 : 0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
coefficient("D", 0.1) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
fluxAndSource(u, "-D*dot(grad(u),normal())", "0") 

#exportCode("fvheat1dcode") # uncomment to export generated code to a file
#importCode("fvheat1dcode") # uncomment to import code from a file

solve(u)

@finalize()

##### Uncomment below to plot
# # The initial condition
# u0 = zeros(n);
# x = Femshop.fv_info.cellCenters[:]
# for i=1:n
#     u0[i] = x[i]<0.5 ? sin(pi*x[i])^4 : 0
# end

# using Plots
# pyplot();
# display(plot([x x], [u0 u.values[:]], markershape=:circle, label=["initial" "t="*string(T)]))
