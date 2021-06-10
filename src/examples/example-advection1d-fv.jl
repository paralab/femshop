if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVadvection1d");

@useLog("FVadvection1dlog")

# Configuration setup
@domain(1)
@solver(FV)
@stepper(RK4)

# Mesh
n = 40 # number of elements
@mesh(LINEMESH, n, 2)

# Variables and BCs
@variable(u)
@variable(v)
@boundary(u, 1, FLUX, "t<0.2 ? 1 : 0")
@boundary(u, 2, NO_BC)
@boundary(v, 1, FLUX, "t<0.2 ? 1 : 0")
@boundary(v, 2, NO_BC)

# Time interval and initial condition
T = 0.5;
@timeInterval(T)
@initial(u, "0")
@initial(v, "0")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
@coefficient(a, 1) # advection velocity
# The "upwind" function applies upwinding to the term (a.n)*u with flow velocity a.
# The optional third parameter is for tuning. Default upwind = 0, central = 1. Choose something between these.
@fluxAndSource([u, v], ["upwind(a,u)", "upwind(a,v,0.75)"], ["0", "0"]) 

#@exportCode("fvad1dcode") # uncomment to export generated code to a file
#@importCode("fvad1dcode") # uncomment to import code from a file

solve([u,v])

@finalize()

##### Uncomment below to compare to exact solution

# # The exact solution with constant velocity v
# a = 1;
# exact = zeros(n);
# x = Femshop.fv_info.cellCenters[:]
# for i=1:n
#     xt = x[i] - a*T;
#     if xt < 0
#         exact[i] = xt > -0.2 ? 1 : 0
#     end
# end

# using Plots
# pyplot();
# display(plot([x x x], [exact u.values[:] v.values[:]], markershape=:circle, label=["exact" "upwind" "alpha=0.75"]))
