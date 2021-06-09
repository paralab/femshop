if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVadvection1d");

@useLog("FVadvection1dlog")

# Configuration setup
@domain(1)
@solver(FV)
@stepper(EULER_EXPLICIT)

# Mesh
n = 30
@mesh(LINEMESH, n)

# Variable
@variable(u)

# Time interval and initial condition
T = 0.5;
@timeInterval(T)
@initial(u, "sin(pi*x)^5")

# The flux and source terms of the conservation equation
# F and S in the following equation:
# Dt(int(u dx)) = int(S dx) - int(F.n ds)
@coefficient(v, 1) # advection velocity
# the "upwind" function applies upwinding to the term v*u with flow velocity v
@fluxAndSource(u, "upwind(v,u)", "0") 

#@exportCode("fvad1dcode")
#@importCode("fvad1dcode")

solve(u)

@finalize()

##### Uncomment below to compare to exact solution

# # The exact solution with constant velocity v
# v = 1;
# exact = zeros(n);
# x = Femshop.fv_info.cellCenters[:]
# for i=1:n
#     exact[i] = sin(pi*(x[i] - v*T))^6
# end

# using Plots
# pyplot();
# display(plot([x x], [u.values[:] exact], markershape=:circle, legend=false))
