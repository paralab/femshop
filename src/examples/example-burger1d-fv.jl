if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("FVburger1d");

useLog("FVburger1dlog", level=3)

# Configuration setup
domain(1)
solverType(FV)
timeStepper(RK4, cfl=.1)

# Mesh
n = 50 # number of elements
mesh(LINEMESH, elsperdim=n, bids=2)

# Variables and BCs
u = variable("u", SCALAR, CELL)
boundary(u, 1, FLUX, 0)
boundary(u, 2, NO_BC)

v = variable("v", SCALAR, CELL)
boundary(v, 1, FLUX, 0)
boundary(v, 2, NO_BC)

# Time interval and initial condition
T = 0.2;
timeInterval(T)
initial(u, "sin(pi*2*x)")
initial(v, "sin(pi*2*x)")
# initial(u, "x>0.1 && x<0.5 ? 1 : 0")
# initial(v, "x>0.1 && x<0.5 ? 1 : 0")

# The flux and source terms of the conservation equation
flux([u,v], ["0.5*burgerGodunov(u,u*u)", "0.5*central(v*v)"]) 

# printLatex(u)

#@exportCode("fvad1dcode") # uncomment to export generated code to a file
#@importCode("fvad1dcode") # uncomment to import code from a file

solve([u,v])

finalize_femshop()

##### Uncomment below to compare to plot

x = Femshop.fv_info.cellCenters[:]
n = length(x);
ic = zeros(n);
for i=1:n
    # ic[i] = (x[i]>0.1 && x[i]<0.5) ? 1 : 0
    ic[i] = sin(2*pi*x[i])
end

using Plots
pyplot();
display(plot([x x x], [u.values[:] v.values[:] ic], markershape=:circle, label=["Godunov" "central" "initial"]))
