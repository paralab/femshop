#=
# 1D advection
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("advection1d");

@useLog("advection1dlog")

@domain(1)
@solver(DG)
@functionSpace(LEGENDRE, 1)
@stepper(EULER_EXPLICIT)

@mesh(LINEMESH, 20)

@variable(u)
@testSymbol(v)

T = 1;
@timeInterval(T)
@initial(u, "abs(x-0.5) < 0.3 ? 1 : 0")

@boundary(u, 1, NO_BC, 0)

@coefficient(b, 0.2)
@weakForm(u, "Dt(u*v) - b*u*grad(v) + surface(b*ave(u)*v + 0.5*b*jump(u)*v)")

solve(u);

# solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Femshop.grid_data.allnodes[1,:], u.values[1,:], reuse=false))

@finalize()
