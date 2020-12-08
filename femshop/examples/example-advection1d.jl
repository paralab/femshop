#=
# 1D advection matching the example in the book.
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("advection1d");

@useLog("advection1dlog")

@domain(1)
@solver(DG)
@functionSpace(LEGENDRE, 3)
@stepper(LSRK4)

@mesh(LINEMESH, 10, 2, [0,2])

@variable(u)
@testSymbol(v)

T = 1.0;
@timeInterval(T)
@initial(u, "sin(x)")

@boundary(u, 1, DIRICHLET, "-sin(2*pi*t)")
@boundary(u, 2, NO_BC)

@coefficient(a, 2*pi)
@coefficient(alpha, 1)
#@weakForm(u, "Dt(u*v) + a*grad(u)*v - surface(a*normal()*u*v - a*normal()*ave(u)*v - 0.5*(1-alpha)*a*jump(u)*v)")
@weakForm(u, "Dt(u*v) + a*grad(u)*v - surface(a*normal()*u*v - a*normal()*ave(u)*v)")

solve(u);

# solution is stored in the variable's "values"
using Plots
pyplot();
display(plot(Femshop.grid_data.allnodes[1,:], u.values[1,:], marker=:circle, reuse=false))

@finalize()
