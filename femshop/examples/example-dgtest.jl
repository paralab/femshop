#=
# DG test
# 1D Poisson using SIPG formulation for DG
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("test");

@useLog("dgtestlog")

@domain(1)
@solver(DG)
@functionSpace(LEGENDRE, 1)

@mesh(LINEMESH, 9)

@variable(u)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)

@coefficient(f, "-pi*pi*sin(pi*x)")
@coefficient(beta, 1)
@weakForm(u, "dot(grad(u),grad(v)) + f*v - surface( ave_normdotgrad(u) * jump(v) - jump(u) * ave_normdotgrad(v) + beta*jump(u)*jump(v) )")

solve(u);

# solution is stored in the variable's "values"
# using Plots
# pyplot();
# display(plot(Femshop.grid_data.allnodes[1,:], u.values[1,:], reuse=false))

@finalize()
