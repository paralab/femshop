#=
# test a simple triangle mesh from a file
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("tritest");

@useLog("tritestlog")

@domain(2)
@functionSpace(LEGENDRE, 1)

@mesh("src/examples/square.msh")

@finalize()
