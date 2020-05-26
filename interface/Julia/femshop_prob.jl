#=
A struct containing problem information
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_prob
    # Domain
    bc_type::String
    bid::Int
    #songzhe: not sure how to add expressions for bc yet

    # Constructor builds a default config.
    Femshop_prob() = new(
        DIRICHLET,
        1
    );
end
