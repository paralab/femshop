#=
A struct containing problem information
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_prob
    mesh_dofs::Int          # DOFs corresponding to the mesh, not variables
    
    # Domain
    bc_type::Array{String,2}
    bid::Array{Int,2}
    bc_func::Array{Any,2}
    
    # Time dependent info
    time_dependent::Bool
    end_time::Float64
    initial                 # an array of GenFunctions, one for each variable
    lhs_time_deriv::Array{Bool,1}
    

    # Constructor builds a default prob.
    Femshop_prob() = new(
        0,
        Array{String,2}(undef,(1,1)),
        Array{Int,2}(undef,(1,1)),
        Array{Any,2}(undef,(1,1)),
        false,
        0,
        [],
        []
    );
end
