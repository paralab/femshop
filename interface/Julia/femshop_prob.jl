#=
A struct containing problem information
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_prob
    mesh_dofs::Int          # DOFs corresponding to the mesh, not variables
    
    # Domain
    bc_type::String
    bid::Int
    
    # Time dependent info
    time_dependent::Bool
    end_time::Float64
    initial                 # an array of GenFunctions, one for each variable
    
    

    # Constructor builds a default prob.
    Femshop_prob() = new(
        0,
        DIRICHLET,
        1,
        false,
        0,
        []
    );
end
