#=
A struct containing problem information
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_prob
    # Boundary condition info
    bc_type::Array{String,2}        # DIRICHLET, etc. for each variable and bid
    bid::Array{Int,2}               # BID for each variable and boundary section
    bc_func::Array{Any,2}           # GenFunctions or numbers for each variable and bid
    
    # Time dependence info
    time_dependent::Bool            # Is this problem time dependent
    end_time::Float64               # If so, what is the final time
    initial::Array{Any,2}           # An array of initial condition GenFunctions, one for each variable
    lhs_time_deriv::Array{Bool,1}   # (may be removed) signals LHS time derivatives for each variable
    
    # Constructor builds an empty prob.
    Femshop_prob() = new(
        Array{String,2}(undef,(0,0)), 
        Array{Int,2}(undef,(0,0)), 
        Array{Any,2}(undef,(0,0)),
        false, 
        0, 
        Array{Any,2}(undef,(0,0)), 
        Array{Bool,1}(undef,(0))
    );
end
