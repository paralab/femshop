#=
# This file should be generated.
# It contains the RHS function parts to be passed to the solver's rhs function
=#
export RHSParams

mutable struct RHSParams
    solve_order::Array{Int,1}       # The order in which to solve equations(array of variable indices)
    dependence::Array{Int,2}        # What variables does each eq. depend on
    rhs_eq::Array                   # Functions corresponding to each eq
    bdry_ind::Array                 # Boundary indices
    bdry_type::Array                # Boundary types for each variable and boundary index
    bdry_func::Array                # Functions for each boundary condition
end

# ############# This is the thing that needs to be generated ##################
# function build_rhs_params()
#     # These are set in the generated file.
#     # solve_order = [1];
#     # dependence = [0 0; 0 0];
#     # rhs_eq = [];
    
#     # Include the generated file
#     include(output_dir*"/femshop_generated_rhs_parameters.jl");
#     #(solve_order, dependence, rhs_eq) = Base.invokelatest(get_generated_parameters, );
    
#     #return RHSParams(solve_order, dependence, rhs_eq, bdry_ind, bdry_type, bdry_func);
#     return Base.invokelatest(get_generated_parameters, );
# end
