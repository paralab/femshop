#=
# This file should be generated.
# It contains the RHS function parts to be passed to the solver's rhs function
=#

struct RHSParams
    solve_order::Array{Int,1}       # The order in which to solve equations(array of variable indices)
    dependence::Array{Int,2}        # What variables does each eq. depend on
    rhs_eq::Array                   # Functions corresponding to each eq
    bdry_ind::Array{Int,1}          # Boundary indices
    bdry_type::Array                # Boundary types for each variable and boundary index
    bdry_func::Array                # Functions for each boundary condition
end

############# This is the thing that needs to be generated ##################
function build_rhs_params()
    solve_order = [2;1];
    dependence = [
        0 1 ;
        1 0
    ]
    function eq_1(rhsv, dif, time)
        return dg_mass_inv_advective(rhsv[:,:,2]) .- dg_surface_int(dif[:,:,2]);
    end
    function eq_2(rhsv, dif, time)
        return dg_mass_inv_advective(rhsv[:,:,1]) .- dg_surface_int(dif[:,:,1]);
    end
    rhs_eq = [
        eq_1 ;
        eq_2
    ]
    
    function bdry_u_1(t)
        return 0;
    end
    function bdry_q_1(t)
        return 0;
    end
    bdry_ind = [1];
    bdry_type = [
        DIRICHLET ;
        NEUMANN
    ]
    bdry_func = [
        bdry_u_1 ;
        bdry_q_1
    ]
    
    return RHSParams(solve_order, dependence, rhs_eq, bdry_ind, bdry_type, bdry_func);
end
