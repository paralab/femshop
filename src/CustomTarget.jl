#=
This is a wrapper for the custom generation functions supplied elsewhere
=#
module CustomTarget

import ..Femshop: JULIA, CPP, MATLAB, DENDRO, HOMG, CUSTOM_GEN_TARGET,
            SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, 
            CG, DG, HDG, FV,
            NODAL, MODAL, CELL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, 
            NONLINEAR_NEWTON, NONLINEAR_SOMETHING, 
            EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, ABM4, 
            DEFAULT_SOLVER, PETSC, 
            VTK, RAW_OUTPUT, CUSTOM_OUTPUT, 
            DIRICHLET, NEUMANN, ROBIN, NO_BC, FLUX,
            MSH_V2, MSH_V4,
            SCALAR, VECTOR, TENSOR, SYM_TENSOR,
            LHS, RHS,
            LINEMESH, QUADMESH, HEXMESH
import ..Femshop: Femshop_config, Femshop_prob, GenFunction, Variable, Coefficient
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, refel, mesh_data, grid_data, genfunctions, variables, coefficients, 
        test_functions, linears, bilinears, time_stepper
import ..Femshop: SymExpression, SymEntity
import ..Femshop: CachesimOut, use_cachesim
import ..Femshop: custom_gen_funcs

# There are three necessary functions. These are temporary placeholders
function external_get_language_elements_function() return (".jl", "#", ["#=", "=#"]) end;
function external_generate_code_layer_function(var, entities, terms, lorr, vors) return ("","") end;
function external_generate_code_files(lhs_vol, lhs_surf, rhs_vol, rhs_surf) return 0 end;

function include_custom_target_file(file)
    include(file);
    
    
end

end # module