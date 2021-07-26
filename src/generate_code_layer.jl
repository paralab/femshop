#=
Use the symbolic layer expressions to generate the code
Redirects to solver and target specific functions.
=#

# lorr = LHS or RHS, vors = volume or surface
function generate_code_layer(symex, var, lorr, vors, solver, language, framework)
    # Count variables, dofs, and store offsets
    multivar = typeof(var) <:Array;
    varcount = 1;
    dofsper = 0;
    offset_ind = [0];
    if multivar
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        dofsper = length(var[1].symvar);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + length(var[i].symvar);
        end
    else
        dofsper = length(var.symvar);
    end
    
    # symex is an array of arrays of SymExpressions which are Expr trees with SymEntities as leaves. (array for variable components, terms)
    # In the case of muliple variables, it's an array of arrays of arrays. (variables, components, terms)
    # The Symexpression contains a list of all leaves that need to be evaluated before combining.
    # First collect all of them and eliminate identical ones.
    entities = []
    if multivar
        for vi=1:length(symex)
            for ci=1:length(symex[vi])
                for ti=1:length(symex[vi][ci])
                    for i=1:length(symex[vi][ci][ti].entities) # loop over entities for this variable/component
                        entity_present = false;
                        for j=1:length(entities) # check against existing entities
                            if is_same_entity(symex[vi][ci][ti].entities[i], entities[j])
                                entity_present = true;
                                break;
                            end
                        end
                        if !entity_present
                            push!(entities, symex[vi][ci][ti].entities[i]); # add it if unique
                        end
                    end
                end
            end
        end
    else # same thing as above, but for symex rather than symex[vi]
        for ci=1:length(symex)
            for ti=1:length(symex[ci])
                for i=1:length(symex[ci][ti].entities) # loop over entities for this variable/component
                    entity_present = false;
                    for j=1:length(entities) # check against existing entities
                        if is_same_entity(symex[ci][ti].entities[i], entities[j])
                            entity_present = true;
                            break;
                        end
                    end
                    if !entity_present
                        push!(entities, symex[ci][ti].entities[i]); # add it if unique
                    end
                end
            end
        end
    end
    
    # Determine if derivative matrices will be required
    need_derivs = false;
    need_deriv_matrix_for_fv = false;
    for i=1:length(entities)
        if length(entities[i].derivs) > 0
            # FV only needs the matrices for coefficient derivatives
            if config.solver_type == FV
                (ctype, cval) = get_coef_val(entities[i]);
                if ctype == 2 # a coefficient with a function
                    need_deriv_matrix_for_fv = true;
                end
            end
            need_derivs = true;
        end
    end
    
    # To make things easier, separate the terms and work with them separately
    terms = process_terms(symex);
    
    log_entry("Terms being sent to the code generator:\n\t\t" * string(terms), 3)
    # println(symex)
    # println(" -> ")
    # println(terms)
    
    ###### Below are target specific code generation functions ############################################
    if language == JULIA || language == 0
        if solver == CG
            handle_input_args_fun = handle_input_args_cg_julia;
            build_derivative_matrices_fun = build_derivative_matrices_cg_julia;
            prepare_needed_values_fun = prepare_needed_values_cg_julia;
            make_elemental_computation_fun = make_elemental_computation_cg_julia;
            
        elseif solver == DG
            handle_input_args_fun = handle_input_args_dg_julia;
            build_derivative_matrices_fun = build_derivative_matrices_dg_julia;
            prepare_needed_values_fun = prepare_needed_values_dg_julia;
            make_elemental_computation_fun = make_elemental_computation_dg_julia;
            
        elseif solver == FV
            handle_input_args_fun = handle_input_args_fv_julia;
            build_derivative_matrices_fun = build_derivative_matrices_fv_julia;
            prepare_needed_values_fun = prepare_needed_values_fv_julia;
            make_elemental_computation_fun = make_elemental_computation_fv_julia;
        end
        
        ### Build the code as a string ###
        
        # The whole body of code is stored in this string.
        code = "";
        
        # Handle input arguments
        code *= handle_input_args_fun(lorr, vors);
        code *= "\n";
        
        # If needed, compute derivative matrices
        if need_derivs
            if config.solver_type == FV
                code *= build_derivative_matrices_fun(lorr, vors, need_deriv_matrix_for_fv);
            else
                code *= build_derivative_matrices_fun(lorr, vors);
            end
            code *= "\n";
        end
        
        # Evaluate or fetch the values for each needed entity.
        code *= prepare_needed_values_fun(entities, var, lorr, vors);
        code *= "\n";
        
        # Form the final elemental calculation
        code *= make_elemental_computation_fun(terms, var, dofsper, offset_ind, lorr, vors);
        
        # For Julia, return both the generated string and a parsed Expr block of code.
        return (code, code_string_to_expr(code));
        
    ### External targets ##############################################################
    else
        code_string = external_generate_code_layer_function(var, entities, terms, lorr, vors);
        return (code_string, code_string);
    end
    
    
end
