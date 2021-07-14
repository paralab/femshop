#=
Code generation functions for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_cg_julia(lorr, vors)
    code = "";
    if vors == "volume"
        code *=
"var =   args[1]; # unknown variables
nodex =     args[2]; # global coords of element's nodes
gbl =   args[3]; # global indices of the nodes
refel = args[4]; # reference element
wdetj = args[5]; # quadrature weights scaled by detJ
J =     args[6]; # Jacobian
time =  args[7]; # time for time dependent coefficients
dt =    args[8]; # dt for time dependent problems
"
        # A trick for uniform grids to avoid repeated work
        if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
            code *= "stiffness = args[9]; # set of stiffness matrices for each dimension\n";
            code *= "mass =      args[10]; # mass matrix\n";
        end
        
    else
        # surface
        # TODO see DG
    end
    
    return code;
end

# If needed, build derivative matrices
function build_derivative_matrices_cg_julia(lorr, vors)
    code = 
"
# Note on derivative matrices:
# RQn are vandermond matrices for the derivatives of the basis functions
# with Jacobian factors. They are made like this.
# |RQ1|   | rx sx tx || Qx |
# |RQ2| = | ry sy ty || Qy |
# |RQ3|   | rz sz tz || Qz |

"
    if config.dimension == 1
        if vors == "volume"
            code *= "(RQ1, RD1) = build_deriv_matrix(refel, J);\n";
            code *= "TRQ1 = RQ1';\n"
        else
            # TODO see DG
        end
        
    elseif config.dimension == 2
        if vors == "volume"
            code *= "(RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);\n";
            code *= "(TRQ1, TRQ2) = (RQ1', RQ2');\n"
        else
            # TODO see DG
        end
        
    elseif config.dimension == 3
        if vors == "volume"
            code *= "(RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);\n";
            code *= "(TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');\n"
        else
            # TODO see DG
        end
    end
    
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_cg_julia(entities, var, lorr, vors)
    code = "";
    for i=1:length(entities)
        cname = make_coef_name(entities[i]);
        if is_test_function(entities[i])
            # Assign it a transpose quadrature matrix
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= cname * " = TRQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of test function\n";
                end
            else
                code *= cname * " = refel.Q'; # test function.\n";
            end
        elseif is_unknown_var(entities[i], var) && lorr == LHS
            if length(entities[i].derivs) > 0
                xyzchar = ["x","y","z"];
                for di=1:length(entities[i].derivs)
                    code *= cname * " = RQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
                end
            else
                code *= cname * " = refel.Q; # trial function.\n";
            end
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt
            elseif ctype == 0
                # It was a number, do nothing?
            elseif ctype == 1 # a constant wrapped in a coefficient
                # This generates something like: coef_k_i = 4;
                if length(entities[i].derivs) > 0
                    code *= cname * " = 0; # NOTE: derivative applied to constant coefficient = 0\n";
                else
                    code *= cname * " = " * string(cval) * ";\n";
                end
                
            elseif ctype == 2 # a coefficient function
                # This generates something like:
                ######################################
                # coef_n_i = zeros(refel.Np);
                # for coefi = 1:refel.Np
                #     coef_k_i[coefi] = (Femshop.genfunctions[cval]).func(x[1,coefi], x[2,coefi],x[3,coefi],time);
                # end
                ######################################
                cargs = "(nodex[coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs = "(nodex[1, coefi], nodex[2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs = "(nodex[1, coefi], nodex[2, coefi], nodex[3, coefi], time)";
                end
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np);\n";
                    code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Femshop.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*" and interpolate at quadrature points.\n";
                        end
                    else
                        code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                else
                    #TODO surface
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = copy((Femshop.variables[1]).values[1, gbl])
                if vors == "volume"
                    code *= cname * " = copy((Femshop.variables["*string(cval)*"]).values["*string(entities[i].index)*", gbl]);\n";
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RQ"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*"\n";
                        end
                    else
                        code *= cname * " = refel.Q * " * cname * "; # Interpolate at quadrature points.\n";
                    end
                else
                    #TODO surface
                end
                
            end
        end # if coefficient
    end # entity loop
    
    return code;
end

function make_elemental_computation_cg_julia(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    code = "";
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "element_vector = zeros(refel.Np * "*string(dofsper)*"); # Allocate the returned vector.\n"
        else
            code *= "element_matrix = zeros(refel.Np * "*string(dofsper)*", refel.Np * "*string(dofsper)*"); # Allocate the returned matrix.\n"
        end
    end
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        # Submatrices or subvectors for each component
        if lorr == LHS
            submatrices = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            submatrices = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(submatrices)
            submatrices[smi] = "";
        end
        
        if typeof(var) <: Array
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[vi][ci][i], var, lorr);
                        
                        # println(terms)
                        # println(terms[vi])
                        # println(terms[vi][ci])
                        # println(terms[vi][ci][i])
                        # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
                        # Find the appropriate submatrix for this term
                        submati = offset_ind[vi] + test_ind;
                        submatj = trial_ind;
                        if lorr == LHS
                            submat_ind = submati + dofsper * (submatj-1);
                        else
                            submat_ind = submati;
                        end
                        
                        
                        if length(submatrices[submat_ind]) > 1
                            submatrices[submat_ind] *= " .+ " * term_result;
                        else
                            submatrices[submat_ind] = term_result;
                        end
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[ci][i], var, lorr);
                    
                    # Find the appropriate submatrix for this term
                    if lorr == LHS
                        submat_ind = test_ind + dofsper * (trial_ind-1);
                    else
                        submat_ind = test_ind;
                    end
                    
                    if length(submatrices[submat_ind]) > 1
                        submatrices[submat_ind] *= " .+ " * term_result;
                    else
                        submatrices[submat_ind] = term_result;
                    end
                end
            end
            
        end
        
        # Put the submatrices together into element_matrix or element_vector
        if lorr == LHS
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(submatrices[emi, emj]) > 1
                        rangei = "(("*string(emi)*"-1)*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        rangej = "(("*string(emj)*"-1)*refel.Np + 1):("*string(emj)*"*refel.Np)";
                        code *= "element_matrix["*rangei*", "*rangej*"] = " * submatrices[emi,emj] * "\n";
                    end
                end
            end
            code *= "return element_matrix;\n"
            
        else # RHS
            for emi=1:dofsper
                if length(submatrices[emi]) > 1
                    rangei = "(("*string(emi)*"-1)*refel.Np + 1):("*string(emi)*"*refel.Np)";
                    code *= "element_vector["*rangei*"] = " * submatrices[emi] * "\n";
                end
            end
        end
        code *= "return element_vector;\n"
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            result = "zeros(refel.Np, refel.Np)";
        else
            result = "zeros(refel.Np)";
        end
        
        #process each term
        for i=1:length(terms)
            (term_result, test_ind, trial_ind) = generate_term_calculation_cg_julia(terms[i], var, lorr);
            
            if i > 1
                result *= " .+ " * term_result;
            else
                result = term_result;
            end
        end
        code *= "return " * result * ";\n";
    end
    
    return code;
end

function generate_term_calculation_cg_julia(term, var, lorr)
    result = "";
    
    if lorr == LHS
        (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term, var);
        # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                    string(replace_entities_with_symbols(coef_part)) * ") * " * 
                    string(replace_entities_with_symbols(trial_part));
        else # no coef_part
            result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                    string(replace_entities_with_symbols(trial_part));
        end
    else
        (test_part, trial_part, coef_part, test_ind, trial_ind) = separate_factors(term);
        # RHS: test_part * (weight_part .* coef_part)
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(test_part)) * " * (wdetj .* " * 
                    string(replace_entities_with_symbols(coef_part)) * ")";
        else
            result = string(replace_entities_with_symbols(test_part)) * " * (wdetj)";
        end
    end
    
    return (result, test_ind, trial_ind);
end
