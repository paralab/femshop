#=
Use the symbolic layer expressions to generate the FEM code
=#

# lorr = LHS or RHS, vors = volume or surface
function generate_code_layer_julia(symex, var, lorr, vors)
    # The whole body of code is stored in this string.
    code = "";
    
    # Extract the args. This must match the arguments passed by the solver.
    if vors == "volume"
        code *=
"var =   args[1]; # unknown variables
x =     args[2]; # global coords of element's nodes
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
        
    end
    
    # Count variables, dofs, and store offsets
    multivar = typeof(var) <:Array;
    varcount = 1;
    dofsper = 0;
    if multivar
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        dofsper = length(var[1].symvar.vals);
        for i=2:length(var)
            offset_ind[i] = dofsper;
            dofsper = dofsper + length(var[i].symvar.vals);
        end
    else
        dofsper = length(var.symvar.vals);
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
    for i=1:length(entities)
        if length(entities[i].derivs) > 0
            need_derivs = true;
        end
    end
    # If needed, compute derivative matrices
    if need_derivs
        code *= 
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
                if lorr == LHS
                    code *= "TRQ1 = RQ1';\n"
                end
            else
                #TODO
            end
            
        elseif config.dimension == 2
            if vors == "volume"
                code *= "(RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);\n";
                if lorr == LHS
                    code *= "(TRQ1, TRQ2) = (RQ1', RQ2');\n"
                end
            else
                #TODO
            end
            
        elseif config.dimension == 3
            if vors == "volume"
                code *= "(RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);\n";
                if lorr == LHS
                    code *= "(TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');\n"
                end
            else
                #TODO
            end
        end
    end
    code *= "\n";
    
    # Then evaluate or fetch the values for each needed entity.
    # Note that test functions and unknowns on the LHS do not need this.
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
                cargs = "(x[coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs = "(x[1, coefi], x[2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs = "(x[1, coefi], x[2, coefi], x[3, coefi], time)";
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
    code *= "\n";
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "element_vector = zeros(refel.Np * "*string(dofsper)*"); # Allocate the returned vector.\n"
        else
            code *= "element_matrix = zeros(refel.Np * "*string(dofsper)*", refel.Np * "*string(dofsper)*"); # Allocate the returned matrix.\n"
        end
    end
    
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term is multiplied by one test function component, and if LHS, involves one unknown component.
    # The submatrix modified by a term is determined by these, so go through the terms and divide them
    # into their submatrix expressions. 
    # Each term will look something like 
    # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
    # RHS: test_part * (weight_part .* coef_part)
    
    # To make things easier, separate the terms and work with them separately
    terms = process_terms(symex);
    
    # print(symex)
    # print(" -> ")
    # println(terms)
    
    # Separate the factors of each term into test, trial, coef and form the calculation
    if dofsper > 1
        if lorr == LHS
            submatrices = Array{String, 2}(undef, dofsper, dofsper);
            for smi=1:length(submatrices)
                submatrices[smi] = "";
            end
            
        else # RHS
            submatrices = Array{String, 1}(undef, dofsper);
            for smi=1:length(submatrices)
                submatrices[smi] = "";
            end
        end
    else # one dof
        terms = terms[1];
        if lorr == LHS
            result = "";
            #process each term
            for i=1:length(terms)
                # find test and trial functions
                (test_part, trial_part, coef_part) = separate_factors(terms[i], var);
                
                if i > 1
                    result *= " .+ ";
                end
                
                # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
                if !(coef_part === nothing)
                    result *= string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
                            string(replace_entities_with_symbols(coef_part)) * ") * " * 
                            string(replace_entities_with_symbols(trial_part));
                else # no coef_part
                    result *= string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
                            string(replace_entities_with_symbols(trial_part));
                end
            end
            code *= "return " * result * ";\n";
            
        else # RHS
            result = "";
            #process each term
            for i=1:length(terms)
                # find test function
                #println(string(i)*" : "*string(terms))
                (test_part, trial_part, coef_part) = separate_factors(terms[i]);
                
                if i > 1
                    result *= " .+ ";
                end
                
                # RHS: test_part * (weight_part .* coef_part)
                if !(coef_part === nothing)
                    result *= string(replace_entities_with_symbols(test_part)) * " * (wdetj .* " * 
                            string(replace_entities_with_symbols(coef_part)) * ")";
                else
                    result *= string(replace_entities_with_symbols(test_part)) * " * (wdetj)";
                end
            end
            code *= "return " * result * ";\n";
        end
    end
    
    return code;
end
