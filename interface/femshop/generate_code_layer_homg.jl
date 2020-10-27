
###############################################################################################################
# homg
###############################################################################################################

# homg version returns a string that will be included in the Linear of Bilinear file
function generate_code_layer_homg(symex, var, lorr)
    code = ""; # the code will be in this string
    
    need_derivative = false;
    term_derivative = [];
    needed_coef = [];
    needed_coef_ind = [];
    test_ind = [];
    trial_ind = [];
    
    # symex is an array of arrays of symengine terms
    # turn each one into an Expr
    terms = [];
    sz = size(symex);
    if length(sz) == 1 # scalar or vector
        for i=1:length(symex)
            for ti=1:length(symex[i])
                push!(terms, Meta.parse(string(symex[i][ti])));
            end
        end
    elseif length(sz) == 2 # matrix
        #TODO
    end
    
    # Process the terms turning them into the code layer
    # cindex = 0;
    # for i=1:length(terms)
    #     (codeterm, der, coe) = process_term_dendro(terms[i], cindex, length(terms), i, var, lorr);
    #     need_derivative = need_derivative || der;
    #     push!(term_derivative, der);
    #     append!(needed_coef, coe);
    #     cindex = length(needed_coef);
    #     push!(code_terms, codeterm);
    # end
    # Process the terms turning them into the code layer
    code_terms = [];
    for i=1:length(terms)
        (codeterm, der, coe, coeind, testi, trialj) = process_term_homg(terms[i], var, lorr);
        if coeind == -1
            # processing failed due to nonlinear term
            printerr("term processing failed for: "*string(terms[i]));
            return nothing;
        end
        need_derivative = need_derivative || der;
        append!(needed_coef, coe);
        append!(needed_coef_ind, coeind);
        # change indices into one number
        
        push!(test_ind, testi);
        push!(trial_ind, trialj);
        push!(code_terms, codeterm);
    end
    
    # If coefficients need to be computed, do so
    # # First remove duplicates
    unique_coef = [];
    unique_coef_ind = [];
    for i=1:length(needed_coef)
        already = false;
        for j=1:length(unique_coef)
            if unique_coef[j] === needed_coef[i] && unique_coef_ind[j] == needed_coef_ind[i]
                already = true;
            end
        end
        if !already
            push!(unique_coef, needed_coef[i]);
            push!(unique_coef_ind, needed_coef_ind[i]);
        end
    end
    needed_coef = unique_coef;
    needed_coef_ind = unique_coef_ind;
    
    # For constant coefficients, this generates something like:
    ######################################
    # double coef_n_i = a_i;
    ######################################
    
    # For variable coefficients, this generates something like:
    ######################################
    # double* coef_n_i = new double[nPe];
    # double val;
    # for(unsigned int i=0; i<nPe; i++){
    #     a_i(coords[i*3+0], coords[i*3+0], coords[i*3+0], &val);
    #     coef_n_i[coefi] = val;
    # }
    ######################################
    coef_alloc = "";
    coef_loop = "";
    if length(needed_coef) > 0
        coef_loop = "gpts = mesh.element_gauss(e, refel);\n";
        for i=1:length(needed_coef)
            if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                cind = get_coef_index(needed_coef[i]);
                if cind < 0
                    # probably a variable
                    cind = string(needed_coef[i]);
                    # This presents a problem in dendro. TODO
                    printerr("HOMG not yet available for multivariate problems. Expect an error.");
                end
                # The string name for this coefficient
                cname = "coef_"*string(cind)*"_"*string(needed_coef_ind[i]);
                
                (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
                
                if ctype == 1
                    # constant coefficient -> coef_n_i = cval
                    coef_alloc *= cname*" = "*string(cval)*";\n";
                    
                elseif ctype == 2
                    # genfunction coefficients -> coef_n_i = coef(idx)
                    coef_loop *= cname*" = "*cval*"(idx);\n";
                    
                elseif ctype == 3
                    # variable values -> coef_n = variable.values
                    #TODO THIS IS AN ERROR. multivariate support needed.
                end
                
            end
        end
    end
    
    dofsper = 1;
    if typeof(var) <: Array
        for vi=1:length(var)
            dofsper += length(var[vi].symvar.vals); # The number of components for this variable
        end
    elseif !(var.type == SCALAR)
        dofsper = length(var.symvar.vals);
    end
    # Not ready
    if dofsper > 1
        printerr("homg not ready for multi dofs per node. Expect errors");
    end
    
    # Put the pieces together
    code = coef_loop*"\n";
    if lorr == LHS
        code *= "elMat = "*code_terms[1];
        for i=2:length(code_terms)
            code *= " + "*code_terms[i];
        end
        code *= ";\n";
    else
        code *= "elVec = "*code_terms[1];
        for i=2:length(code_terms)
            code *= " + "*code_terms[i];
        end
        code *= ";\n";
    end
    
    return code;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term_homg(sterm, var, lorr)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    needed_coef_ind = [];
    
    test_part = "";
    trial_part = "";
    coef_part = "";
    test_component = 0;
    trial_component = 0;
    deriv_dir = 0;
    
    # extract each of the factors.
    factors = separate_factors(term);
    
    # strip off all negatives, combine and reattach at the end
    neg = false;
    for i=1:length(factors)
        if typeof(factors[i]) == Expr && factors[i].args[1] === :- && length(factors[i].args) == 2
            neg = !neg;
            factors[i] = factors[i].args[2];
        end
    end
    
    # Separate factors into test/trial/coefficient parts
    coef_facs = [];
    coef_inds = [];
    for i=1:length(factors)
        (index, v, mods) = extract_symbols(factors[i]);
        
        if is_test_func(v)
            test_component = index; # the vector index
            if length(mods) > 0
                # TODO more than one derivative mod
                need_derivative = true;
                deriv_dir = parse(Int, mods[1][2]);
                if deriv_dir == 1
                    if config.dimension == 2
                        test_part = "([diag(Jac.rx) diag(Jac.sx)] * [refel.Qx; refel.Qy])\'";
                    else
                        test_part = "([diag(Jac.rx) diag(Jac.sx) diag(Jac.tx)] * [refel.Qx; refel.Qy; refel.Qz])\'";
                    end
                elseif deriv_dir == 2
                    if config.dimension == 2
                        test_part = "([diag(Jac.ry) diag(Jac.sy)] * [refel.Qx; refel.Qy])\'";
                    else
                        test_part = "([diag(Jac.ry) diag(Jac.sy) diag(Jac.ty)] * [refel.Qx; refel.Qy; refel.Qz])\'";
                    end
                elseif deriv_dir == 3
                    if config.dimension == 2
                        test_part = "([diag(Jac.rz) diag(Jac.sz)] * [refel.Qx; refel.Qy])\'";
                    else
                        test_part = "([diag(Jac.rz) diag(Jac.sz) diag(Jac.tz)] * [refel.Qx; refel.Qy; refel.Qz])\'";
                    end
                elseif deriv_dir == 4
                    # This will eventually be a time derivative
                    printerr("Derivative index problem in "*string(factors[i]));
                else
                    printerr("Derivative index problem in "*string(factors[i]));
                end
            else
                # no derivative mods
                test_part = "refel.Q\'";
            end
        elseif is_unknown_var(v, var)
            if !(trial_part == "")
                # Two unknowns multiplied in this term. Nonlinear. abort.
                printerr("Nonlinear term. Code layer incomplete.");
                return (-1, -1, -1, -1, -1, -1);
            end
            trial_component = index;
            if length(mods) > 0
                # TODO more than one derivative mod
                need_derivative = true;
                deriv_dir = parse(Int, mods[1][2]);
                if deriv_dir == 1
                    if config.dimension == 2
                        trial_part = "[diag(Jac.rx) diag(Jac.sx)] * [refel.Qx; refel.Qy]";
                    else
                        trial_part = "[diag(Jac.rx) diag(Jac.sx) diag(Jac.tx)] * [refel.Qx; refel.Qy; refel.Qz]";
                    end
                elseif deriv_dir == 2
                    if config.dimension == 2
                        trial_part = "[diag(Jac.ry) diag(Jac.sy)] * [refel.Qx; refel.Qy]";
                    else
                        trial_part = "[diag(Jac.ry) diag(Jac.sy) diag(Jac.ty)] * [refel.Qx; refel.Qy; refel.Qz]";
                    end
                elseif deriv_dir == 3
                    if config.dimension == 2
                        trial_part = "[diag(Jac.rz) diag(Jac.sz)] * [refel.Qx; refel.Qy]";
                    else
                        trial_part = "[diag(Jac.rz) diag(Jac.sz) diag(Jac.tz)] * [refel.Qx; refel.Qy; refel.Qz]";
                    end
                elseif deriv_dir == 4
                    # This will eventually be a time derivative
                    printerr("Derivative index problem in "*string(factors[i]));
                else
                    printerr("Derivative index problem in "*string(factors[i]));
                end
            else
                # no derivative mods
                trial_part = "refel.Q";
            end
        else
            if length(index) == 1
                ind = index[1];
            end
            push!(coef_facs, v);
            push!(coef_inds, ind);
        end
    end
    
    # If rhs, change var into var.values and treat as a coefficient
    if lorr == RHS && !(trial_part === nothing)
        # tmpv = :(a.values[gbl]);
        # tmpv.args[1].args[1] = var.symbol; #TODO, will not work for var=array
        # push!(coef_facs, tmpv);
        # push!(coef_inds, trial_component);
    end
    
    # If there was no trial part, it's an RHS and we need to finish the quadrature with this
    if trial_part == ""
        trial_part = "refel.Q";
    end
    
    # build weight/coefficient parts
    weight_part = "refel.W .* detJ";
    
    # If term is negative, apply it here
    if neg
        weight_part = "-"*weight_part;
    end
    
    # coefficients
    if length(coef_facs) > 0
        for j=1:length(coef_facs)
            tmp = coef_facs[j];
            if typeof(tmp) == Symbol && !(tmp ===:dt)
                push!(needed_coef, tmp);
                push!(needed_coef_ind, coef_inds[j]);
                ind = get_coef_index(coef_facs[j]);
                if ind >= 0
                    tmp = "coef_"*string(ind)*"_"*string(coef_inds[j]);
                else
                    tmp = "coef_"*string(tmp)*"_"*string(coef_inds[j]);
                end
            else
                tmp = string(tmp);
            end
            if j>1
                coef_part = coef_part*".*"*tmp;
            else
                coef_part = coef_part*tmp;
            end
        end
    end
    
    # Put the pieces togetherd
    if !(coef_part === "")
        if lorr == LHS
            weight_part = weight_part*" .* "*coef_part;
            weight_part = "diag("*weight_part*")";
            term = test_part * " * " *  weight_part * " * " *  trial_part;
        else # RHS
            weight_part = "diag("*weight_part*")";
            term = test_part * " * " *  weight_part * " * (" *  trial_part *" * "* coef_part *")";
        end
    else
        weight_part = "diag("*weight_part*")";
        term = test_part * " * " *  weight_part * " * " *  trial_part;
    end
    
    return (term, need_derivative, needed_coef, needed_coef_ind, test_component, trial_component);
end
