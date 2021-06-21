#=
Use the symbolic layer expressions to generate the FEM code
=#

# External code gen in these similar files
include("generate_code_layer_dendro.jl");
include("generate_code_layer_homg.jl");
include("generate_code_layer_matlab.jl");
include("generate_code_layer_cachesim.jl");

# Surface integrals are generated separately
include("generate_code_layer_surface.jl");

function generate_code_layer(ex, var, lorr)
    if config.solver_type == FV
        return generate_code_layer_fv(ex, var, lorr)
    end
    
    if use_cachesim
        if language == 0 || language == JULIA
            return generate_code_layer_cachesim(ex, var, lorr);
        else
            printerr("Using cachesim. Not generating code to solve.")
        end
    else
        if language == 0 || language == JULIA
            return generate_code_layer_julia(ex, var, lorr);
        elseif language == CPP
            if framework == DENDRO
                return generate_code_layer_dendro(ex, var, lorr);
            else
                printerr("Plain C++ is not ready for code layer gen.")
            end
            
        elseif language == MATLAB
            if framework == HOMG
                return generate_code_layer_homg(ex, var, lorr);
            else
                return generate_code_layer_matlab(ex, var, lorr);
            end
            
        elseif framework == CUSTOM_GEN_TARGET
            return custom_code_layer_fun(ex, var, lorr);
        end
    end
end

# Julia and utils are in this file

###############################################################################################################
# julia
###############################################################################################################

# Julia version returns an expression for the generated function for linear or bilinear term
function generate_code_layer_julia(symex, var, lorr)
    # This is the basic info passed in "args"
    code = Expr(:block);
    push!(code.args, :(var = args[1]));  # unknown variables
    push!(code.args, :(x = args[2]));    # global coords of element's nodes
    push!(code.args, :(gbl = args[3]));  # global indices of the nodes
    push!(code.args, :(refel = args[4]));# reference element
    push!(code.args, :(wdetj = args[5]));# weights * detJ
    push!(code.args, :(J = args[6]));    # Jacobian
    push!(code.args, :(borl = args[7])); # bilinear or linear? lhs or rhs?
    push!(code.args, :(time = args[8])); # time for time dependent coefficients
    if prob.time_dependent
        push!(code.args, :(dt = args[9])); # dt for time dependent problems
    end
    
    # A trick for uniform grids to avoid repeated work
    if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
        push!(code.args, :(stiffness = args[10])); # dt for time dependent problems
        push!(code.args, :(mass = args[11])); # dt for time dependent problems
    end
    
    # push!(code.args, :((detJ, J) = geometric_factors(refel, x)));
    # push!(code.args, :(wdetj = refel.wg .* detJ));
    
    need_derivative = false;
    needed_coef = [];
    needed_coef_name = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    test_ind = [];
    trial_ind = [];
    
    # For multi variables, count them
    multivar = typeof(var) <:Array;
    varcount = 1;
    if multivar
        varcount = length(var);
        offset_ind = zeros(Int, varcount);
        tmp = length(var[1].symvar.vals);
        for i=2:length(var)
            offset_ind[i] = tmp;
            tmp = tmp + length(var[i].symvar.vals);
        end
    end
    
    # symex is an array of arrays of symengine terms, or array of arrays of arrays for multivar
    # turn each one into an Expr for translation purposes
    terms = [];
    if multivar
        for vi=1:varcount
            push!(terms, terms_to_expr(symex[vi]));
        end
    else
        terms = terms_to_expr(symex);
    end
    
    # Process the terms turning them into the code layer
    if multivar
        for vi=1:varcount
            subtest_ind = [];
            subtrial_ind = [];
            for i=1:length(terms[vi])
                (codeterm, der, coe, coenam, coeind, coederiv, testi, trialj) = process_term_julia(terms[vi][i], var, lorr, offset_ind);
                if coeind == -1
                    # processing failed due to nonlinear term
                    printerr("term processing failed for: "*string(terms[vi][i])*" , possible nonlinear term?");
                    return nothing;
                end
                need_derivative = need_derivative || der;
                append!(needed_coef, coe);
                append!(needed_coef_name, coenam);
                append!(needed_coef_ind, coeind);
                append!(needed_coef_deriv, coederiv);
                # change indices into one number
                
                push!(subtest_ind, testi);
                push!(subtrial_ind, trialj);
                terms[vi][i] = codeterm;
            end
            push!(test_ind, subtest_ind);
            push!(trial_ind, subtrial_ind);
        end
    else
        for i=1:length(terms)
            (codeterm, der, coe, coenam, coeind, coederiv, testi, trialj) = process_term_julia(terms[i], var, lorr);
            if coeind == -1
                # processing failed due to nonlinear term
                printerr("term processing failed for: "*string(terms[i])*" , possible nonlinear term");
                return nothing;
            end
            need_derivative = need_derivative || der;
            append!(needed_coef, coe);
            append!(needed_coef_name, coenam);
            append!(needed_coef_ind, coeind);
            append!(needed_coef_deriv, coederiv);
            # change indices into one number
            
            push!(test_ind, testi);
            push!(trial_ind, trialj);
            terms[i] = codeterm;
        end
    end
    
    # If derivatives are needed, prepare the appropriate matrices
    if need_derivative
        if config.dimension == 1
            push!(code.args, :((RQ1,RD1) = build_deriv_matrix(refel, J)));
            push!(code.args, :(TRQ1 = RQ1'));
        elseif config.dimension == 2
            push!(code.args, :((RQ1,RQ2,RD1,RD2) = build_deriv_matrix(refel, J)));
            push!(code.args, :((TRQ1,TRQ2) = (RQ1',RQ2')));
        elseif config.dimension == 3
            push!(code.args, :((RQ1,RQ2,RQ3,RD1,RD2,RD3) = build_deriv_matrix(refel, J)));
            push!(code.args, :((TRQ1,TRQ2,TRQ3) = (RQ1',RQ2',RQ3')));
        end
    end
    
    # If coefficients need to be computed, do so
    # # First remove duplicates
    #println("coef-"*string(length(needed_coef))*": "*string(needed_coef));
    #println("ind-"*string(length(needed_coef_ind))*": "*string(needed_coef_ind));
    #println("der-"*string(length(needed_coef_deriv))*": "*string(needed_coef_deriv));
    unique_coef = [];
    unique_coef_name = [];
    unique_coef_ind = [];
    unique_coef_deriv = [];
    for i=1:length(needed_coef)
        already = false;
        for j=1:length(unique_coef)
            if unique_coef[j] === needed_coef[i] && unique_coef_ind[j] == needed_coef_ind[i] && unique_coef_deriv[j] == needed_coef_deriv[i] && unique_coef_name[j] == needed_coef_name[i]
                already = true;
            end
        end
        if !already
            #println("coef: "*string(needed_coef[i])*" , ind: "*string(needed_coef_ind[i])*" , deriv: "*string(needed_coef_deriv[i]));
            push!(unique_coef, needed_coef[i]);
            push!(unique_coef_name, needed_coef_name[i]);
            push!(unique_coef_ind, needed_coef_ind[i]);
            push!(unique_coef_deriv, needed_coef_deriv[i]);
        end
    end
    needed_coef = unique_coef;
    needed_coef_name = unique_coef_name;
    needed_coef_ind = unique_coef_ind;
    needed_coef_deriv = unique_coef_deriv;
    # println("coef-"*string(length(needed_coef))*": "*string(needed_coef));
    # println("ind-"*string(length(needed_coef_ind))*": "*string(needed_coef_ind));
    # println("der-"*string(length(needed_coef_deriv))*": "*string(needed_coef_deriv));
    # println("nam-"*string(length(needed_coef_name))*": "*string(needed_coef_name));
    
    # For constant coefficients, this generates something like:
    ######################################
    # coef_n_i = a.value[i];
    ######################################
    
    # For variable coefficients, this generates something like:
    ######################################
    # coef_n_i = zeros(refel.Nqp);
    # for coefi = 1:refel.Nqp
    #     coef_n_i[coefi] = a.value[i].func(x[1,coefi], x[2,coefi],x[3,coefi],time);
    # end
    ######################################
    got_coef_quadrature_points = false; # compute quadrature point coordinates only if needed
    coef_at_nodes = zeros(Bool,length(needed_coef)); # set to true if the coef values are known at nodes(not quadrature points)
    
    if length(needed_coef) > 0
        cloop = :(for coefi=1:refel.Np end);
        cloopin = Expr(:block);
        cargs = [:(x[coefi]); 0; 0; :time];
        if config.dimension == 2
            cargs = [:(x[1,coefi]); :(x[2,coefi]); 0; :time];
        elseif config.dimension == 3
            cargs = [:(x[1,coefi]); :(x[2,coefi]); :(x[3,coefi]); :time];
        end
        
        for i=1:length(needed_coef)
            if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                tmps = needed_coef_name[i];
                tmpc = Symbol(tmps);
                
                (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
                if ctype == 1
                    # constant coefficient -> coef_n = cval
                    tmpn = cval;
                    push!(code.args, Expr(:(=), tmpc, tmpn));
                elseif ctype == 2
                    coef_at_nodes[i] = true;
                    # genfunction coefficients -> coef_n_i = coef.value[i].func(cargs)
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc;
                    tmpn = :(Femshop.genfunctions[$cval]); # Femshop.genfunctions[cval]
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    
                    push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                elseif ctype == 3
                    coef_at_nodes[i] = true;
                    
                    # variable values -> coef_n = variable.values
                    if variables[cval].type == SCALAR
                        tmpb = :(copy(Femshop.variables[$cval].values[gbl])); 
                    else
                        compo = needed_coef_ind[i];
                        tmpb = :(copy(Femshop.variables[$cval].values[$compo, gbl]));
                    end
                    
                    push!(code.args, Expr(:(=), tmpc, tmpb));
                end
                
            end# number?
        end# coef loop
        
        # Write the loop that computes coefficient values
        if length(cloopin.args) > 0
            cloop.args[2] = cloopin;
            push!(code.args, cloop); # add loop to code
        end
        
        # Apply derivatives after the initializing loop
        # Will look like: coef_i_j = RDn * coef_i_j
        for i=1:length(needed_coef_deriv)
            if length(needed_coef_deriv[i][2]) > 0 && !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                tmps = needed_coef_name[i];
                tmpc = Symbol(tmps);
                
                dmat = Symbol("RD"*needed_coef_deriv[i][3]);
                tmpb= :(length($tmpc) == 1 ? 0 : $dmat * $tmpc);
                
                push!(code.args, Expr(:(=), tmpc, tmpb));
            else
                # derivatives of constant coefficients should be zero
                # TODO
            end
        end
        
        # Interpolate the nodal coefficients at quadrature points
        if lorr == LHS
            for i=1:length(needed_coef)
                if coef_at_nodes[i] == true
                    tmps = needed_coef_name[i];
                    tmpc = Symbol(tmps);
                    tmpb = :(refel.Q * $tmpc);
                    push!(code.args, Expr(:(=), tmpc, tmpb));
                end
            end
        end
        
    end# needed_coef loop
    
    # finally add the code expression
    # For multiple dofs per node it will be like:
    # [A11  A12  A13 ...] where the indices are from test_ind and trial_ind (vector components)
    # [A21  A22  A23 ...] Note: things will be rearranged when inserted into the global matrix/vector
    # [...              ]
    # [...              ]
    #
    # [b1 ] from test_ind
    # [b2 ]
    # [...]
    
    # Allocate if needed
    dofsper = 0;
    if typeof(var) <: Array
        for vi=1:length(var)
            dofsper = dofsper + length(var[vi].symvar.vals); # The number of components for this variable
        end
    else
        dofsper = length(var.symvar.vals);
    end
    
    if dofsper > 1
        if lorr == RHS
            push!(code.args, Expr(:(=), :element_vector, :(zeros(refel.Np*$dofsper)))); # allocate vector
        else
            push!(code.args, Expr(:(=), :element_matrix, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper)))); # allocate matrix
        end
    end
    
    result = nothing; # Will hold the returned expression
    if typeof(var) <: Array # multivar
        for vi=1:length(var)
            # Add terms into full matrix according to testind/trialind
            # Each component/dof should have one expression so that submatrix is only modified once.
            if lorr == LHS
                #comps = length(var.symvar.vals);
                comps = dofsper;
                submatrices = Array{Any,2}(undef, comps, comps);
                for smi=1:length(submatrices)
                    submatrices[smi] = nothing;
                end
                for i=1:length(terms[vi])
                    ti = test_ind[vi][i][1] + offset_ind[vi];
                    tj = trial_ind[vi][i][1];
                    
                    if submatrices[ti, tj] === nothing
                        submatrices[ti, tj] = terms[vi][i];
                    else
                        addexpr = :(a+b);
                        addexpr.args[2] = submatrices[ti, tj];
                        addexpr.args[3] = terms[vi][i];
                        submatrices[ti, tj] = addexpr;
                    end
                end
                
                for cj=1:comps
                    for ci=1:comps
                        if !(submatrices[ci, cj] === nothing)
                            ti = ci-1;
                            tj = cj-1;
                            sti = :($ti*refel.Np + 1);
                            eni = :(($ti + 1)*refel.Np);
                            stj = :($tj*refel.Np + 1);
                            enj = :(($tj + 1)*refel.Np);
                            
                            push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                        end
                    end
                end
                
                result = :element_matrix;
                
            else #RHS
                #comps = length(var.symvar.vals);
                comps = dofsper;
                submatrices = Array{Any,1}(undef, comps);
                for smi=1:length(submatrices)
                    submatrices[smi] = nothing;
                end
                for i=1:length(terms[vi])
                    ti = test_ind[vi][i][1] + offset_ind[vi];
                    
                    if submatrices[ti] === nothing
                        submatrices[ti] = terms[vi][i];
                    else
                        addexpr = :(a+b);
                        addexpr.args[2] = submatrices[ti];
                        addexpr.args[3] = terms[vi][i];
                        submatrices[ti] = addexpr;
                    end
                end
                
                for ci=1:comps
                    if !(submatrices[ci] === nothing)
                        ti = ci-1;
                        sti = :($ti*refel.Np + 1);
                        eni = :(($ti + 1)*refel.Np);
                        
                        push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), submatrices[ci]));
                    end
                end
                
                result = :element_vector;
                
            end
        end
        
    else
        if length(terms) > 1
            if var.type == SCALAR # Only one component
                tmp = :(a.+b);
                tmp.args = [:.+];
                for i=1:length(terms)
                    push!(tmp.args, terms[i]);
                end
                result = tmp;
            else # More than one component
                # Add terms into full matrix according to testind/trialind
                # Each component/dof should have one expression so that submatrix is only modified once.
                if lorr == LHS
                    comps = length(var.symvar.vals);
                    submatrices = Array{Any,2}(undef, comps, comps);
                    for smi=1:length(submatrices)
                        submatrices[smi] = nothing;
                    end
                    for i=1:length(terms)
                        ti = test_ind[i][1];
                        tj = trial_ind[i][1];
                        
                        if submatrices[ti, tj] === nothing
                            submatrices[ti, tj] = terms[i];
                        else
                            addexpr = :(a.+b);
                            addexpr.args[2] = submatrices[ti, tj];
                            addexpr.args[3] = terms[i];
                            submatrices[ti, tj] = addexpr;
                        end
                    end
                    
                    for cj=1:comps
                        for ci=1:comps
                            if !(submatrices[ci, cj] === nothing)
                                ti = ci-1;
                                tj = cj-1;
                                sti = :($ti*refel.Np + 1);
                                eni = :(($ti + 1)*refel.Np);
                                stj = :($tj*refel.Np + 1);
                                enj = :(($tj + 1)*refel.Np);
                                
                                push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                            end
                        end
                    end
                    
                    result = :element_matrix;
                    
                else #RHS
                    comps = length(var.symvar.vals);
                    submatrices = Array{Any,1}(undef, comps);
                    for smi=1:length(submatrices)
                        submatrices[smi] = nothing;
                    end
                    for i=1:length(terms)
                        ti = test_ind[i][1];
                        
                        if submatrices[ti] === nothing
                            submatrices[ti] = terms[i];
                        else
                            addexpr = :(a.+b);
                            addexpr.args[2] = submatrices[ti];
                            addexpr.args[3] = terms[i];
                            submatrices[ti] = addexpr;
                        end
                    end
                    
                    for ci=1:comps
                        if !(submatrices[ci] === nothing)
                            ti = ci-1;
                            sti = :($ti*refel.Np + 1);
                            eni = :(($ti + 1)*refel.Np);
                            
                            push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), submatrices[ci]));
                        end
                    end
                    
                    result = :element_vector;
                    
                end
                
            end
        elseif length(terms) == 1# one term (one variable)
            result = terms[1];
            
        else # there were no terms. Just return arrays of zeros
            if lorr == LHS
                result = :(zeros(refel.Np, refel.Np));
            else
                result = :(zeros(refel.Np));
            end
        end
    end
    
    
    # At this point everything is packed into terms[1]
    push!(code.args, Expr(:return, result));
    return code;
end

#=
Changes the symbolic layer term into a code layer term.
An example term with variable u, coefficient c
Q is a quadrature matrix
wgdetj is weights and geometric factors like: (weights*detJ)

2*_c_1*_u_1*_v_1    ->  Q' * diagm(2 .* coef_0_1 .* wgdetj) * Q       on LHS
                    ->  Q' * (wgdetj .* coef_0_1 .* (Q * coef_u_1))   on RHS
                    
Also records the needed derivative matrices like: J.rx * refel.Ddr
D1 -> D1x
D2 -> D1y
D1_D1 -> D2x
etc.

Returns:
- expression for term
- needed coefficients
- coefficient indices(vector components)
- assigned names for coefficients (_c_1 -> coef_n_1 when c is the nth coefficient in coefficients array)
- needed derivative matrices

=#
function process_term_julia(sterm, var, lorr, offset_ind=0)
    term = copy(sterm);
    need_derivative = false;
    need_derivative_for_coefficient = false;
    needed_coef = [];
    needed_coef_name = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    weight_part = :wdetj;
    test_component = 0;
    trial_component = 0;
    
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
    coef_derivs = [];
    coef_expr_facs = [];
    for i=1:length(factors)
        if typeof(factors[i]) <: Number
            push!(coef_facs, factors[i]);
            push!(coef_inds, -1);
            push!(coef_derivs, [nothing, "", ""]);
            
        elseif typeof(factors[i]) == Expr && factors[i].head === :call
            # These should both be purely coefficient/known expressions. 
            if factors[i].args[1] === :./
                # The second arg should be 1, the third should not contain an unknown or test symbol
                # The denominator expression needs to be processed completely
                (piece, nd, nc, ncm, nci, ncd) = process_known_expr_julia(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                factors[i].args[3] = piece;
                push!(coef_expr_facs, [factors[i], nc, ncm, nci, ncd]);
                
            elseif factors[i].args[1] === :.^
                # The second arg is the thing raised
                (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_julia(factors[i].args[2]);
                need_derivative = need_derivative || nd;
                factors[i].args[2] = piece1;
                # Do the same for the power just in case
                (piece2, nd2, nc2, ncm2, nci2, ncd2) = process_known_expr_julia(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                append!(nc, nc2);
                append!(ncm, ncm2);
                append!(nci, nci2);
                append!(ncd, ncd2);
                factors[i].args[3] = piece2;
                
                push!(coef_expr_facs, [factors[i], nc, ncm, nci, ncd]);
                
            elseif factors[i].args[1] === :sqrt
                factors[i].args[1] = :.^
                # The second arg is the thing sqrted
                (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_julia(factors[i].args[2]);
                need_derivative = need_derivative || nd;
                factors[i].args[2] = piece1;
                # add a 1/2 power argument
                push!(factors[i].args, 1/2);
                
                push!(coef_expr_facs, [factors[i], nc, ncm, nci, ncd]);
            end
            
        else
            (index, BbB, mods) = extract_symbols(factors[i]);
            
            if is_test_func(BbB)
                test_component = index; # the vector index
                if length(mods) > 0
                    # TODO more than one derivative mod
                    need_derivative = true;
                    dmat = Symbol("TRQ"*mods[1][2]);
                    test_part = dmat;
                else
                    # no derivative mods
                    test_part = :(refel.Q');
                end
            elseif is_unknown_var(BbB, var) && lorr == LHS # If rhs, treat as a coefficient
                if !(trial_part === nothing)
                    # Two unknowns multiplied in this term. Nonlinear. abort.
                    printerr("Nonlinear term. Code layer incomplete.");
                    return (-1, -1, -1, -1, -1, -1);
                end
                trial_component = index;
                #offset component for multivar
                trial_var = var;
                if typeof(var) <:Array
                    for vi=1:length(var)
                        if BbB === var[vi].symbol
                            trial_component = trial_component .+ offset_ind[vi];
                            trial_var = var[vi];
                        end
                    end
                end
                if length(mods) > 0
                    # TODO more than one derivative mod
                    need_derivative = true;
                    dmatr = Symbol("RQ"*mods[1][2]);
                    trial_part = dmatr;
                else
                    # no derivative mods
                    trial_part = :(refel.Q);
                end
            else # coefficients
                if length(index) == 1
                    ind = index[1];
                end
                # Check for derivative mods
                if typeof(BbB) == Symbol && !(BbB ===:dt)
                    if length(mods) > 0
                        need_derivative = true;
                        need_derivative_for_coefficient = true;
                        
                        push!(needed_coef_deriv, [BbB, mods[1], mods[1][2]]);
                        
                    else
                        push!(needed_coef_deriv, [BbB, "", ""]);
                    end
                    push!(needed_coef, BbB);
                    push!(needed_coef_name, make_coef_name(BbB, BbB, needed_coef_deriv[end][2], ind););
                    push!(needed_coef_ind, ind);
                    
                    push!(coef_derivs, needed_coef_deriv[end]);
                    
                else
                    push!(coef_derivs, [nothing, "", ""]);
                end
                
                push!(coef_facs, BbB);
                push!(coef_inds, ind);
            end
        end
        
    end # factors loop
    
    # If there's no trial part, need to do this
    if trial_part === nothing
        trial_part = :(refel.Q);
    end
    
    # build coefficient parts
    name_ind = 0;
    if length(coef_facs) > 0
        for j=1:length(coef_facs)
            tmp = coef_facs[j];
            #println("coef_facs: "*string(tmp)*" : "*string(typeof(tmp)));
            if typeof(tmp) == Symbol && !(tmp ===:dt)
                name_ind += 1;
                
                tmps = needed_coef_name[name_ind];
                tmp = Symbol(tmps);
                
            end
            if j>1
                coef_part = :($coef_part .* $tmp);
            else
                coef_part = tmp;
            end
        end
    end
    
    if length(coef_expr_facs) > 0
        for j=1:length(coef_expr_facs)
            tmp = coef_expr_facs[j][1];
            if !(coef_part === nothing)
                coef_part = :($coef_part .* $tmp);
            else
                coef_part = tmp;
            end
            
            append!(needed_coef, coef_expr_facs[j][2]);
            append!(needed_coef_name, coef_expr_facs[j][3]);
            append!(needed_coef_ind, coef_expr_facs[j][4]);
            append!(needed_coef_deriv, coef_expr_facs[j][5]);
        end
    end
    
    # If there's no test part this is probably a denominator expression being processed and should only contain coefficients/knowns
    if test_part === nothing
        #
        #
    else
        term = test_part;
        if !(coef_part === nothing)
            if lorr == LHS
                term = :($test_part * (diagm($weight_part .* $coef_part) * $trial_part));
            else # RHS
                # Stiffness and mass are precomputed for uniform grid meshes
                if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
                    if test_part == :(refel.Q') && weight_part === :wdetj && trial_part == :(refel.Q)
                        term = :(mass * $coef_part);
                        need_derivative = need_derivative_for_coefficient;
                    elseif test_part === :TRQ1 && weight_part === :wdetj && trial_part === :RQ1
                        term = :(stiffness[1] * $coef_part);
                        need_derivative = need_derivative_for_coefficient;
                    elseif test_part === :TRQ2 && weight_part === :wdetj && trial_part === :RQ2
                        term = :(stiffness[2] * $coef_part);
                        need_derivative = need_derivative_for_coefficient;
                    elseif test_part === :TRQ3 && weight_part === :wdetj && trial_part === :RQ3
                        term = :(stiffness[3] * $coef_part);
                        need_derivative = need_derivative_for_coefficient;
                    else
                        term = :($test_part * ($weight_part .* ($trial_part * $coef_part)));
                    end
                else
                    term = :($test_part * ($weight_part .* ($trial_part * $coef_part)));
                end
            end
            
        else
            # Stiffness and mass are precomputed for uniform grid meshes
            if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
                if test_part == :(refel.Q') && weight_part === :wdetj && trial_part == :(refel.Q)
                    term = :(mass);
                    need_derivative = false;
                elseif test_part === :TRQ1 && weight_part === :wdetj && trial_part === :RQ1
                    term = :(stiffness[1]);
                    need_derivative = false;
                elseif test_part === :TRQ2 && weight_part === :wdetj && trial_part === :RQ2
                    term = :(stiffness[2]);
                    need_derivative = false;
                elseif test_part === :TRQ3 && weight_part === :wdetj && trial_part === :RQ3
                    term = :(stiffness[3]);
                    need_derivative = false;
                else
                    term = :($test_part * diagm($weight_part) * $trial_part);
                end
            else
                term = :($test_part * diagm($weight_part) * $trial_part);
            end
        end
        
        if neg
            negex = :(-a);
            negex.args[2] = copy(term);
            term = negex;
        end
    end
    
    return (term, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv, test_component, trial_component);
end

# Special processing for sub expressions of known things.(denominators, sqrt, etc.)
# It should only contain knowns/coeficients, so just replace symbols (f -> coef_0_1)
# And return the needed coefficient info
function process_known_expr_julia(ex)
    need_derivative = false;
    needed_coef = [];
    needed_coef_name = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    
    # Work recursively through the expression
    if typeof(ex) <: Number
        return (ex, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
        
    elseif typeof(ex) == Symbol
        # turn arithmetic ops into dotted versions
        if ex === :+ || ex === :.+
            return (:.+ , false, [], [], [], []);
        elseif ex === :- || ex === :.-
            return (:.- , false, [], [], [], []);
        elseif ex === :* || ex === :.*
            return (:.* , false, [], [], [], []);
        elseif ex === :/ || ex === :./
            return (:./ , false, [], [], [], []);
        elseif ex === :^ || ex === :.^
            return (:.^ , false, [], [], [], []);
        end
        
        (index, v, mods) = extract_symbols(ex);
        if length(index) == 1
            ind = index[1];
        end
        
        if !(v ===:dt)
            if v === :DGNORMAL
                tmp = :(normal[$ind]);
            elseif v === :DGNORMAL1
                tmp = :(normal[$ind]);
            elseif v === :DGNORMAL2
                tmp = :(-normal[$ind]);
            else
                # Check for derivative mods
                if length(mods) > 0 && typeof(v) == Symbol
                    need_derivative = true;
                    
                    push!(needed_coef_deriv, [v, mods[1], mods[1][2]]);
                    
                else
                    push!(needed_coef_deriv, [v, "", ""]);
                end
                
                push!(needed_coef, v);
                push!(needed_coef_name, make_coef_name(v, v, needed_coef_deriv[end][2], ind););
                push!(needed_coef_ind, ind);
                
                tmp = Symbol(needed_coef_name[end]); # The symbol to return
            end
        else
            tmp = ex;
            
        end
        
        return (tmp, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
        
    elseif typeof(ex) == Expr
        newex = copy(ex);
        for i=1:length(ex.args)
            (piece, nd, nc, ncm, nci, ncd) = process_known_expr_julia(ex.args[i]);
            newex.args[i] = piece;
            need_derivative = need_derivative || nd;
            append!(needed_coef, nc);
            append!(needed_coef_name, ncm);
            append!(needed_coef_ind, nci);
            append!(needed_coef_deriv, ncd);
        end
        
        return (newex, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
    end
    
end

###############################################################################################################
# utils
###############################################################################################################

# Separates expressions into terms.
# Assumes the expressions are expanded as they should be coming from the parser.
# 2*a + b*c*5 - w*2 -> [2*a, b*c*5, w*2]
function separate_terms(ex)
    terms = [ex];
    if typeof(ex) == Expr && ex.head === :call
        if ex.args[1] === :+ || (ex.args[1] === :- && length(ex.args) > 2)
            terms = [];
            for i=2:length(ex.args)
                append!(terms, separate_terms(ex.args[i]));
            end
        end
    end
    
    return terms;
end

# Parses symengine terms into julia expressions
function terms_to_expr(symex)
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
        printerr("sorry. still need to implement code layer for tensors.")
    end
    
    # for i=1:length(terms)
    #     convert_sqrt_to_power(terms[i]);
    # end
    
    return terms;
end

# Separates terms into factors.
# Assumes the term is only multiplied symbols or numbers
# 2*a*thing -> [2, a, thing]
function separate_factors(ex)
    factors::Array{Any,1} = [ex];
    if typeof(ex) == Expr && ex.head === :call
        if ex.args[1] === :* || ex.args[1] === :.*
            factors = [];
            for i=2:length(ex.args)
                append!(factors, separate_factors(ex.args[i]));
            end
            
        elseif ex.args[1] === :^ || ex.args[1] === :.^
            factors = [];
            power = ex.args[3];
            if power == 1
                factors = [ex.args[2]]; #a^1 = a
            else
                ex.args[1] = :.^ ;
                factors = [ex]; # the power is handled later
            end
            
        elseif ex.args[1] === :/ || ex.args[1] === :./
            factors = [];
            append!(factors, separate_factors(ex.args[2])); # a/b = a * 1/b
            divex = :(1 ./ a);
            divex.args[3] = ex.args[3];
            push!(factors, divex);
            
        elseif ex.args[1] === :sqrt
            ex.args[1] = :.^ ;
            push!(ex.args, 0.5);
            factors = [ex];
            
        elseif ex.args[1] === :- && length(ex.args) == 2
            # strip off negetive and place on first factor
            subex = ex.args[2];
            factors = separate_factors(subex);
            negex = :(-a);
            negex.args[2] = factors[1];
            factors[1] = negex;
        end
    end
    
    return factors;
end

# Checks if the coefficient has constant value.
# If so, also returns the value.
function is_constant_coef(c)
    isit = false;
    val = 0;
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            isit = (typeof(coefficients[i].value[1]) <: Number);
            if isit
                val = coefficients[i].value[1];
            else
                val = coefficients[i].value[1].name;
            end
        end
    end
    
    return (isit, val);
end

# Checks the type of coefficient: constant, genfunction, or variable
# Returns: (type, val)
# constant: type=1, val=number
# genfunction: type=2, val= index in genfunctions array
# variable: type=3, val=index in variables array
function get_coef_val(c, comp)
    type = 0;
    val = 0;
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            isit = (typeof(coefficients[i].value[comp]) <: Number);
            if isit
                type = 1;
                val = coefficients[i].value[comp];
            else
                type = 2;
                name = coefficients[i].value[comp].name;
                for j=1:length(genfunctions)
                    if name == genfunctions[j].name
                        val = j;
                    end
                end
            end
        end
    end
    if type == 0
        for i=1:length(variables)
            if c === variables[i].symbol
                type = 3;
                val = variables[i].index;
            end
        end
    end
    
    return (type, val);
end

function get_coef_index(c)
    ind = -1;
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            ind = coefficients[i].index;
        end
    end
    
    return ind;
end

function make_coef_name(c, tmptag, deriv, cind)
    ind = get_coef_index(c);
    if ind >= 0
        tag = string(ind);
    else
        tag = string(tmptag);
    end
    
    tag = deriv * tag;
    str = "coef_"*tag*"_"*string(cind);
    return str;
end

function is_test_func(v)
    for t in test_functions
        if t.symbol === v
            return true;
        end
    end
    return false;
end

function is_unknown_var(v, vars)
    if typeof(vars) == Variable
        return v===vars.symbol;
    end
    for t in vars
        if t.symbol === v
            return true;
        end
    end
    return false;
end

# Extract meaning from the symbolic object name
# The format of the input symbol should look like this
#   MOD1_MOD2_..._var_n
# There could be any number of mods, _var_n is the symvar symbol (n is the vector index, or 1 for scalar)
# Returns ([n], var, [MOD1,MOD2,...]) all as strings
function extract_symbols(ex)
    str = string(ex);
    #println("extracting from: "*str);
    index = [];
    var = nothing;
    mods = [];
    l = lastindex(str);
    e = l; # end of variable name
    b = l; # beginning of variable name
    
    # dt is a special symbol that will be assigned a number value in the generated function.
    if str == "dt"
        return([0], ex, []);
    end
    
    if occursin("TEST", str)
        # index = [];
        var = :TEST;
        for i=l:-1:0
            if e==l
                if str[i] == '_'
                    e = i-1;
                else
                    # These are the indices at the end. Parse one digit at a time.
                    index = [parse(Int, str[i]); index] # The indices on the variable
                end
            end
        end
        e = e-5;
        b = e;
    else
        for i=l:-1:0
            if e==l
                if str[i] == '_'
                    e = i-1;
                else
                    # These are the indices at the end. Parse one digit at a time.
                    try
                        index = [parse(Int, str[i]); index] # The indices on the variable
                    catch
                        return ([0],ex,[]);
                    end
                    
                end
            elseif b==l && i>0
                if str[i] == '_'
                    b = i+1;
                end
                
            else
                # At this point we know b and e
                if var === nothing
                    var = Symbol(SubString(str, b, e));
                    b = b-2;
                    e = b;
                end
            end
        end
    end
    
    # Change index from [a, b] to [a*d + b]
    newindex = index[end];
    for i=1:(length(index)-1)
        newindex = newindex + (index[i]-1)*config.dimension^(length(index)-i);
    end
    index = [newindex];
    
    # extract the modifiers like D1_ -> :Dx
    if b>1
        e = b-1;
        for i=e:-1:1
            if str[i] == '_'
                push!(mods, SubString(str, b, e));
                
                e = i-1;
                
            elseif i == 1
                push!(mods, SubString(str, 1, e));
            end
            b = i;
        end
    end
    
    #println("got: "*string(mods)*" "*string(var)*" "*string(index));
    
    return (index, var, mods);
end

# Separate the derivative mods(Dn_) from the other mods and return an array 
# of derivative indices and the remaining mods.
function get_derivative_mods(mods)
    derivs = [];
    others = [];
    
    for i=1:length(mods)
        if length(mods[i]) == 2 && mods[i][1] == 'D'
            try
                index = parse(Int, mods[i][2])
                push!(derivs, index);
            catch
                printerr("Unexpected modifier: "*mods[i]*" see get_derivative_mods() in generate_code_layer")
                push!(others, mods[i]);
            end
        else
            push!(others, mods[i]);
        end
    end
    
    return (derivs, others);
end