###############################################################################################################
# generate for cachesim
###############################################################################################################
#=
Want

l 2342
s 512 8
l 512 8

to be turned into

cs.load(2342)  # Loads one byte from address 2342
cs.store(512, length=8)  # Stores 8 bytes to addresses 512-519
cs.load(512, length=8)  # Loads from address 512 until (exclusive) 520 (eight bytes)
...
cs.force_write_back()
cs.print_stats()

The arrays are indexed as such:
1 A
2 b
3 u
4 Ak
5 bk
6 Q
7 RQ1
8 RQ2
9 RQ3
10 TRQ1
11 TRQ2
12 TRQ3
13 RD1
14 RD2
15 RD3
17 wdetJ
18+ any other needed arrays
=#

function generate_code_layer_cachesim(symex, var, lorr)
    # This is the basic info passed in "args"
    code = Expr(:block);
    push!(code.args, :(var = args[1]));
    push!(code.args, :(x = args[2]));    # global coords of element's nodes
    push!(code.args, :(gbl = args[3]));  # global indices of the nodes
    push!(code.args, :(refel = args[4]));# 
    push!(code.args, :(borl = args[5])); # bilinear or linear? code or rhs?
    push!(code.args, :(time = args[6])); # time for time dependent coefficients
    if prob.time_dependent
        push!(code.args, :(dt = args[7])); # dt for time dependent problems
    end
    #push!(code.args, :((detJ, J) = geometric_factors(refel, x)));
    
    #push!(code.args, :(wgdetj = refel.wg .* detJ));
    push!(code.args, :(Femshop.cachesim_store_range(13))); # cachesim
    push!(code.args, :(Femshop.cachesim_store_range(14))); # cachesim
    push!(code.args, :(Femshop.cachesim_store_range(15))); # cachesim
    push!(code.args, :(Femshop.cachesim_store_range(16))); # cachesim
    
    need_derivative = false;
    needed_coef = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    test_ind = [];
    trial_ind = [];
    
    # For multi variables
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
                (codeterm, der, coe, coeind, coederiv, testi, trialj) = process_term_cachesim(terms[vi][i], var, lorr, offset_ind);
                if coeind == -1
                    # processing failed due to nonlinear term
                    printerr("term processing failed for: "*string(terms[vi][i])*" , possible nonlinear term?");
                    return nothing;
                end
                need_derivative = need_derivative || der;
                append!(needed_coef, coe);
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
            (codeterm, der, coe, coeind, coederiv, testi, trialj) = process_term_cachesim(terms[i], var, lorr);
            if coeind == -1
                # processing failed due to nonlinear term
                printerr("term processing failed for: "*string(terms[i])*" , possible nonlinear term?");
                return nothing;
            end
            need_derivative = need_derivative || der;
            append!(needed_coef, coe);
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
        push!(code.args, :(Femshop.cachesim_load_range(13))); # cachesim
        push!(code.args, :(Femshop.cachesim_load_range(14))); # cachesim
        push!(code.args, :(Femshop.cachesim_load_range(15))); # cachesim
        push!(code.args, :(Femshop.cachesim_load_range(16))); # cachesim
        if config.dimension == 1
            #push!(code.args, :((RQ1,RD1) = build_deriv_matrix(refel, J)));
            #push!(code.args, :(TRQ1 = RQ1'));
            push!(code.args, :(Femshop.cachesim_load_range(7))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(10))); # cachesim
        elseif config.dimension == 2
            # push!(code.args, :((RQ1,RQ2,RD1,RD2) = build_deriv_matrix(refel, J)));
            # push!(code.args, :((TRQ1,TRQ2) = (RQ1',RQ2')));
            push!(code.args, :(Femshop.cachesim_load_range(7))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(8))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(10))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(11))); # cachesim
        elseif config.dimension == 3
            # push!(code.args, :((RQ1,RQ2,RQ3,RD1,RD2,RD3) = build_deriv_matrix(refel, J)));
            # push!(code.args, :((TRQ1,TRQ2,TRQ3) = (RQ1',RQ2',RQ3')));
            push!(code.args, :(Femshop.cachesim_load_range(7))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(8))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(9))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(10))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(11))); # cachesim
            push!(code.args, :(Femshop.cachesim_load_range(12))); # cachesim
        end
    end
    
    # If coefficients need to be computed, do so
    # # First remove duplicates
    unique_coef = [];
    unique_coef_ind = [];
    unique_coef_deriv = [];
    for i=1:length(needed_coef)
        already = false;
        for j=1:length(unique_coef)
            if unique_coef[j] === needed_coef[i] && unique_coef_ind[j] == needed_coef_ind[i] && unique_coef_deriv[j] == needed_coef_deriv[i]
                already = true;
            end
        end
        if !already
            push!(unique_coef, needed_coef[i]);
            push!(unique_coef_ind, needed_coef_ind[i]);
            push!(unique_coef_deriv, needed_coef_deriv[i]);
        end
    end
    needed_coef = unique_coef;
    needed_coef_ind = unique_coef_ind;
    needed_coef_deriv = unique_coef_deriv;
    
    # For constant coefficients, this generates something like:
    ######################################
    # coef_n_i = a.value[i];
    ######################################
    
    # For variable coefficients, this generates something like:
    ######################################
    # coef_n_i = zeros(refel.Np);
    # for coefi = 1:refel.Np
    #     coef_n_i[coefi] = a.value[i].func(x[coefi,1], x[coefi,2],x[coefi,3],time);
    # end
    ######################################
    if length(needed_coef) > 0
        cloop = :(for coefi=1:refel.Np end);
        cloopin = Expr(:block);
        cargs = [:(x[coefi]); 0; 0; :time];
        if config.dimension == 2
            cargs = [:(x[coefi,1]); :(x[coefi,2]); 0; :time];
        elseif config.dimension == 3
            cargs = [:(x[coefi,1]); :(x[coefi,2]); :(x[coefi,3]); :time];
        end
        
        for i=1:length(needed_coef)
            if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                cind = get_coef_index(needed_coef[i]);
                if cind >= 0
                    tag = string(cind);
                else
                    tag = string(needed_coef[i]);
                end
                #derivatives of coefficients
                tag = needed_coef_deriv[i][2] * tag;
                tmps = "coef_"*tag*"_"*string(needed_coef_ind[i]);
                tmpc = Symbol(tmps);
                
                (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
                if ctype == 1
                    # constant coefficient -> coef_n = cval
                    tmpn = cval;
                    push!(code.args, Expr(:(=), tmpc, tmpn));
                elseif ctype == 2
                    # genfunction coefficients -> coef_n_i = coef.value[i].func(cargs)
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc;
                    #tmpn = :(a.value[1]);
                    #tmpn.args[2] = needed_coef_ind[i];
                    #tmpn.args[1].args[1] = :(Femshop.genfunctions[$cval]); # Femshop.genfunctions[cval]
                    tmpn = :(Femshop.genfunctions[$cval]); # Femshop.genfunctions[cval]
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    #push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    cssymb = Symbol(tmps*"csa");
                    push!(code.args, :($cssymb = Femshop.add_cachesim_tmp_array())); # cachesim
                    push!(code.args, :(Femshop.cachesim_store_range($cssymb))); # cachesim
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                elseif ctype == 3
                    # variable values -> coef_n = variable.values
                    csarrayind = 17+cval;
                    if variables[cval].type == SCALAR
                        tmpb = :(copy(Femshop.variables[$cval].values[gbl])); 
                        push!(code.args, :(Femshop.cachesim_load_range($csarrayind, gbl))); # cachesim
                    else
                        compo = needed_coef_ind[i];
                        tmpb = :(copy(Femshop.variables[$cval].values[gbl, $compo]));
                        push!(code.args, :(Femshop.cachesim_load_range($csarrayind, gbl, [$compo]))); # cachesim
                    end
                    
                    #push!(code.args, Expr(:(=), tmpc, tmpb));
                end
                
            end# number?
        end# coef loop
        
        # Write the loop that computes coefficient values
        if length(cloopin.args) > 0
            cloop.args[2] = cloopin;
            #push!(code.args, cloop); # add loop to code
        end
        
        # Apply derivatives after the initializing loop
        # Will look like: coef_i_j = Rnmatrix * Qnmatrix * coef_i_j
        for i=1:length(needed_coef_deriv)
            if length(needed_coef_deriv[i][2]) > 0 && !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                cind = get_coef_index(needed_coef[i]);
                
                if cind >= 0
                    tag = string(cind);
                else
                    tag = string(needed_coef[i]);
                end
                #derivatives of coefficients
                tag = needed_coef_deriv[i][2] * tag;
                tmps = "coef_"*tag*"_"*string(needed_coef_ind[i]);
                tmpc = Symbol(tmps);
                
                dmat = Symbol("RD"*needed_coef_deriv[i][3]);
                tmpb= :(length($tmpc) == 1 ? 0 : $dmat * $tmpc);
                #dmatr = Symbol("R"*needed_coef_deriv[i][3]*"matrix");
                #dmatq = Symbol("D"*needed_coef_deriv[i][3]*"matrix");
                #tmpb= :(length($tmpc) == 1 ? 0 : $dmatr * $dmatq * $tmpc);
                
                #push!(code.args, Expr(:(=), tmpc, tmpb));
            else
                # derivatives of constant coefficients should be zero
                # TODO
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
            #push!(code.args, Expr(:(=), :element_vector, :(zeros(refel.Np*$dofsper)))); # allocate vector
            push!(code.args, :(Femshop.cachesim_store_range(5))); # cachesim
        else
            #push!(code.args, Expr(:(=), :element_matrix, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper)))); # allocate matrix
            push!(code.args, :(Femshop.cachesim_store_range(4))); # cachesim
        end
    end
    
    result = nothing; # Will hold the returned expression
    if typeof(var) <: Array # multivar
        for vi=1:length(var)
            # add terms into full matrix according to testind/trialind
            for i=1:length(terms[vi])
                ti = test_ind[vi][i][1]-1 + offset_ind[vi];
                tj = trial_ind[vi][i][1]-1;
                
                sti = :($ti*refel.Np + 1);
                eni = :(($ti + 1)*refel.Np);
                stj = :($tj*refel.Np + 1);
                enj = :(($tj + 1)*refel.Np);
                
                if lorr == LHS
                    #push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), terms[vi][i]));
                    #push!(code.args, :(Femshop.cachesim_load_range(4, $sti:$eni, $stj:$enj))); # cachesim
                    #push!(code.args, :(Femshop.cachesim_store_range(4, $sti:$eni, $stj:$enj))); # cachesim
                else
                    #push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), terms[vi][i]));
                    #push!(code.args, :(Femshop.cachesim_load_range(5, $sti:$eni))); # cachesim
                    #push!(code.args, :(Femshop.cachesim_store_range(5, $sti:$eni))); # cachesim
                end
            end
        end
        if lorr == LHS
            result = :element_matrix;
        else
            result = :element_vector;
        end
    else
        if length(terms) > 1
            if var.type == SCALAR # Only one component
                tmp = :(a+b);
                tmp.args = [:+];
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
                        
                        if submatrices[i][j] === nothing
                            submatrices[i][j] = terms[i];
                        else
                            addexpr = :(a+b);
                            addexpr.args[2] = submatrices[i][j];
                            addexpr.args[3] = terms[i];
                        end
                    end
                    
                    for cj=1:comps
                        for ci=1:comps
                            if !(submatrices[ci][cj] === nothing)
                                ti = ci-1;
                                tj = cj-1;
                                sti = :($ti*refel.Np + 1);
                                eni = :(($ti + 1)*refel.Np);
                                stj = :($tj*refel.Np + 1);
                                enj = :(($tj + 1)*refel.Np);
                                
                                #push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), submatrices[ci][cj]));
                                push!(code.args, :(Femshop.cachesim_load_range(4, $sti:$eni, $stj:$enj))); # cachesim
                                push!(code.args, :(Femshop.cachesim_store_range(4, $sti:$eni, $stj:$enj))); # cachesim
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
                        
                        if submatrices[i] === nothing
                            submatrices[i] = terms[i];
                        else
                            addexpr = :(a+b);
                            addexpr.args[2] = submatrices[i];
                            addexpr.args[3] = terms[i];
                        end
                    end
                    
                    for ci=1:comps
                        if !(submatrices[ci] === nothing)
                            ti = ci-1;
                            sti = :($ti*refel.Np + 1);
                            eni = :(($ti + 1)*refel.Np);
                            
                            #push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), submatrices[ci]));
                            push!(code.args, :(Femshop.cachesim_load_range(5, $sti:$eni))); # cachesim
                            push!(code.args, :(Femshop.cachesim_store_range(5, $sti:$eni))); # cachesim
                        end
                    end
                    
                    result = :element_vector;
                    
                end
            end
        else# one term (one variable)
            result = terms[1];
        end
    end
    
    push!(code.args, :(Femshop.remove_all_tmp_arrays()));
    
    # At this point everything is packed into terms[1]
    result = 0;
    push!(code.args, Expr(:return, result));
    return code;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term_cachesim(sterm, var, lorr, offset_ind=0)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    weight_part = :wgdetj;
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
    for i=1:length(factors)
        if typeof(factors[i]) <: Number
            push!(coef_facs, factors[i]);
            push!(coef_inds, -1);
            
        elseif typeof(factors[i]) == Expr && factors[i].head === :call
            # These should both be purely coefficient/known expressions. 
            if factors[i].args[1] === :./
                # The second arg should be 1, the third should not contain an unknown or test symbol
                # The denominator expression needs to be processed completely
                (piece, nd, nc, nci, ncd) = process_known_expr_cachesim(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[3] = piece;
                push!(coef_facs, factors[i]);
                push!(coef_inds, 0);
                
            elseif factors[i].args[1] === :.^
                # The second arg is the thing raised
                (piece1, nd, nc, nci, ncd) = process_known_expr_cachesim(factors[i].args[2]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[2] = piece1;
                # Do the same for the power just in case
                (piece2, nd, nc, nci, ncd) = process_known_expr_cachesim(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[3] = piece2;
                
                push!(coef_facs, factors[i]);
                push!(coef_inds, 0);
            end
            
        else
            (index, v, mods) = extract_symbols(factors[i]);
            
            if is_test_func(v)
                test_component = index; # the vector index
                if length(mods) > 0
                    # TODO more than one derivative mod
                    need_derivative = true;
                    dmat = Symbol("TRQ"*mods[1][2]);
                    test_part = dmat;
                    # dmatr = Symbol("R"*mods[1][2]*"matrix");
                    # dmatq = Symbol("Q"*mods[1][2]*"matrix");
                    # test_part = :(transpose($dmatr * $dmatq));
                else
                    # no derivative mods
                    test_part = :(refel.Q');
                end
            elseif is_unknown_var(v, var) && lorr == LHS # If rhs, treat as a coefficient
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
                        if v === var[vi].symbol
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
                    # dmatr = Symbol("R"*mods[1][2]*"matrix");
                    # dmatq = Symbol("Q"*mods[1][2]*"matrix");
                    # trial_part = :($dmatr * $dmatq);
                else
                    # no derivative mods
                    trial_part = :(refel.Q);
                    # if lorr == RHS # If rhs, change var into var.values and treat as a coefficient
                    #     if length(trial_var.symvar.vals) > 1
                    #         # need a component index in there
                    #         tmpv = :(a.values[gbl,$trial_component]);
                    #     else
                    #         tmpv = :(a.values[gbl]);
                    #     end
                    #     tmpv.args[1].args[1] = trial_var.symbol;
                    #     push!(coef_facs, tmpv);
                    #     push!(coef_inds, trial_component);
                    # else
                    #     trial_part = :(refel.Q);
                    # end
                end
            else # coefficients
                if length(index) == 1
                    ind = index[1];
                end
                # Check for derivative mods
                if length(mods) > 0 && typeof(v) == Symbol && !(v ===:dt)
                    need_derivative = true;
                    
                    push!(needed_coef_deriv, [v, mods[1], mods[1][2]]);
                    
                elseif !(v ===:dt)
                    push!(needed_coef_deriv, [v, "", ""]);
                end
                
                push!(coef_facs, v);
                push!(coef_inds, ind);
            end
        end
        
    end # factors loop
    
    # If there's no trial part, need to do this
    if trial_part === nothing
        trial_part = :(refel.Q);
    end
    
    # build coefficient parts
    if length(coef_facs) > 0
        for j=1:length(coef_facs)
            tmp = coef_facs[j];
            #println("coef_facs: "*string(tmp)*" : "*string(typeof(tmp)));
            if typeof(tmp) == Symbol && !(tmp ===:dt)
                push!(needed_coef, tmp);
                push!(needed_coef_ind, coef_inds[j]);
                ind = get_coef_index(coef_facs[j]);
                if ind >= 0
                    tag = string(ind);
                else
                    tag = string(tmp);
                end
                #derivatives of coefficients
                tag = needed_coef_deriv[length(needed_coef)][2] * tag;
                tmps = "coef_"*tag*"_"*string(coef_inds[j]);
                tmp = Symbol(tmps);
                
            end
            if j>1
                coef_part = :($coef_part .* $tmp);
            else
                coef_part = tmp;
            end
        end
    end
    
    # If there's no test part this is probably a denominator expression being processed and should only contain coefficients/knowns
    if test_part === nothing
        
        
    else
        term = test_part;
        if !(coef_part === nothing)
            if lorr == LHS
                term = :($test_part * (diagm($weight_part .* $coef_part) * $trial_part));
            else # RHS
                term = :($test_part * ($weight_part .* ($trial_part * $coef_part)));
                #term = :($test_part * (diagm($weight_part) * ($trial_part * $coef_part)));
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
    
    return (term, need_derivative, needed_coef, needed_coef_ind, needed_coef_deriv, test_component, trial_component);
end

# Special processing for sub expressions of known things.(denominators, sqrt, etc.)
# It should only contain knowns/coeficients, so just replace symbols (f -> coef_0_1)
# And return the needed coefficient info
function process_known_expr_cachesim(ex)
    need_derivative = false;
    needed_coef = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    
    # Work recursively through the expression
    if typeof(ex) <: Number
        return (ex, need_derivative, needed_coef, needed_coef_ind, needed_coef_deriv);
        
    elseif typeof(ex) == Symbol
        # turn arithmetic ops into dotted versions
        if ex === :+ || ex === :.+
            return (:.+ , false, [], [], []);
        elseif ex === :- || ex === :.-
            return (:.- , false, [], [], []);
        elseif ex === :* || ex === :.*
            return (:.* , false, [], [], []);
        elseif ex === :/ || ex === :./
            return (:./ , false, [], [], []);
        elseif ex === :^ || ex === :.^
            return (:.^ , false, [], [], []);
        end
        
        (index, v, mods) = extract_symbols(ex);
        if length(index) == 1
            ind = index[1];
        end
        
        if !(v ===:dt)
            # Check for derivative mods
            if length(mods) > 0 && typeof(v) == Symbol
                need_derivative = true;
                
                push!(needed_coef_deriv, [v, mods[1], mods[1][2]]);
                
            else
                push!(needed_coef_deriv, [v, "", ""]);
            end
            
            push!(needed_coef, v);
            push!(needed_coef_ind, ind);
            
            cind = get_coef_index(v);
            if cind >= 0
                tag = string(cind);
            else
                tag = string(v);
            end
            #derivatives of coefficients
            tag = needed_coef_deriv[length(needed_coef)][2] * tag;
            tmps = "coef_"*tag*"_"*string(ind);
            tmp = Symbol(tmps); # The symbol to return
            
        else
            tmp = ex;
            
        end
        
        return (tmp, need_derivative, needed_coef, needed_coef_ind, needed_coef_deriv);
        
    elseif typeof(ex) == Expr
        newex = copy(ex);
        for i=1:length(ex.args)
            (piece, nd, nc, nci, ncd) = process_known_expr_cachesim(ex.args[i]);
            newex.args[i] = piece;
            need_derivative = need_derivative || nd;
            append!(needed_coef, nc);
            append!(needed_coef_ind, nci);
            append!(needed_coef_deriv, ncd);
        end
        
        return (newex, need_derivative, needed_coef, needed_coef_ind, needed_coef_deriv);
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

# Parses symengine terms into cachesim expressions
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
                    try
                        index = [parse(Int, str[i]); index] # The indices on the variable
                    catch
                        return ([0],ex,[]);
                    end
                    
                end
            elseif b==l
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

