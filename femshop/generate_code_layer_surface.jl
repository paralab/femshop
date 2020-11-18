#=
This generates for the surface integrals
=#

function generate_code_layer_surface(ex, var, lorr)
    if use_cachesim
        if language == 0 || language == JULIA
            printerr("surface integrals not ready for cachesim")
            return nothing;
        elseif language == CPP
            printerr("surface integrals not ready for cachesim")
            return "";
        elseif language == MATLAB
            printerr("surface integrals not ready for cachesim")
            return "";
        end
    else
        if language == 0 || language == JULIA
            return generate_code_layer_julia_surface(ex, var, lorr);
        elseif language == CPP
            #return generate_code_layer_dendro_surface(ex, var, lorr);
            printerr("surface integrals only ready for Julia");
            return "";
        elseif language == MATLAB
            #return generate_code_layer_homg_surface(ex, var, lorr);
            printerr("surface integrals only ready for Julia");
            return "";
        end
    end
end

###############################################################################################################
# julia
###############################################################################################################

# Julia version returns an expression for the generated function for linear or bilinear term
function generate_code_layer_julia_surface(symex, var, lorr)
    # This is the basic info passed in "args"
    # args = (var, val_s1, node_s1, normal_s1, xe1, val_s2, node_s2, normal_s2, xe2, refel, face_refel, RHS/LHS, t, dt);
    code = Expr(:block);
    push!(code.args, :(var = args[1]));         # list of unknown variables for this expression
    push!(code.args, :(val_s1 = args[2]));      # ??
    push!(code.args, :(node_s1 = args[3]));     # global indices of the nodes on side 1
    push!(code.args, :(normal_s1 = args[4]));   # normal vector on side 1
    push!(code.args, :(xf1 = args[5]));         # global coords of element's nodes on side 1
    push!(code.args, :(val_s2 = args[6]));      # ??
    push!(code.args, :(node_s2 = args[7]));     # global indices of the nodes on side 2
    push!(code.args, :(normal_s2 = args[8]));   # normal vector on side 2
    push!(code.args, :(xf2 = args[9]));         # global coords of element's nodes on side 1
    push!(code.args, :(refel = args[10]));      # reference element for volume
    push!(code.args, :(face_refel = args[11])); # reference element for face
    push!(code.args, :(borl = args[12]));       # bilinear or linear? lhs or rhs?
    push!(code.args, :(time = args[13]));       # time for time dependent coefficients
    push!(code.args, :(dt = args[14]));         # dt for time dependent problems
    
    # Build geometric factors for both volume and surface
    push!(code.args, :((detJ, J) = geometric_factors(refel, x)));
    push!(code.args, :((face_detJ, face_J) = geometric_factors(face_refel, xf1)));
    push!(code.args, :(wgdetj = refel.wg .* detJ));
    push!(code.args, :(face_wgdetj = face_refel.wg .* face_detJ));
    
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
                (codeterm, der, coe, coeind, coederiv, testi, trialj) = process_term_julia(terms[vi][i], var, lorr, offset_ind);
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
            (codeterm, der, coe, coeind, coederiv, testi, trialj) = process_term_julia(terms[i], var, lorr);
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
            #println("coef: "*string(needed_coef[i])*" , ind: "*string(needed_coef_ind[i])*" , deriv: "*string(needed_coef_deriv[i]));
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
                    push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                elseif ctype == 3
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
                
                #push!(code.args, :(println( $tmpc )));
                push!(code.args, Expr(:(=), tmpc, tmpb));
                #push!(code.args, :(println( $tmpc )));
                #push!(code.args, :(println("-----------------------------------------------------------------------------")));
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
        #     # add terms into full matrix according to testind/trialind
        #     for i=1:length(terms[vi])
        #         ti = test_ind[vi][i][1]-1 + offset_ind[vi];
        #         tj = trial_ind[vi][i][1]-1;
                
        #         sti = :($ti*refel.Np + 1);
        #         eni = :(($ti + 1)*refel.Np);
        #         stj = :($tj*refel.Np + 1);
        #         enj = :(($tj + 1)*refel.Np);
        #         #submat = Symbol("element_matrix_"*string(var[vi].symbol));
        #         #subvec = Symbol("element_vector_"*string(var[vi].symbol));
                
        #         if lorr == LHS
        #             push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), terms[vi][i]));
        #         else
        #             push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), terms[vi][i]));
        #         end
        #     end
        # end
        # if lorr == LHS
        #     result = :element_matrix;
        # else
        #     result = :element_vector;
        # end
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
                        
                        if submatrices[ti, tj] === nothing
                            submatrices[ti, tj] = terms[i];
                        else
                            addexpr = :(a+b);
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
                            addexpr = :(a+b);
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
                # # add terms into full matrix according to testind/trialind
                # tmp = :(a+b);
                # tmp.args = [:+];
                # for i=1:length(terms)
                #     ti = test_ind[i][1]-1;
                #     tj = trial_ind[i][1]-1;
                    
                #     st = :($ti*refel.Np + 1);
                #     en = :(($ti + 1)*refel.Np);
                    
                #     if lorr == LHS
                #         push!(code.args, Expr(:(+=), :(element_matrix[$st:$en, $st:$en]), terms[i]));
                #     else
                #         push!(code.args, Expr(:(+=), :(element_vector[$st:$en]), terms[i]));
                #     end
                    
                # end
                # if lorr == LHS
                #     result = :element_matrix;
                # else
                #     result = :element_vector;
                # end
                
            end
        elseif length(terms) == 1# one term (one variable)
            result = terms[1];
            
        else # there were no terms. Just return arrays of zeros
            if lorr == LHS
                return :(return zeros(face_refel.Np, face_refel.Np));
            else
                return :(return zeros(face_refel.Np));
            end
        end
    end
    
    
    # At this point everything is packed into terms[1]
    push!(code.args, Expr(:return, result));
    return code;
end
