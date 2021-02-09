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
    ################## TEMP
    # if lorr==LHS
    #     return :(A=zeros(args[3].Np, args[3].Np); return (A, A, A, A));
    # end
    # Now we can safely assume RHS
    #########################
    
    # This is the basic info passed in "args"
    # args = (var, fid, frefelind, facenodes, face2glb, normal, faceBID, fdetJ, fJ, face_wgdetj, RHS, t, dt);
    code = Expr(:block);
    
    push!(code.args, :(var =        args[1]));  # list of unknown variables for this expression
    push!(code.args, :(refel =      args[2]));  # reference element for volume
    push!(code.args, :(loc2glb =    args[3]));  # loc2glb for full elements
    push!(code.args, :(fid =        args[4]));  # face index
    push!(code.args, :(frefelind =     args[5]));  # reference element indices for faces
    push!(code.args, :(fnodes =     args[6]));  # global coords of face nodes 
    #push!(code.args, :(flocal =     args[7]));  # local coords of face nodes within element refel
    push!(code.args, :(face2glb =   args[7]));  # Global indices of inside face nodes
    push!(code.args, :(normal =     args[8]));  # Global indices of outside face nodes
    push!(code.args, :(faceBID =    args[9])); # Coordinates of face nodes
    push!(code.args, :(face_detJ =  args[10])); # geometric factor for face
    push!(code.args, :(face_J =     args[11])); # geometric factor for face
    push!(code.args, :(vol_J1 =     args[12])); # geometric factor for el1
    push!(code.args, :(vol_J2 =     args[13])); # geometric factor for el2
    push!(code.args, :(face_wgdetj =args[14])); # quadrature weights*detJ
    push!(code.args, :(borl =       args[15])); # bilinear or linear? lhs or rhs?
    push!(code.args, :(time =       args[16])); # time for time dependent coefficients
    push!(code.args, :(dt =         args[17])); # dt for time dependent problems
    
    # push!(code.args, :(Q1 = zeros(size(fnodes,2))));
    # push!(code.args, :(Q2 = zeros(size(fnodes,2))));
    # Qloop = :(for i=1:size(refel.surf_Q[frefelind[1]], 2) end);
    # Qloopin = Expr(:block);
    # push!(Qloopin.args, :(Q1 = Q1 + refel.surf_Q[frefelind[1]][:,i]));
    # push!(Qloopin.args, :(Q2 = Q2 + refel.surf_Q[frefelind[2][:,i]));
    # Qloop.args[2] = Qloopin;
    # push!(code.args, Qloop);
    
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
    
    # 4 sets
    terms = [copy(terms),copy(terms),copy(terms),copy(terms)];
    
    # Process the terms turning them into the code layer
    if multivar
        # for vi=1:varcount
        #     subtest_ind = [];
        #     subtrial_ind = [];
        #     for i=1:length(terms[vi])
        #         (codeterm, der, coe, coeind, coederiv, testi, trialj) = process_surface_term_julia(terms[vi][i], var, lorr, offset_ind);
        #         if coeind == -1
        #             # processing failed due to nonlinear term
        #             printerr("term processing failed for: "*string(terms[vi][i])*" , possible nonlinear term?");
        #             return nothing;
        #         end
        #         need_derivative = need_derivative || der;
        #         append!(needed_coef, coe);
        #         append!(needed_coef_ind, coeind);
        #         append!(needed_coef_deriv, coederiv);
        #         # change indices into one number
                
        #         push!(subtest_ind, testi);
        #         push!(subtrial_ind, trialj);
        #         terms[vi][i] = codeterm;
        #     end
        #     push!(test_ind, subtest_ind);
        #     push!(trial_ind, subtrial_ind);
        # end
    else
        for i=1:length(terms[1])
            (codeterms, der, coe, coeind, coederiv, testi, trialj) = process_surface_term_julia(terms[1][i], var, lorr);
            if coeind == -1
                # processing failed due to nonlinear term
                printerr("term processing failed for: "*string(terms[1][i])*" , possible nonlinear term?");
                return nothing;
            end
            need_derivative = need_derivative || der;
            append!(needed_coef, coe);
            append!(needed_coef_ind, coeind);
            append!(needed_coef_deriv, coederiv);
            
            push!(test_ind, testi);
            push!(trial_ind, trialj);
            terms[1][i] = codeterms[1];
            terms[2][i] = codeterms[2];
            terms[3][i] = codeterms[3];
            terms[4][i] = codeterms[4];
        end
    end
    
    # If derivatives are needed, prepare the appropriate matrices
    if need_derivative
        if config.dimension == 1
            push!(code.args, :(RQ1_1 = refel.surf_Qr[frefelind[1]] .* vol_J1.rx[1]));
            push!(code.args, :(RQ2_1 = refel.surf_Qr[frefelind[2]] .* vol_J2.rx[1]));
            push!(code.args, :(TRQ1_1 = RQ1_1'));
            push!(code.args, :(TRQ2_1 = RQ2_1'));
        elseif config.dimension == 2
            push!(code.args, :((RQ1_1,RQ1_2,RD1_1,RD1_2) = build_face_deriv_matrix(refel, frefelind[1], vol_J1)));
            push!(code.args, :((RQ2_1,RQ2_2,RD2_1,RD2_2) = build_face_deriv_matrix(refel, frefelind[2], vol_J2)));
            push!(code.args, :((TRQ1_1,TRQ1_2) = (RQ1_1',RQ1_2')));
            push!(code.args, :((TRQ2_1,TRQ2_2) = (RQ2_1',RQ2_2')));
        elseif config.dimension == 3
            push!(code.args, :((RQ1_1,RQ1_2,RQ1_3,RD1_1,RD1_2,RD1_3) = build_face_deriv_matrix(refel, frefel[1], vol_J1)));
            push!(code.args, :((RQ2_1,RQ2_2,RQ2_3,RD2_1,RD2_2,RD2_3) = build_face_deriv_matrix(refel, frefel[2], vol_J2)));
            push!(code.args, :((TRQ1_1,TRQ1_2,TRQ1_3) = (RQ1_1',RQ1_2',RQ1_3')));
            push!(code.args, :((TRQ2_1,TRQ2_2,TRQ2_3) = (RQ2_1',RQ2_2',RQ2_3')));
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
            if unique_coef[j] === needed_coef[i] && unique_coef_ind[j] == needed_coef_ind[i]
                if unique_coef_deriv[j] == needed_coef_deriv[i] || (needed_coef_deriv[i][3] == "" && unique_coef_deriv[j][3] == "")
                    already = true;
                end
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
    # coef_n_i = zeros(face_refel.Np);
    # for coefi = 1:face_refel.Np
    #     coef_n_i[coefi] = a.value[i].func(x[coefi,1], x[coefi,2],x[coefi,3],time);
    # end
    ######################################
    if length(needed_coef) > 0
        cloop = :(for coefi=1:refel.Nfp[frefelind[1]] end);
        cloopin = Expr(:block);
        cargs = [:(fnodes[coefi]); 0; 0; :time];
        if config.dimension == 2
            cargs = [:(fnodes[coefi,1]); :(fnodes[coefi,2]); 0; :time];
        elseif config.dimension == 3
            cargs = [:(fnodes[coefi,1]); :(fnodes[coefi,2]); :(fnodes[coefi,3]); :time];
        end
        
        for i=1:length(needed_coef)
            if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt || needed_coef[i] === :DGNORMAL || needed_coef[i] === :DGNORMAL1 || needed_coef[i] === :DGNORMAL2)
                cind = get_coef_index(needed_coef[i]);
                if cind >= 0
                    tag = string(cind);
                else
                    tag = string(needed_coef[i]);
                end
                #derivatives of coefficients
                tag = needed_coef_deriv[i][2] * tag;
                #tmps = "coef_"*tag*"_"*string(needed_coef_ind[i]);
                tmps1 = "coef_DGSIDE1"*tag*"_"*string(needed_coef_ind[i]);
                tmps2 = "coef_DGSIDE2"*tag*"_"*string(needed_coef_ind[i]);
                #tmpc = Symbol(tmps);
                tmpc1 = Symbol(tmps1);
                tmpc2 = Symbol(tmps2);
                
                (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
                if ctype == 1
                    # constant coefficient -> coef_n = cval
                    tmpn = cval;
                    push!(code.args, Expr(:(=), tmpc1, tmpn));
                    push!(code.args, Expr(:(=), tmpc2, tmpn));
                    #push!(code.args, Expr(:(=), tmpc, tmpn));
                elseif ctype == 2
                    # genfunction coefficients -> coef_n_i = coef.value[i].func(cargs)
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc1;
                    tmpn = :(Femshop.genfunctions[$cval]); # Femshop.genfunctions[cval]
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    push!(code.args, Expr(:(=), tmpc1, :(zeros(refel.Nfp[frefelind[1]])))); # allocate coef_n
                    push!(code.args, Expr(:(=), tmpc2, :(zeros(refel.Nfp[frefelind[1]]))));
                    #push!(code.args, Expr(:(=), tmpc, tmpc1));
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                    
                elseif ctype == 3
                    tag = string(needed_coef[i]);
                    tmps1 = "coef_DGSIDE1"*tag*"_"*string(needed_coef_ind[i]);
                    tmps2 = "coef_DGSIDE2"*tag*"_"*string(needed_coef_ind[i]);
                    tmpc1 = Symbol(tmps1);
                    tmpc2 = Symbol(tmps2);
                    # variable values -> coef_n = variable.values
                    if variables[cval].type == SCALAR
                        tmpb1 = :(copy(Femshop.variables[$cval].values[loc2glb[1]]));
                        tmpb2 = :(copy(Femshop.variables[$cval].values[loc2glb[2]]));
                    else
                        compo = needed_coef_ind[i];
                        tmpb1 = :(copy(Femshop.variables[$cval].values[$compo, loc2glb[1]]));
                        tmpb2 = :(copy(Femshop.variables[$cval].values[$compo, loc2glb[2]]));
                    end
                    tmpindex = needed_coef_ind[i];
                    
                    push!(code.args, Expr(:(=), tmpc1, tmpb1));
                    push!(code.args, Expr(:(=), tmpc2, tmpb2));
                    
                    # bdrycond = :(if faceBID > 0 statement end);
                    # bdrycondin = Expr(:block);
                    # vind = 0;
                    # for varind=1:length(variables)
                    #     if variables[varind].symbol === needed_coef[i]
                    #         vind = variables[varind].index;
                    #         break;
                    #     end
                    # end
                    # push!(bdrycondin.args, Expr(:(=), tmpcave, :(Femshop.DGSolver.face_boundary_condition(Femshop.variables[$vind], faceBID, face2glb[:,1], fnodes, time))));
                    # bdrycond.args[2] = bdrycondin;
                    # push!(code.args, bdrycond);
                end
                
            end# number?
        end# coef loop
        
        # Write the loop that computes coefficient values
        if length(cloopin.args) > 0
            cloop.args[2] = cloopin;
            push!(code.args, cloop); # add loop to code
        end
        
        # Apply derivatives and do ave and jump after the initializing loop
        # Derivs will look like: coef_i_j = RDn * coef_i_j
        # Ave will look like: coef_DGAVEu_j = (coef_u_j_s1 + coef_u_j_s2).*0.5
        for i=1:length(needed_coef_deriv)
            if length(needed_coef_deriv[i][2]) > 0
                cind = get_coef_index(needed_coef[i]);
                
                if cind >= 0
                    tag = string(cind);
                else
                    tag = string(needed_coef[i]);
                end
                #modifications of coefficients
                tag = needed_coef_deriv[i][2] * tag;
                #tmps = "coef_"*tag*"_"*string(needed_coef_ind[i]);
                tmps1 = "coef_DGSIDE1"*tag*"_"*string(needed_coef_ind[i]);
                tmps2 = "coef_DGSIDE2"*tag*"_"*string(needed_coef_ind[i]);
                #tmpc = Symbol(tmps);
                tmpc1 = Symbol(tmps1);
                tmpc2 = Symbol(tmps2);
                
                if length(needed_coef_deriv[i][3]) > 0 && !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                    
                    
                    dmat1 = Symbol("RDL"*needed_coef_deriv[i][3]);
                    dmat2 = Symbol("RDR"*needed_coef_deriv[i][3]);
                    tmpb1= :(length($tmpc1) == 1 ? 0 : $dmat1 * $tmpc1);
                    tmpb2= :(length($tmpc2) == 1 ? 0 : $dmat2 * $tmpc2);
                    
                    #push!(code.args, :(println( $tmpc )));
                    push!(code.args, Expr(:(=), tmpc1, tmpb1));
                    push!(code.args, Expr(:(=), tmpc2, tmpb2));
                    #push!(code.args, Expr(:(=), tmpc, tmpb1)); # ?? This is not useful. surface terms should not have this
                    #push!(code.args, :(println( $tmpc )));
                    #push!(code.args, :(println("-----------------------------------------------------------------------------")));
                    
                else  # Could be some DG thing. TODO
                    # if occursin("DGAVE", needed_coef_deriv[i][2])
                    #     # {{u}} -> (fmap_s1(Q) + fmap_s2(Q))*0.5
                    #     #tmpb = :((refel.Q[fmap_s1,:]*$tmpc1 + refel.Q[fmap_s2,:]*$tmpc2) .* 0.5);
                    #     # tmpb1 = :((refel.Q[fmap_s1,:]*$tmpc1) .* 0.5);
                    #     # tmpb2 = :((refel.Q[fmap_s2,:]*$tmpc2) .* 0.5);
                    #     # tmpb1 = :($tmpc1 .* 0.5);
                    #     # tmpb2 = :($tmpc2 .* 0.5);
                    # elseif occursin("DGJUMP", needed_coef_deriv[i][2])
                    #     # [[u]] -> (fmap_s1(Q) - fmap_s2(Q))
                    #     #tmpb = :((refel.Q[fmap_s1,:]*$tmpc1 - refel.Q[fmap_s2,:]*$tmpc2));
                    #     # tmpb1 = :((refel.Q[fmap_s1,:]*$tmpc1));
                    #     # tmpb2 = :((-refel.Q[fmap_s2,:]*$tmpc2));
                    #     # tmpb1 = :($tmpc1);
                    #     # tmpb2 = :(-$tmpc2);
                        
                    # else
                    #     # tmpb1 = :(refel.Q[fmap_s1,:]*$tmpc1); # This should not happen?
                    #     # tmpb2 = :(refel.Q[fmap_s2,:]*$tmpc2);
                    #     # tmpb1 = :($tmpc1); # This should not happen?
                    #     # tmpb2 = :($tmpc2);
                    # end
                    
                    # #push!(code.args, Expr(:(=), tmpc, tmpb));
                    # # push!(code.args, Expr(:(=), tmpc1, tmpb1));
                    # # push!(code.args, Expr(:(=), tmpc2, tmpb2));
                    # # push!(code.args, :(println("tmpc1="*string($tmpc1))));
                    # # push!(code.args, :(println("tmpc2="*string($tmpc2))));
                end
                
            else
                # nothing?
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
    
    #if dofsper > 1
        if lorr == RHS
            push!(code.args, Expr(:(=), :element_vectorL, :(zeros(refel.Np*$dofsper)))); # allocate vector
            push!(code.args, Expr(:(=), :element_vectorR, :(zeros(refel.Np*$dofsper))));
        elseif dofsper > 1
            push!(code.args, Expr(:(=), :element_matrixLL, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper)))); # allocate matrix
            push!(code.args, Expr(:(=), :element_matrixLR, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper))));
            push!(code.args, Expr(:(=), :element_matrixRL, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper))));
            push!(code.args, Expr(:(=), :element_matrixRR, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper))));
        end
    #end
    
    # remove any :EMPTY terms
    newterms = [[], [], [], []];
    for i=1:length(terms[1])
        for tind=1:4
            if !(terms[tind][i] === nothing || terms[tind][i] === :EMPTY)
                push!(newterms[tind], terms[tind][i]);
            end
        end
    end
    terms = newterms;
    
    # If it was empty, just return zeros without doing any work
    if length(terms[1])==0 && length(terms[2])==0 && length(terms[3])==0 && length(terms[4])==0
        if lorr == LHS
            code = Expr(:block);
            push!(code.args, :(A=zeros(args[2].Np*$dofsper, args[2].Np*$dofsper)));
            push!(code.args, :(return [A, A, A, A]));
            return code;
        else
            code = Expr(:block);
            push!(code.args, :(b=zeros(args[2].Np*$dofsper)));
            push!(code.args, :(return [b, b]));
            return code;
        end
    end
    
    result = Array{Any,1}(undef,4); # Will hold the returned expression
    for i=1:4
        result[i] = nothing;
    end
    if typeof(var) <: Array # multivar
        println("multivar not ready");
        # for vi=1:length(var)
        #     # Add terms into full matrix according to testind/trialind
        #     # Each component/dof should have one expression so that submatrix is only modified once.
        #     if lorr == LHS
        #         #comps = length(var.symvar.vals);
        #         comps = dofsper;
        #         submatrices = Array{Any,2}(undef, comps, comps);
        #         for smi=1:length(submatrices)
        #             submatrices[smi] = nothing;
        #         end
        #         for i=1:length(terms[vi])
        #             ti = test_ind[vi][i][1] + offset_ind[vi];
        #             tj = trial_ind[vi][i][1];
                    
        #             if submatrices[ti, tj] === nothing
        #                 submatrices[ti, tj] = terms[vi][i];
        #             else
        #                 addexpr = :(a+b);
        #                 addexpr.args[2] = submatrices[ti, tj];
        #                 addexpr.args[3] = terms[vi][i];
        #                 submatrices[ti, tj] = addexpr;
        #             end
        #         end
                
        #         for cj=1:comps
        #             for ci=1:comps
        #                 if !(submatrices[ci, cj] === nothing)
        #                     ti = ci-1;
        #                     tj = cj-1;
        #                     sti = :($ti*refel.Np + 1);
        #                     eni = :(($ti + 1)*refel.Np);
        #                     stj = :($tj*refel.Np + 1);
        #                     enj = :(($tj + 1)*refel.Np);
                            
        #                     push!(code.args, Expr(:(+=), :(element_matrix[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
        #                 end
        #             end
        #         end
                
        #         result = :element_matrix;
                
        #     else #RHS
        #         #comps = length(var.symvar.vals);
        #         comps = dofsper;
        #         submatrices = Array{Any,1}(undef, comps);
        #         for smi=1:length(submatrices)
        #             submatrices[smi] = nothing;
        #         end
        #         for i=1:length(terms[vi])
        #             ti = test_ind[vi][i][1] + offset_ind[vi];
                    
        #             if submatrices[ti] === nothing
        #                 submatrices[ti] = terms[vi][i];
        #             else
        #                 addexpr = :(a+b);
        #                 addexpr.args[2] = submatrices[ti];
        #                 addexpr.args[3] = terms[vi][i];
        #                 submatrices[ti] = addexpr;
        #             end
        #         end
                
        #         for ci=1:comps
        #             if !(submatrices[ci] === nothing)
        #                 ti = ci-1;
        #                 sti = :($ti*refel.Np + 1);
        #                 eni = :(($ti + 1)*refel.Np);
                        
        #                 push!(code.args, Expr(:(+=), :(element_vector[$sti:$eni]), submatrices[ci]));
        #             end
        #         end
                
        #         result = :element_vector;
                
        #     end
        # end
    else
        for termind=1:4
            if length(terms[termind]) > 1
                if var.type == SCALAR # Only one component
                    tmp1 = :(a.+b);
                    tmp1.args = [:.+];
                    for i=1:length(terms[termind])
                        push!(tmp1.args, terms[termind][i]);
                    end
                    if lorr == LHS
                        if termind == 1
                            push!(code.args, Expr(:(=), :(element_matrixLL), tmp1));
                            result[1] = :element_matrixLL;
                        elseif termind == 2
                            push!(code.args, Expr(:(=), :(element_matrixLR), tmp1));
                            result[2] = :element_matrixLR;
                        elseif termind == 3
                            push!(code.args, Expr(:(=), :(element_matrixRL), tmp1));
                            result[3] = :element_matrixRL;
                        elseif termind == 4
                            push!(code.args, Expr(:(=), :(element_matrixRR), tmp1));
                            result[4] = :element_matrixRR;
                        end
                    else
                        if termind == 1
                            push!(code.args, Expr(:(.+=), :(element_vectorL), tmp1));
                            result[1] = :element_vectorL;
                        elseif termind == 2
                            push!(code.args, Expr(:(.+=), :(element_vectorL), tmp1));
                        elseif termind == 3
                            push!(code.args, Expr(:(.+=), :(element_vectorR), tmp1));
                        elseif termind == 4
                            push!(code.args, Expr(:(.+=), :(element_vectorR), tmp1));
                            result[4] = :element_vectorR;
                        end
                    end
                    
                else # More than one component
                    # Add terms into full matrix according to testind/trialind
                    # Each component/dof should have one expression so that submatrix is only modified once.
                    if lorr == LHS
                        comps = length(var.symvar.vals);
                        submatrices = Array{Any,2}(undef, comps, comps);
                        for smi=1:length(submatrices)
                            submatrices[smi] = nothing;
                        end
                        for i=1:length(terms[termind])
                            ti = test_ind[i][1];
                            tj = trial_ind[i][1];
                            
                            if submatrices[ti, tj] === nothing
                                submatrices[ti, tj] = terms[termind][i];
                            else
                                addexpr = :(a.+b);
                                addexpr.args[2] = submatrices[ti, tj];
                                addexpr.args[3] = terms[termind][i];
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
                                    
                                    if termind == 1
                                        push!(code.args, Expr(:(+=), :(element_matrixLL[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                                    elseif termind == 2
                                        push!(code.args, Expr(:(+=), :(element_matrixLR[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                                    elseif termind == 3
                                        push!(code.args, Expr(:(+=), :(element_matrixRL[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                                    elseif termind == 4
                                        push!(code.args, Expr(:(+=), :(element_matrixRR[$sti:$eni, $stj:$enj]), submatrices[ci, cj]));
                                    end
                                end
                            end
                        end
                        
                        if termind == 1
                            result[termind]= :element_matrixLL;
                        elseif termind == 2
                            result[termind]= :element_matrixLR;
                        elseif termind == 3
                            result[termind]= :element_matrixRL;
                        elseif termind == 4
                            result[termind]= :element_matrixRR;
                        end
                        
                    else #RHS
                        comps = length(var.symvar.vals);
                        submatrices = Array{Any,1}(undef, comps);
                        for smi=1:length(submatrices)
                            submatrices[smi] = nothing;
                        end
                        for i=1:length(terms[termind])
                            ti = test_ind[i][1];
                            
                            if submatrices[ti] === nothing
                                submatrices[ti] = terms[termind][i];
                            else
                                addexpr = :(a.+b);
                                addexpr.args[2] = submatrices[ti];
                                addexpr.args[3] = terms[termind][i];
                                submatrices[ti] = addexpr;
                            end
                        end
                        
                        for ci=1:comps
                            if !(submatrices[ci] === nothing)
                                ti = ci-1;
                                sti = :($ti*refel.Np + 1);
                                eni = :(($ti + 1)*refel.Np);
                                
                                push!(code.args, Expr(:(+=), :(element_vector_s1[$sti:$eni]), submatrices[ci]));
                                if termind == 1
                                    push!(code.args, Expr(:(+=), :(element_vectorL[$sti:$eni]), submatrices[ci]));
                                elseif termind == 4
                                    push!(code.args, Expr(:(+=), :(element_vectorR[$sti:$eni]), submatrices[ci]));
                                end
                            end
                        end
                        
                        if termind == 1
                            result[termind] = :element_vectorL;
                        elseif termind == 4
                            result[termind] = :element_vectorR;
                        end
                        
                    end
                end
            elseif length(terms[termind]) == 1# one term (one variable)
                #result[termind] = terms[termind][1];
                if lorr == LHS
                    if termind == 1
                        push!(code.args, Expr(:(=), :(element_matrixLL), terms[1][1]));
                        result[1] = :element_matrixLL;
                    elseif termind == 2
                        push!(code.args, Expr(:(=), :(element_matrixLR), terms[2][1]));
                        result[2] = :element_matrixLR;
                    elseif termind == 3
                        push!(code.args, Expr(:(=), :(element_matrixRL), terms[3][1]));
                        result[3] = :element_matrixRL;
                    elseif termind == 4
                        push!(code.args, Expr(:(=), :(element_matrixRR), terms[4][1]));
                        result[4] = :element_matrixRR;
                    end
                    
                else
                    if termind == 1
                        push!(code.args, Expr(:(=), :(element_vectorL), terms[1][1]));
                        result[1] = :element_vectorL;
                    elseif termind == 4
                        push!(code.args, Expr(:(=), :(element_vectorR), terms[4][1]));
                        result[4] = :element_vectorR;
                    end
                    
                end
            end
        end
        
    end
    if lorr == LHS
        if (result[1] === nothing || result[1] == 0) result[1] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper,refel.Nfp[frefelind[1]]*$dofsper)); end
        if (result[2] === nothing || result[2] == 0) result[2] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper,refel.Nfp[frefelind[1]]*$dofsper)); end
        if (result[3] === nothing || result[3] == 0) result[3] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper,refel.Nfp[frefelind[1]]*$dofsper)); end
        if (result[4] === nothing || result[4] == 0) result[4] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper,refel.Nfp[frefelind[1]]*$dofsper)); end
        
        r1 = result[1];
        r2 = result[2];
        r3 = result[3];
        r4 = result[4];
        finalresult = :([$r1, $r2, $r3, $r4]);
    else
        if (result[1] === nothing || result[1] == 0) result[1] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper)); end
        if (result[4] === nothing || result[4] == 0) result[4] = :(zeros(refel.Nfp[frefelind[1]]*$dofsper)); end
        
        r1 = result[1];
        r4 = result[4];
        finalresult = :([$r1, $r4]);
    end
    
    push!(code.args, Expr(:return, finalresult));
    return code;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_surface_term_julia(sterm, var, lorr, offset_ind=0)
    terms = Array{Any,1}(undef,4);
    terms[1] = copy(sterm);
    terms[2] = copy(sterm);
    terms[3] = copy(sterm);
    terms[4] = copy(sterm);
    #terms = [copy(sterm),copy(sterm),copy(sterm),copy(sterm)];
    need_derivative = false;
    need_derivative_for_coefficient = false;
    needed_coef = [];
    needed_coef_ind = [];
    needed_coef_deriv = [];
    
    # test_parts = [nothing,nothing,nothing,nothing];
    # coef_parts = [nothing,nothing,nothing,nothing];
    # trial_parts = [nothing,nothing,nothing,nothing];
    test_parts = Array{Any,1}(undef,4);
    coef_parts = Array{Any,1}(undef,4);
    trial_parts = Array{Any,1}(undef,4);
    for i=1:4
        test_parts[i] = nothing;
        coef_parts[i] = nothing;
        trial_parts[i] = nothing;
    end
    weight_parts = [:face_wgdetj,:face_wgdetj,:face_wgdetj,:face_wgdetj];
    test_component = 0;
    trial_component = 0;
    
    # extract each of the factors.
    factors = separate_factors(terms[1]);
    
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
            if factors[i].args[1] === :./ || factors[i].args[1] === :/
                # The second arg should be 1, the third should not contain an unknown or test symbol
                # The denominator expression needs to be processed completely
                (piece, nd, nc, nci, ncd) = process_known_expr_julia(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[3] = piece;
                push!(coef_expr_facs, factors[i]);
                #push!(coef_inds, 0);
                
            elseif factors[i].args[1] === :.^ || factors[i].args[1] === :^
                factors[i].args[1] = :.^ ;
                # The second arg is the thing raised
                (piece1, nd, nc, nci, ncd) = process_known_expr_julia(factors[i].args[2]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[2] = piece1;
                # Do the same for the power just in case
                (piece2, nd, nc, nci, ncd) = process_known_expr_julia(factors[i].args[3]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[3] = piece2;
                
                push!(coef_expr_facs, factors[i]);
                #push!(coef_inds, 0);
            elseif factors[i].args[1] === :sqrt
                factors[i].args[1] = :.^
                # The second arg is the thing sqrted
                (piece1, nd, nc, nci, ncd) = process_known_expr_julia(factors[i].args[2]);
                need_derivative = need_derivative || nd;
                append!(needed_coef, nc);
                append!(needed_coef_ind, nci);
                append!(needed_coef_deriv, ncd);
                factors[i].args[2] = piece1;
                # add a 1/2 power argument
                push!(factors[i].args, 1/2);
                
                push!(coef_expr_facs, factors[i]);
                #push!(coef_inds, 0);
            end
            
        else
            (index, v, mods) = extract_symbols(factors[i]);
            
            #println("factor: "*string(factors[i])*", mods: "*string(mods));
            
            if is_test_func(v)
                test_component = index; # the vector index
                if length(mods) > 0
                    # Could be Dn for derivative, DGJUMP for jump, DGAVE for ave or versions with NORMDOTGRAD
                    if occursin("DGSIDE1", mods[1]) || (length(mods)>1 && occursin("DGSIDE1", mods[2]))
                        if length(mods)>1 && occursin("D", mods[1])
                            # TODO more than one derivative mod
                            need_derivative = true;
                            tp1 = Symbol("TRQ1_"*mods[1][2]);
                            test_parts[1] = :($tp1);
                            test_parts[2] = :($tp1);
                            test_parts[3] = :EMPTY;
                            test_parts[4] = :EMPTY;
                        else
                            test_parts[1] = :(refel.surf_Q[frefelind[1]]');
                            test_parts[2] = :(refel.surf_Q[frefelind[1]]');
                            test_parts[3] = :EMPTY;
                            test_parts[4] = :EMPTY;
                        end
                        
                    elseif occursin("DGSIDE2", mods[1]) || (length(mods)>1 && occursin("DGSIDE2", mods[2]))
                        if length(mods)>1 && occursin("D", mods[1])
                            # TODO more than one derivative mod
                            need_derivative = true;
                            tp1 = Symbol("TRQ2_"*mods[1][2]);
                            test_parts[1] = :EMPTY;
                            test_parts[2] = :EMPTY;
                            test_parts[3] = :($tp1);
                            test_parts[4] = :($tp1);
                        else
                            test_parts[1] = :EMPTY;
                            test_parts[2] = :EMPTY;
                            test_parts[3] = :(refel.surf_Q[frefelind[2]]');
                            test_parts[4] = :(refel.surf_Q[frefelind[2]]');
                        end
                        
                    elseif occursin("D", mods[1])
                        # TODO more than one derivative mod
                        need_derivative = true;
                        test_parts[1] = Symbol("TRQ1_"*mods[1][2]);
                        test_parts[2] = Symbol("TRQ1_"*mods[1][2]);
                        test_parts[3] = Symbol("TRQ2_"*mods[1][2]);
                        test_parts[4] = Symbol("TRQ2_"*mods[1][2]);
                    end
                else
                    # no mods
                    test_parts[1] = :(refel.surf_Q[frefelind[1]]');
                    test_parts[2] = :(refel.surf_Q[frefelind[1]]');
                    test_parts[3] = :(refel.surf_Q[frefelind[2]]');
                    test_parts[4] = :(refel.surf_Q[frefelind[2]]');
                end
            elseif is_unknown_var(v, var) && lorr == LHS # If rhs, treat as a coefficient
                # if !(trial_parts[1] === nothing)
                #     # Two unknowns multiplied in this term. Nonlinear. abort.
                #     printerr("Nonlinear term. Code layer incomplete.");
                #     return (-1, -1, -1, -1, -1, -1);
                # end
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
                    # Could be Dn for derivative, DGJUMP for jump, DGAVE for ave or versions with NORMDOTGRAD
                    if occursin("DGSIDE1", mods[1]) || (length(mods)>1 && occursin("DGSIDE1", mods[2]))
                        if length(mods)>1 && occursin("D", mods[1])
                            # TODO more than one derivative mod
                            need_derivative = true;
                            tp1 = Symbol("RQ1_"*mods[1][2]);
                            trial_parts[1] = :($tp1);
                            trial_parts[2] = :EMPTY;
                            trial_parts[3] = :($tp1);
                            trial_parts[4] = :EMPTY;
                            
                        else
                            trial_parts[1] = :(refel.surf_Q[frefelind[1]]);
                            trial_parts[2] = :EMPTY;
                            trial_parts[3] = :(refel.surf_Q[frefelind[1]]);
                            trial_parts[4] = :EMPTY;
                        end
                        
                    elseif occursin("DGSIDE2", mods[1]) || (length(mods)>1 && occursin("DGSIDE2", mods[2]))
                        if length(mods)>1 && occursin("D", mods[1])
                            # TODO more than one derivative mod
                            need_derivative = true;
                            tp1 = Symbol("RQ2_"*mods[1][2]);
                            trial_parts[1] = :EMPTY;
                            trial_parts[2] = :($tp1);
                            trial_parts[3] = :EMPTY;
                            trial_parts[4] = :($tp1);
                            
                        else
                            trial_parts[1] = :EMPTY;
                            trial_parts[2] = :(refel.surf_Q[frefelind[2]]);
                            trial_parts[3] = :EMPTY;
                            trial_parts[4] = :(refel.surf_Q[frefelind[2]]);
                        end
                        
                    elseif occursin("D", mods[1])
                        # TODO more than one derivative mod
                        need_derivative = true;
                        trial_parts[1] = Symbol("RQ1_"*mods[1][2]);
                        trial_parts[2] = :EMPTY;
                        trial_parts[3] = :EMPTY;
                        trial_parts[4] = Symbol("RQ2_"*mods[1][2]);
                    end
                else
                    # no mods
                    trial_parts[1] = :(refel.surf_Q[frefelind[1]]);
                    trial_parts[2] = :EMPTY;
                    trial_parts[3] = :EMPTY;
                    trial_parts[4] = :(refel.surf_Q[frefelind[2]]);
                end
                
            else # coefficients
                if length(index) == 1
                    ind = index[1];
                end
                # Check for derivative mods
                if typeof(v) == Symbol && !(v ===:dt)
                    if length(mods) > 0
                        # Could be Dn for derivative, DGJUMP for jump, DGAVE for ave or versions with NORMDOTGRAD
                        if occursin("DGSIDE1", mods[1])
                            # {{n.grad(u)}}
                            push!(needed_coef_deriv, [v, mods[1], ""]);
                        elseif occursin("DGSIDE2", mods[1])
                            # {{n.grad(u)}}
                            push!(needed_coef_deriv, [v, mods[1], ""]);
                        elseif occursin("D", mods[1])
                            need_derivative = true;
                            need_derivative_for_coefficient = false;
                            
                            push!(needed_coef_deriv, [v, mods[1], mods[1][2]]);
                        end
                        
                    else
                        push!(needed_coef_deriv, [v, "", ""]);
                    end
                    push!(needed_coef, v);
                    push!(needed_coef_ind, ind);
                    
                    push!(coef_derivs, needed_coef_deriv[end]);
                    
                else
                    push!(coef_derivs, [nothing, "", ""]);
                end
                
                push!(coef_facs, v);
                push!(coef_inds, ind);
            end
        end
        
    end # factors loop
    
    # If there's no trial part, need to do this
    # if trial_part === nothing
    #     trial_part = :(refel.Q[fmap_s1,:]);
    # end
    
    # build coefficient parts
    if length(coef_facs) > 0
        for j=1:length(coef_facs)
            tmp1 = coef_facs[j];
            tmp2 = coef_facs[j];
            tmp3 = coef_facs[j];
            tmp4 = coef_facs[j];
            #println("coef_facs: "*string(tmp)*" : "*string(typeof(tmp)));
            if typeof(tmp1) == Symbol && !(tmp1 ===:dt)
                ind = get_coef_index(coef_facs[j]);
                if ind >= 0
                    tag = string(ind);
                else
                    tag = string(tmp1);
                end
                if tmp1 === :DGNORMAL
                    ind = coef_inds[j];
                    tmp1 = :(normal[$ind]);
                    tmp2 = :(-normal[$ind]);
                    tmp3 = :(normal[$ind]);
                    tmp4 = :(-normal[$ind]);
                elseif tmp1 === :DGNORMAL1
                    ind = coef_inds[j];
                    tmp1 = :(normal[$ind]);
                    tmp2 = :(normal[$ind]);
                    tmp3 = :(normal[$ind]);
                    tmp4 = :(normal[$ind]);
                elseif tmp1 === :DGNORMAL2
                    ind = coef_inds[j];
                    tmp1 = :(-normal[$ind]);
                    tmp2 = :(-normal[$ind]);
                    tmp3 = :(-normal[$ind]);
                    tmp4 = :(-normal[$ind]);
                else
                    #derivatives of coefficients
                    tag = coef_derivs[j][2] * tag;
                    if occursin("DGSIDE1", tag)
                        tmps1 = "coef_"*tag*"_"*string(coef_inds[j]);
                        tmps2 = "EMPTY";
                        tmps3 = "coef_"*tag*"_"*string(coef_inds[j]);
                        tmps4 = "EMPTY";
                    elseif occursin("DGSIDE2", tag)
                        tmps1 = "EMPTY";
                        tmps2 = "coef_"*tag*"_"*string(coef_inds[j]);
                        tmps3 = "EMPTY";
                        tmps4 = "coef_"*tag*"_"*string(coef_inds[j]);
                    else
                        tmps1 = "coef_DGSIDE1"*tag*"_"*string(coef_inds[j]);
                        tmps2 = "coef_DGSIDE2"*tag*"_"*string(coef_inds[j]);
                        tmps3 = "coef_DGSIDE1"*tag*"_"*string(coef_inds[j]);
                        tmps4 = "coef_DGSIDE2"*tag*"_"*string(coef_inds[j]);
                    end
                    tmp1 = Symbol(tmps1);
                    tmp2 = Symbol(tmps2);
                    tmp3 = Symbol(tmps3);
                    tmp4 = Symbol(tmps4);
                end
            end
            if j>1
                c1=coef_parts[1];
                c2=coef_parts[2];
                c3=coef_parts[3];
                c4=coef_parts[4];
                if (tmp1 === :EMPTY)||(c1 === :EMPTY) 
                    coef_parts[1] = :EMPTY;
                else
                    coef_parts[1] = :($c1 .* $tmp1);
                end
                if (tmp2 === :EMPTY)||(c2 === :EMPTY) 
                    coef_parts[2] = :EMPTY;
                else
                    coef_parts[2] = :($c2 .* $tmp2);
                end
                if (tmp3 === :EMPTY)||(c3 === :EMPTY) 
                    coef_parts[3] = :EMPTY;
                else
                    coef_parts[3] = :($c3 .* $tmp3);
                end
                if (tmp4 === :EMPTY)||(c4 === :EMPTY) 
                    coef_parts[4] = :EMPTY;
                else
                    coef_parts[4] = :($c4 .* $tmp4);
                end
                # coef_parts[2] = :($c2 .* $tmp2);
                # coef_parts[3] = :($c3 .* $tmp1);
                # coef_parts[4] = :($c4 .* $tmp2);
            else
                coef_parts[1] = tmp1;
                coef_parts[2] = tmp2;
                coef_parts[3] = tmp3;
                coef_parts[4] = tmp4;
            end
        end
    end
    
    if length(coef_expr_facs) > 0
        for j=1:length(coef_expr_facs)
            tmpc = coef_expr_facs[j];
            for cind=1:4
                if !(coef_parts[cind] === nothing)
                    cp = coef_parts[cind];
                    coef_parts[cind] = :($cp .* $tmpc);
                else
                    coef_parts[cind] = tmpc;
                end
            end
        end
    end
    
    for tind=1:4
        # If there's no trial part, need to do this
        if trial_parts[tind] === nothing
            if tind == 1 || tind == 3
                trial_parts[tind] = :(refel.surf_Q[frefelind[1]]);
            else
                trial_parts[tind] = :(refel.surf_Q[frefelind[2]]);
            end
        end
        
        # If there's no test part this is probably a denominator expression being processed and should only contain coefficients/knowns
        if test_parts[tind] === nothing
            #
            #
        elseif test_parts[tind] === :EMPTY || trial_parts[tind] === :EMPTY
            # This term will be zero, so just make it EMPTY
            terms[tind] = :EMPTY;
        else
            terms[tind] = test_parts[tind];
            tp = test_parts[tind];
            wp = weight_parts[tind];
            trp = trial_parts[tind];
            if !(coef_parts[tind] === nothing)
                cp = coef_parts[tind];
                if lorr == LHS
                    terms[tind] = :($tp * (diagm($wp .* $cp) * $trp));
                else # RHS
                    if cp === :EMPTY
                        terms[tind] = :EMPTY;
                    else
                        terms[tind] = :($tp * ($wp .* ($trp * $cp)));
                    end
                end
                
            else
                terms[tind]= :($tp * diagm($wp) * $trp);
            end
            
            if neg
                negex = :(-a);
                negex.args[2] = copy(terms[tind]);
                terms[tind] = negex;
            end
        end
    end
    
    return (terms, need_derivative, needed_coef, needed_coef_ind, needed_coef_deriv, test_component, trial_component);
end
