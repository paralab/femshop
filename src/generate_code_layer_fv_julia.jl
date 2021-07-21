#=
Code generation functions for FV for Julia.
=#

# Extract the input args. This must match the arguments passed by the solver.
function handle_input_args_fv_julia(lorr, vors)
    code = "";
    
    if vors == "volume" # source (integrated over the cell volume)
        code =
"
var =    args[1];  # unknown variables
el =     args[2];  # This element index
nodex =  args[3];  # global coords of element's nodes
loc2glb= args[4];  # global indices of the nodes
refel =  args[5];  # reference element
detj =   args[6];  # quadrature weights * detJ
J =      args[7];  # Jacobian
time =   args[8];  # time for time dependent coefficients
dt =     args[9];  # dt for time dependent problems
"
        
    elseif vors == "surface" # flux (integrated over the cell surface)
        code =
"
var =        args[1];  # list of unknown variables for this expression
els =        args[2];  # elements on both sides of the face
refel =      args[3];  # reference element for volume
loc2glb =    args[4];  # loc2glb for full elements
nodex =      args[5];  # global coords of element's nodes
cellx =      args[6];  # coordinates of cell center
frefelind =  args[7];  # reference element indices for faces
facex =      args[8];  # global coords of face nodes
face2glb =   args[9];  # Global indices of face nodes
normal =     args[10]; # Normal vector from e1 to e2
face_detJ =  args[11]; # geometric factor for face
area =       args[12]; # area face
pvol_J =     args[13]; # jacobian for both elements
ptime =      args[14]; # time for time dependent coefficients
dt =         args[15]; # dt for time dependent problems
"
    end
    
    return code;
end

# If needed, build derivative matrices
function build_derivative_matrices_fv_julia(lorr, vors, need_matrix)
    if need_matrix
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
                #TODO
            end
            
        elseif config.dimension == 2
            if vors == "volume"
                code *= "(RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);\n";
                code *= "(TRQ1, TRQ2) = (RQ1', RQ2');\n"
            else
                #TODO
            end
            
        elseif config.dimension == 3
            if vors == "volume"
                code *= "(RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);\n";
                code *= "(TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');\n"
            else
                #TODO
            end
        end
        
    else
        code = ""
    end
    
    # This is computed for all
    if vors == "surface"
        code *= "dxyz = norm(cellx[2] - cellx[1]) .* normal; # normal scaled by distance between cell centers\n"
    end
    
    return code;
end

# Allocate, compute, or fetch all needed values
function prepare_needed_values_fv_julia(entities, var, lorr, vors)
    code = "";
    needs_integration = [];
    for i=1:length(entities)
        cname = make_coef_name(entities[i]);
        if is_unknown_var(entities[i], var) && lorr == LHS
            # TODO
            # if length(entities[i].derivs) > 0
            #     xyzchar = ["x","y","z"];
            #     for di=1:length(entities[i].derivs)
            #         code *= cname * " = RQ"*string(entities[i].derivs[di])*"; # d/d"*xyzchar[entities[i].derivs[di]]*" of trial function\n";
            #     end
            # else
            #     code *= cname * " = refel.Q; # trial function.\n";
            # end
        else
            # Is coefficient(number or function) or variable(array)?
            (ctype, cval) = get_coef_val(entities[i]);
            if ctype == -1
                # It was a special symbol like dt or FACENORMAL
                if entities[i].name == "FACENORMAL1"
                    code *= cname * " = normal["*string(entities[i].index)*"]; # normal vector component\n"
                else entities[i].name == "FACENORMAL2"
                    code *= cname * " = -normal["*string(entities[i].index)*"]; # reverse normal vector component\n"
                end
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
                # coef_n_i = zeros(refel.Np); # allocate
                # for coefi = 1:refel.Np
                #     coef_k_i[coefi] = (Femshop.genfunctions[cval]).func(x[1,coefi], x[2,coefi],x[3,coefi],time); # evaluate at nodes
                # end
                ######################################
                if vors == "surface" && length(entities[i].derivs) == 0
                    nodesymbol = "facex"
                else
                    nodesymbol = "nodex"
                end
                cargs = "("*nodesymbol*"[coefi], 0, 0, time)";
                if config.dimension == 2
                    cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], 0, time)";
                elseif config.dimension == 3
                    cargs = "("*nodesymbol*"[1, coefi], "*nodesymbol*"[2, coefi], "*nodesymbol*"[3, coefi], time)";
                end
                if vors == "volume"
                    code *= cname * " = zeros(refel.Np);\n";
                    code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Femshop.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                    end
                    # integrate over cell
                    code *= cname * " = (refel.wg .* detj)' * refel.Q * " * cname * "; # integrate over cell\n";
                    
                else # surface
                    # Apply any needed derivative operators. Interpolate at quadrature points.
                    if length(entities[i].derivs) > 0
                        code *= cname * " = zeros(refel.Np);\n";
                        code *= "for coefi = 1:refel.Np " * cname * "[coefi] = (Femshop.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                        xyzchar = ["x","y","z"];
                        for di=1:length(entities[i].derivs)
                            code *= cname * " = RD"*string(entities[i].derivs[di])*" * " * cname * 
                                    "; # Apply d/d"*xyzchar[entities[i].derivs[di]]*".\n";
                        end
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]; # extract face values only.";
                    else # no derivatives, only need surface nodes
                        code *= cname * " = zeros(refel.Nfp[frefelind[1]]);\n";
                        code *= "for coefi = 1:refel.Nfp[frefelind[1]] " * cname * "[coefi] = (Femshop.genfunctions["*string(cval)*"]).func" * cargs * " end\n";
                    end
                    # integrate over face
                    if config.dimension == 1
                        # in 1d there is only one face node
                        code *= cname * " = " * cname * "[1]\n";
                    else
                        code *= cname * " = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:, refel.face2local[frefelind[1]]] * " * cname * " / area; # integrate over face\n";
                    end
                end
                
            elseif ctype == 3 # a known variable value
                # This generates something like: coef_u_1 = copy((Femshop.variables[1]).values[1, gbl])
                dgside = 0; # 0 means no side flag
                for flagi=1:length(entities[i].flags)
                    if occursin("DGSIDE1", entities[i].flags[flagi])
                        dgside = 1;
                    elseif occursin("DGSIDE2", entities[i].flags[flagi])
                        dgside = 2;
                    end
                end
                if vors == "surface"
                    l2gsymbol = "els[1]"
                else
                    l2gsymbol = "el"
                end
                if dgside == 1
                    l2gsymbol = "els[1]"
                elseif dgside == 2
                    l2gsymbol = "els[2]"
                end
                if vors == "volume"
                    # Apply any needed derivative operators.
                    if length(entities[i].derivs) > 0
                        # Need a derivative here...
                        # TODO
                    else
                        code *= cname * " = Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", "*l2gsymbol*"];\n";
                    end
                else
                    if length(entities[i].derivs) > 0
                        code *= cname * " = Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", els[2]] - Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", els[1]];\n";
                        code *= cname * " = (els[1] != els[2] && abs(normal["*string(entities[i].derivs[1])*"]) > 1e-10) ? "*cname*" ./ dxyz["*string(entities[i].derivs[1])*"]  : 0\n"
                    else
                        code *= cname * " = Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", "*l2gsymbol*"];\n";
                    end
                end
                
            end
        end # if coefficient
    end # entity loop
    
    return code;
end

function make_elemental_computation_fv_julia(terms, var, dofsper, offset_ind, lorr, vors)
    # Here is where I make some assumption about the form of the expression.
    # Since it was expanded by the parser it should look like a series of terms: t1 + t2 + t3...
    # Where each term, if LHS, involves one unknown component and possibly some coefficients.
    code = "";
    
    # If there were no terms, 
    if length(terms) < 1
        code = "return nothing # There were no terms to compute here";
    end
    
    # Allocate the vector or matrix to be returned if needed
    if dofsper > 1
        if lorr == RHS
            code *= "cell_average = zeros("*string(dofsper)*"); # Allocate for returned values.\n"
        else
            if vors == "volume"
                code *= "cell_matrix = zeros(refel.Np * "*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned matrix.\n"
            else
                code *= "cell_matrix = zeros(refel.Nfp[frefelind[1]] * "*string(dofsper)*", "*string(dofsper)*"); # Allocate for returned matrix.\n"
            end
        end
    end
    
    # Separate the factors of each term and form the calculation
    if dofsper > 1
        # Subvectors for each component
        if lorr == LHS
            subvector = Array{String, 2}(undef, dofsper, dofsper);
        else # RHS
            subvector = Array{String, 1}(undef, dofsper);
        end
        for smi=1:length(subvector)
            subvector[smi] = "";
        end
        
        if typeof(var) <: Array
            for vi=1:length(var) # variables
                # Process the terms for this variable
                for ci=1:length(terms[vi]) # components
                    for i=1:length(terms[vi][ci])
                        (term_result, var_ind) = generate_term_calculation_fv_julia(terms[vi][ci][i], var, lorr);
                        
                        # println(terms)
                        # println(terms[vi])
                        # println(terms[vi][ci])
                        # println(terms[vi][ci][i])
                        # println(term_result * " : "*string(test_ind)*", "*string(trial_ind))
                        
                        # Find the appropriate subvector for this term
                        subveci = offset_ind[vi] + ci;
                        subvecj = var_ind;
                        if lorr == LHS
                            subvec_ind = subveci + dofsper * (subvecj-1);
                        else
                            subvec_ind = subveci;
                        end
                        
                        if length(subvector[subvec_ind]) > 1
                            subvector[subvec_ind] *= " .+ " * term_result;
                        else
                            subvector[subvec_ind] = term_result;
                        end
                    end
                end
                
            end # vi
            
        else # only one variable
            # Process the terms for this variable
            for ci=1:length(terms) # components
                for i=1:length(terms[ci])
                    (term_result, var_ind) = generate_term_calculation_fv_julia(terms[ci][i], var, lorr);
                    
                    # Find the appropriate submatrix for this term
                    subvec_ind = ci;
                    
                    if length(subvector[subvec_ind]) > 1
                        subvector[subvec_ind] *= " .+ " * term_result;
                    else
                        subvector[subvec_ind] = term_result;
                    end
                end
            end
            
        end
        
        # Put the subvector together into cell_average or cell_matrix
        if lorr == LHS
            for emi=1:dofsper
                for emj=1:dofsper
                    if length(subvector[emi]) > 1
                        if vors == "volume"
                            rangei = "("*string(emi-1)*"*refel.Np + 1):("*string(emi)*"*refel.Np)";
                        else
                            rangei = "("*string(emi-1)*"*refel.Nfp[frefelind[1]] + 1):("*string(emi)*"*refel.Nfp[frefelind[1]])";
                        end
                        
                        code *= "cell_matrix["*rangei*"] = " * subvector[emi, emj] * "\n";
                    end
                end
            end
            code *= "return cell_matrix;\n"
            
        else # RHS
            for emi=1:dofsper
                if length(subvector[emi]) > 1
                    code *= "cell_average["*string(emi)*"] = " * subvector[emi] * "\n";
                end
            end
            code *= "return cell_average;\n"
        end
        
        
    else # one dof
        terms = terms[1];
        if lorr == LHS
            if vors == "volume"
                result = "zeros(refel.Np)";
            else
                result = "zeros(refel.Nfp[frefelind[1]])";
            end
            
        else
            result = "0";
        end
        
        #process each term
        first_term = true;
        for i=1:length(terms)
            (term_result, var_ind) = generate_term_calculation_fv_julia(terms[i], var, lorr);
            
            if !(term_result == "0")
                if first_term
                    result = term_result;
                    first_term = false;
                else
                    result *= " .+ " * term_result;
                end
            end
        end
        code *= "return " * result * ";\n";
    end
    
    return code;
end

function generate_term_calculation_fv_julia(term, var, lorr)
    result = "";
    # Note: separate_factors return test and trial info, but FV will not have any test functions.
    # trial_part refers to the unknown variable part.
    # example:
    # 0.1*_c_1*D1__u_1*_FACENORMAL1_1    ->  ???               on LHS
    #                             ->  0.1 * (coef_D1xu_1 .* normal[1])     on RHS
    if lorr == LHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term, var);
        # # LHS: test_part * diagm(weight_part .* coef_part) * trial_part
        # if !(coef_part === nothing)
        #     result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj .* " * 
        #             string(replace_entities_with_symbols(coef_part)) * ") * " * 
        #             string(replace_entities_with_symbols(trial_part));
        # else # no coef_part
        #     result = string(replace_entities_with_symbols(test_part)) * " * diagm(wdetj) * " * 
        #             string(replace_entities_with_symbols(trial_part));
        # end
    else
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term);
        # RHS: coef_part
        if !(coef_part === nothing)
            result = string(replace_entities_with_symbols(coef_part));
        else
            result = "0";
        end
    end
    
    return (result, var_ind);
end





####################################################################################################################
####################################################################################################################











# #=
# Use the symbolic layer expressions to generate the FV code.
# Instead of bilinear and linear expressions, these are flux and source expressions.
# LHS->flux, RHS->source
# =#

# function generate_code_layer_fv(ex, var, lorr, vors)
#     if use_cachesim
#         printerr("Cachesim not ready for FV.")
#     else
#         if language == 0 || language == JULIA
#             return generate_code_layer_fv_julia(ex, var, lorr, vors);
#         else
#             printerr("Exernal code gen not ready for FV.")
#         end
#     end
# end

# ###############################################################################################################
# # julia
# ###############################################################################################################

# #=
# Generates expressions to be turned into functions called in the elemental loop.
# symex is an array of SymEngine terms that are added together to form the expression.
# Translate them each into Julia expressions to be inserted in the code.

# - symex: symbolic expression(array of Basic)
# - var: variable or array of variables
# - lorr: LHS for flux, RHS for source
# =#
# function generate_code_layer_fv_julia(symex, var, lorr, vors)
#     # This is the basic info passed in "args"
#     code = Expr(:block);
#     if vors == "volume" # source (integrated over the cell volume)
#         # handle_args = quote
#         #     var =   args[1]; # unknown variables
#         #     el =    args[2]; # This element index
#         #     nodex = args[3]; # global coords of element's nodes
#         #     loc2glb=args[4]; # global indices of the nodes
#         #     refel = args[5]; # reference element
#         #     detj =  args[6]; # quadrature weights * detJ
#         #     J =     args[7]; # Jacobian
#         #     time =  args[8]; # time for time dependent coefficients
#         #     dt =    args[9]; # dt for time dependent problems
#         # end
#         push!(code.args, :(var =    args[1]));  # unknown variables
#         push!(code.args, :(el =     args[2]));  # This element index
#         push!(code.args, :(nodex =  args[3]));    # global coords of element's nodes
#         push!(code.args, :(loc2glb= args[4]));  # global indices of the nodes
#         push!(code.args, :(refel =  args[5]));# reference element
#         push!(code.args, :(detj =   args[6]));# quadrature weights * detJ
#         push!(code.args, :(J =      args[7]));    # Jacobian
#         push!(code.args, :(time =   args[8])); # time for time dependent coefficients
#         push!(code.args, :(dt =     args[9]));   # dt for time dependent problems
        
#         # # Quadrature matrix
#         # push!(code.args, :(Q = (refel.wg .* detj)' * refel.Q));
        
#     elseif vors == "surface" # flux (integrated over the cell surface)
#         # handle_args = quote
#         #     var =        args[1];  # list of unknown variables for this expression
#         #     els =        args[2];  # elements on both sides of the face
#         #     refel =      args[3];  # reference element for volume
#         #     loc2glb =    args[4];  # loc2glb for full elements
#         #     nodex =      args[5];  # global coords of element's nodes
#         #     cellx =      args[6];  # coordinates of cell center
#         #     frefelind =  args[7];  # reference element indices for faces
#         #     facex =      args[8];  # global coords of face nodes
#         #     face2glb =   args[9];  # Global indices of face nodes
#         #     normal =     args[10]; # Normal vector from e1 to e2
#         #     face_detJ =  args[11]; # geometric factor for face
#         #     vol_J =      args[12]; # jacobian for both elements
#         #     time =       args[13]; # time for time dependent coefficients
#         #     dt =         args[14]; # dt for time dependent problems
#         # end
#         push!(code.args, :(var =        args[1]));  # list of unknown variables for this expression
#         push!(code.args, :(els =        args[2]));  # elements on both sides of the face
#         push!(code.args, :(refel =      args[3]));  # reference element for volume
#         push!(code.args, :(loc2glb =    args[4]));  # loc2glb for full elements
#         push!(code.args, :(nodex =      args[5]));  # global coords of element's nodes
#         push!(code.args, :(cellx =      args[6]));  # coordinates of cell center
#         push!(code.args, :(frefelind =  args[7]));  # reference element indices for faces
#         push!(code.args, :(facex =      args[8]));  # global coords of face nodes
#         push!(code.args, :(face2glb =   args[9]));  # Global indices of face nodes
#         push!(code.args, :(normal =     args[10]));  # Normal vector from e1 to e2
#         push!(code.args, :(face_detJ =  args[11])); # geometric factor for face
#         push!(code.args, :(area =       args[12])); # area face
#         push!(code.args, :(vol_J =      args[13])); # jacobian for both elements
#         push!(code.args, :(time =       args[14])); # time for time dependent coefficients
#         push!(code.args, :(dt =         args[15])); # dt for time dependent problems
        
#         # # Quadrature matrix
#         # push!(code.args, :(surf_Q1 = (refel.surf_wg[frefelind[1]] .* face_detj)' * refel.surf_Q[frefelind[1]]));
#         # push!(code.args, :(surf_Q2 = (refel.surf_wg[frefelind[2]] .* face_detj)' * refel.surf_Q[frefelind[2]]));
#     end
    
#     # Keep track of what coefficients need to be evaluated and what derivatives need to be applied.
#     needed_coef = [];
#     needed_coef_name = [];
#     needed_coef_ind = [];
#     needed_coef_deriv = [];
    
#     var_ind = [];
    
#     # For multi variables, count them
#     multivar = typeof(var) <:Array;
#     varcount = 1;
#     if multivar
#         varcount = length(var);
#         offset_ind = zeros(Int, varcount);
#         tmp = length(var[1].symvar.vals);
#         for i=2:length(var)
#             offset_ind[i] = tmp;
#             tmp = tmp + length(var[i].symvar.vals);
#         end
#     end
    
#     # symex is an array of arrays of symengine terms, or array of arrays of arrays for multivar
#     # turn each one into an Expr for translation purposes
#     terms = [];
#     if multivar
#         for vi=1:varcount
#             push!(terms, terms_to_expr(symex[vi]));
#         end
#     else
#         terms = terms_to_expr(symex);
#     end
    
#     # Process the terms turning them into the code layer
#     if vors == "volume"
#         process_term_function = process_source_term_fv_julia;
#     else
#         process_term_function = process_flux_term_fv_julia;
#     end
#     if multivar
#         for vi=1:varcount
#             subvar_ind = [];
#             for i=1:length(terms[vi])
#                 (codeterm, coe, coeind, coenam, coederiv, varj) = process_term_function(terms[vi][i], var, lorr, offset_ind);
#                 if coeind == -1
#                     # processing failed
#                     printerr("term processing failed for: "*string(terms[vi][i]));
#                     return nothing;
#                 end
#                 append!(needed_coef, coe);
#                 append!(needed_coef_name, coenam);
#                 append!(needed_coef_ind, coeind);
#                 append!(needed_coef_deriv, coederiv);
                
#                 push!(subvar_ind, varj);
#                 terms[vi][i] = codeterm;
#             end
#             push!(var_ind, subvar_ind);
#         end
#     else
#         for i=1:length(terms)
#             (codeterm, coe, coeind, coenam, coederiv, varj) = process_term_function(terms[i], var, lorr, vors);
#             if coeind == -1
#                 # processing failed
#                 printerr("term processing failed for: "*string(terms[i]));
#                 return nothing;
#             end
#             append!(needed_coef, coe);
#             append!(needed_coef_name, coenam);
#             append!(needed_coef_ind, coeind);
#             append!(needed_coef_deriv, coederiv);
            
#             push!(var_ind, varj);
#             terms[i] = codeterm;
#         end
#     end
    
#     # remove any zero terms
#     if multivar
#         newterms = [];
#         for vi=1:varcount
#             push!(newterms, []);
#             for i=1:length(terms[vi])
#                 if !(terms[vi][i] == 0)
#                     push!(newterms[vi], terms[vi][i]);
#                 end
#             end
#         end
#     else
#         newterms = [];
#         for i=1:length(terms)
#             if !(terms[i] == 0)
#                 push!(newterms, terms[i]);
#             end
#         end
#     end
#     terms = newterms;
    
#     # If coefficients need to be computed, do so
#     # # First remove duplicates
#     unique_coef = [];
#     unique_coef_name = [];
#     unique_coef_ind = [];
#     unique_coef_deriv = [];
#     for i=1:length(needed_coef)
#         already = false;
#         for j=1:length(unique_coef)
#             if unique_coef[j] === needed_coef[i] && unique_coef_ind[j] == needed_coef_ind[i] && unique_coef_deriv[j] == needed_coef_deriv[i] && unique_coef_name[j] == needed_coef_name[i]
#                 already = true;
#             end
#         end
#         if !already
#             #println("coef: "*string(needed_coef[i])*" , ind: "*string(needed_coef_ind[i])*" , deriv: "*string(needed_coef_deriv[i]));
#             push!(unique_coef, needed_coef[i]);
#             push!(unique_coef_name, needed_coef_name[i]);
#             push!(unique_coef_ind, needed_coef_ind[i]);
#             push!(unique_coef_deriv, needed_coef_deriv[i]);
#         end
#     end
#     needed_coef = unique_coef;
#     needed_coef_name = unique_coef_name;
#     needed_coef_ind = unique_coef_ind;
#     needed_coef_deriv = unique_coef_deriv;
    
#     # For constant coefficients, this generates something like:
#     ######################################
#     # coef_n_i = a.value[i];
#     ######################################
    
#     # For function coefficients, this generates something like:
#     ######################################
#     # coef_n_i = zeros(refel.Np);
#     # for coefi = 1:refel.Np
#     #     coef_n_i[coefi] = a.value[i].func(nodex[1,coefi], nodex[2,coefi],nodex[3,coefi],time);
#     # end
#     ######################################
    
#     # For known variables, this generates something like:
#     ######################################
#     # coef_u_1 = copy((Femshop.variables[1]).values[loc2glb])
#     ######################################
    
#     if length(needed_coef) > 0
#         if vors == "volume"
#             cloop = :(for coefi=1:refel.Np end);
#         else
#             cloop = :(for coefi=1:refel.Nfp[frefelind[1]] end);
#         end
        
#         cloopin = Expr(:block);
        
#         coef_deriv_type = zeros(Int, length(needed_coef)); # 0:constant coef, 1:function, 2:variable values
#         need_integration = [];
#         for i=1:length(needed_coef)
#             if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
#                 tmps = needed_coef_name[i];
#                 tmpc = Symbol(tmps);
#                 dgside = 0;
#                 if vors == "surface"
#                     l2gsym = :(els[1])
#                     nodesym = :(facex)
#                 else
#                     l2gsym = :el
#                     nodesym = :nodex
#                 end
#                 if occursin("DGSIDE1", tmps)
#                     dgside = 1;
#                     l2gsym = :(els[1])
#                     nodesym = :(facex)
#                 elseif occursin("DGSIDE2", tmps)
#                     dgside = 2;
#                     l2gsym = :(els[2])
#                     nodesym = :(facex)
#                 end
#                 cargs = [:($nodesym[coefi]); 0; 0; :time];
#                 if config.dimension == 2
#                     cargs = [:($nodesym[1,coefi]); :($nodesym[2,coefi]); 0; :time];
#                 elseif config.dimension == 3
#                     cargs = [:($nodesym[1,coefi]); :($nodesym[2,coefi]); :($nodesym[3,coefi]); :time];
#                 end
                
#                 (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
#                 if ctype == 1
#                     # constant coefficient -> coef_n = cval
#                     tmpn = cval;
#                     push!(code.args, Expr(:(=), tmpc, tmpn));
#                 elseif ctype == 2
#                     # genfunction coefficients -> coef_n_i = coef.value[i].func(cargs)
#                     tmpv = :(a[coefi]);
#                     tmpv.args[1] = tmpc;
#                     tmpn = :(Femshop.genfunctions[$cval]); # Femshop.genfunctions[cval]
#                     tmpb = :(a.func());
#                     tmpb.args[1].args[1]= tmpn;
#                     append!(tmpb.args, cargs);
                    
#                     if vors == "volume"
#                         push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
#                     else
#                         push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Nfp[frefelind[1]])))); # allocate coef_n
#                     end
#                     push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
#                     coef_deriv_type[i] = 1;
#                     push!(need_integration, tmpc);
                    
#                 elseif ctype == 3
#                     # variable values -> coef_n = variable.values
#                     if variables[cval].type == SCALAR
#                         if occursin("D1", tmps) && vors == "surface"
#                             tmpb = :(Femshop.variables[$cval].values[els[2]] - Femshop.variables[$cval].values[els[1]]); 
#                         else
#                             tmpb = :(Femshop.variables[$cval].values[$l2gsym]); 
#                         end
#                     else
#                         compo = needed_coef_ind[i];
#                         tmpb = :(Femshop.variables[$cval].values[$compo, $l2gsym]);
#                     end
                    
#                     push!(code.args, Expr(:(=), tmpc, tmpb));
#                     coef_deriv_type[i] = 2;
#                 end
                
#             end# number?
#         end# coef loop
        
#         # Write the loop that computes coefficient values
#         if length(cloopin.args) > 0
#             cloop.args[2] = cloopin;
#             push!(code.args, cloop); # add loop to code
#         end
        
#         # Apply derivatives after the initializing loop
#         # Will look like: coef_i_j = D1x * coef_i_j
#         made_dxdydz = false;
#         for i=1:length(needed_coef_deriv)
#             if length(needed_coef_deriv[i][2]) > 0 && !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
#                 tmps = needed_coef_name[i];
#                 tmpc = Symbol(tmps);
                
#                 dmatname = make_deriv_matrix_name_fv(needed_coef_deriv[i][2]);
#                 dmat = Symbol(dmatname);
#                 if coef_deriv_type[i] == 0 # derivatives of constant coefficients should be zero
#                     tmpb= 0;
#                     push!(code.args, Expr(:(=), tmpc, tmpb));
#                 elseif coef_deriv_type[i] == 1 # derivative of function evaluated at nodes
#                     tmpb= :($dmat * $tmpc);
#                     push!(code.args, Expr(:(=), tmpc, tmpb));
#                 elseif coef_deriv_type[i] == 2 # derivative of variable values defined on cells
#                     # make dxdydz if needed
#                     if !made_dxdydz
#                         push!(code.args, :(dxyz = norm(cellx[2] - cellx[1]) .* normal));
#                     end
#                     if occursin("D1x", dmatname)
#                         push!(code.args, :(deriv_index = 1));
#                     elseif occursin("D1y", dmatname)
#                         push!(code.args, :(deriv_index = 2));
#                     elseif occursin("D1z", dmatname)
#                         push!(code.args, :(deriv_index = 3));
#                     end
#                     push!(code.args, :($tmpc = (els[1] != els[2] && abs(normal[deriv_index]) > 1e-10) ? $tmpc ./ dxyz[deriv_index] : 0));
#                 end
                
#             else
#                 # derivatives of constant coefficients should be zero
#             end
#         end
        
#         # Integrate the coefficient over the cell
#         # coef = Q * coef
#         if length(need_integration) > 0
#             # Quadrature matrix
#             push!(code.args, :(surf_Q = (refel.surf_wg[frefelind[1]] .* face_detJ)' * refel.surf_Q[frefelind[1]][:,refel.face2local[frefelind[1]]]));
#             for ici=1:length(need_integration)
#                 tmpc = need_integration[ici];
#                 push!(code.args, :($tmpc = (surf_Q * $tmpc) / area));
#             end
#         end
        
#     end# needed_coef loop
    
#     #=
#     # finally add the code expression
#     # For multiple dofs per node it will be like:
#     LHS:    [row quadrature vector(Np*dofs length)]
#             [one row per dof                      ]
            
#     RHS:    [a1] cell averages for each dof in a vector
#             [a2]
            
#     =#
    
#     # Allocate if needed
#     dofsper = 0;
#     if typeof(var) <: Array
#         for vi=1:length(var)
#             dofsper = dofsper + length(var[vi].symvar.vals); # The number of components for this variable
#         end
#     else
#         dofsper = length(var.symvar.vals);
#     end
    
#     if dofsper > 1
#         if lorr == RHS
#             push!(code.args, Expr(:(=), :cell_average, :(zeros($dofsper)))); # allocate result vector
#         else
#             push!(code.args, Expr(:(=), :cell_Q_matrix, :(zeros($dofsper, 2*$dofsper)))); # allocate quadrature matrix
#         end
#     end
    
#     result = nothing; # Will hold the returned expression
#     if typeof(var) <: Array # multivar
#         total_terms = 0;
#         for vi=1:length(var)
#             total_terms += length(terms[vi]);
#             # Add terms
#             # Each component/dof should have one expression so that submatrix is only modified once.
#             if lorr == LHS
#                 comps = dofsper;
#                 submatrices = Array{Any,2}(undef, comps, comps);
#                 for smi=1:length(submatrices)
#                     submatrices[smi] = nothing;
#                 end
#                 for i=1:length(terms[vi])
#                     ti = offset_ind[vi];
#                     tj = var_ind[vi][i][1];
                    
#                     if submatrices[ti, tj] === nothing
#                         submatrices[ti, tj] = terms[vi][i];
#                     else
#                         addexpr = :(a+b);
#                         addexpr.args[2] = submatrices[ti, tj];
#                         addexpr.args[3] = terms[vi][i];
#                         submatrices[ti, tj] = addexpr;
#                     end
#                 end
                
#                 for cj=1:comps
#                     for ci=1:comps
#                         if !(submatrices[ci, cj] === nothing)
#                             ti = ci;
#                             tj = cj-1;
#                             stj = :($tj*refel.Np + 1);
#                             enj = :(($tj + 1)*refel.Np);
                            
#                             push!(code.args, Expr(:(+=), :(cell_Q_matrix[$ti, $stj:$enj]), submatrices[ci, cj]));
#                         end
#                     end
#                 end
                
#                 result = :cell_Q_matrix;
                
#             else #RHS
#                 comps = dofsper;
#                 submatrices = Array{Any,1}(undef, comps);
#                 for smi=1:length(submatrices)
#                     submatrices[smi] = nothing;
#                 end
#                 for i=1:length(terms[vi])
#                     ti = 1 + offset_ind[vi];
                    
#                     if submatrices[ti] === nothing
#                         submatrices[ti] = terms[vi][i];
#                     else
#                         addexpr = :(a+b);
#                         addexpr.args[2] = submatrices[ti];
#                         addexpr.args[3] = terms[vi][i];
#                         submatrices[ti] = addexpr;
#                     end
#                 end
                
#                 for ci=1:comps
#                     if !(submatrices[ci] === nothing)
                        
#                         push!(code.args, Expr(:(+=), :(cell_average[$ci]), submatrices[ci]));
#                     end
#                 end
                
#                 result = :cell_average;
                
#             end
#         end
        
#         if total_terms == 0
#             # there were no terms. Just return 0
#             if lorr == LHS
#                 return :(return nothing);
#             else
#                 return :(return zeros($dofsper));
#             end
#         end
        
#     else
#         if length(terms) > 1
#             if var.type == SCALAR # Only one component
#                 tmp = :(a.+b);
#                 tmp.args = [:.+];
#                 for i=1:length(terms)
#                     push!(tmp.args, terms[i]);
#                 end
#                 result = tmp;
#             else # More than one component
#                 # Add terms
#                 # Each component/dof should have one expression so that submatrix is only modified once.
#                 if lorr == LHS
#                     comps = length(var.symvar.vals);
#                     submatrices = Array{Any,2}(undef, comps, comps);
#                     for smi=1:length(submatrices)
#                         submatrices[smi] = nothing;
#                     end
#                     for i=1:length(terms)
#                         ti = var_ind[i][1];
#                         tj = var_ind[i][1];
                        
#                         if submatrices[ti, tj] === nothing
#                             submatrices[ti, tj] = terms[i];
#                         else
#                             addexpr = :(a.+b);
#                             addexpr.args[2] = submatrices[ti, tj];
#                             addexpr.args[3] = terms[i];
#                             submatrices[ti, tj] = addexpr;
#                         end
#                     end
                    
#                     for cj=1:comps
#                         for ci=1:comps
#                             if !(submatrices[ci, cj] === nothing)
#                                 ti = ci;
#                                 tj = cj-1;
#                                 stj = :($tj*refel.Np + 1);
#                                 enj = :(($tj + 1)*refel.Np);
                                
#                                 push!(code.args, Expr(:(+=), :(cell_Q_matrix[ti, $stj:$enj]), submatrices[ci, cj]));
#                             end
#                         end
#                     end
                    
#                     result = :cell_Q_matrix;
                    
#                 else #RHS
#                     comps = length(var.symvar.vals);
#                     submatrices = Array{Any,1}(undef, comps);
#                     for smi=1:length(submatrices)
#                         submatrices[smi] = nothing;
#                     end
#                     for i=1:length(terms)
#                         ti = var_ind[i][1];
                        
#                         if submatrices[ti] === nothing
#                             submatrices[ti] = terms[i];
#                         else
#                             addexpr = :(a.+b);
#                             addexpr.args[2] = submatrices[ti];
#                             addexpr.args[3] = terms[i];
#                             submatrices[ti] = addexpr;
#                         end
#                     end
                    
#                     for ci=1:comps
#                         if !(submatrices[ci] === nothing)
                            
#                             push!(code.args, Expr(:(+=), :(cell_average[$ci]), submatrices[ci]));
#                         end
#                     end
                    
#                     result = :cell_average;
                    
#                 end
                
#             end
#         elseif length(terms) == 1# one term (one variable)
#             result = terms[1];
            
#         else # there were no terms. Just return 0
#             if lorr == LHS
#                 return :(return nothing);
#             else
#                 return :(return 0);
#             end
#         end
#     end
    
#     # At this point everything is packed into result
#     push!(code.args, Expr(:return, result));
#     return code;
# end

# #=
# Changes the symbolic layer source term into a code layer term.
# An example term with known variable u, coefficient c
# Q is a quadrature vector like: (weights*detJ)' * refel.Q

# 2*_c_1*_u_1     ->  2 * (Q * (coef_0_1 .* coef_u_1))

# 2*D2__c_1*D1_D1__u_1    ->  2 * (Q * ((D1y * coef_0_1) .* (D2x * coef_u_1)))

# Note: These result in a single value for the cell.

# Also records the needed derivative matrices like: J.rx * refel.Ddr
# D1 -> D1x
# D2 -> D1y
# D1_D1 -> D2x
# etc.

# Returns:
# - expression for term
# - needed coefficients
# - coefficient indices(vector components)
# - assigned names for coefficients (_c_1 -> coef_n_1 when c is the nth coefficient in coefficients array)
# - needed derivative matrices

# =#
# function process_source_term_fv_julia(sterm, var, lorr, offset_ind=0)
#     if typeof(sterm) == Symbol
#         term = sterm;
#     else
#         term = copy(sterm);
#     end
    
#     need_derivative = false;
#     need_derivative_for_coefficient = false;
#     needed_coef = [];
#     needed_coef_name = [];
#     needed_coef_ind = [];
#     needed_coef_deriv = [];
#     needed_derivatives = [];
    
#     var_part = nothing;
#     var_derivs = [];
#     coef_part = nothing;
#     const_part = 1;
#     var_component = 0;
    
#     # extract each of the factors.
#     factors = separate_factors(term);
    
#     # strip off all negatives, combine and reattach at the end
#     neg = false;
#     for i=1:length(factors)
#         if typeof(factors[i]) == Expr && factors[i].args[1] === :- && length(factors[i].args) == 2
#             neg = !neg;
#             factors[i] = factors[i].args[2];
#         end
#     end
    
#     # Separate factors into variable/coefficient parts
#     coef_facs = [];
#     coef_inds = [];
#     coef_derivs = [];
#     coef_mods = [];
#     coef_expr_facs = [];
#     for i=1:length(factors)
#         if typeof(factors[i]) <: Number
#             push!(coef_facs, factors[i]);
#             push!(coef_inds, -1);
#             push!(coef_derivs, [nothing, "", ""]);
            
#         elseif typeof(factors[i]) == Expr && factors[i].head === :call
#             # These should both be purely coefficient/known expressions. 
#             if factors[i].args[1] === :./
#                 # The second arg should be 1, the third should not contain an unknown
#                 # The denominator expression needs to be processed completely
#                 (piece, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[3]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[3] = piece;
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
                
#             elseif factors[i].args[1] === :.^
#                 # The second arg is the thing raised
#                 (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[2]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[2] = piece1;
#                 # Do the same for the power just in case
#                 (piece2, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[3]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[3] = piece2;
                
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
#             elseif factors[i].args[1] === :sqrt
#                 factors[i].args[1] = :.^
#                 # The second arg is the thing sqrted
#                 (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[2]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[2] = piece1;
#                 # add a 1/2 power argument
#                 push!(factors[i].args, 1/2);
                
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
#             end
            
#         else
#             (index, v, mods) = extract_symbols(factors[i]);
            
#             if is_unknown_var(v, var) && lorr == LHS # If rhs, treat as a coefficient
#                 if !(var_part === nothing)
#                     # Two unknowns multiplied in this term. Nonlinear LHS. abort.
#                     printerr("Nonlinear term found on LHS. Code layer incomplete.");
#                     return (-1, -1, -1, -1, -1, -1);
#                 end
#                 var_component = index;
#                 #offset component for multivar
#                 mvar_var = var;
#                 if typeof(var) <:Array
#                     for vi=1:length(var)
#                         if v === var[vi].symbol
#                             var_component = var_component .+ offset_ind[vi];
#                             mvar_var = var[vi];
#                         end
#                     end
#                 end
#                 if length(mods) > 0
#                     need_derivative = true;
#                     (derivs, others) = get_derivative_mods_fv(mods);
#                     dmatname = make_deriv_matrix_name_fv(derivs);
#                     push!(needed_derivatives, derivs);
#                     var_part = dmatname; # LHS part is just the matrix
#                 else
#                     # no derivative mods
#                     var_part = :(Iamnotsure);
#                 end
#             else # coefficients or RHS vars
#                 if length(index) == 1 # This should always be the case
#                     ind = index[1];
#                 end
#                 # Check for derivative mods
#                 if typeof(v) == Symbol && !(v ===:dt)
#                     if length(mods) > 0
#                         need_derivative = true;
#                         need_derivative_for_coefficient = true;
#                         (derivs, others) = get_derivative_mods_fv(mods);
#                         dmatname = make_deriv_matrix_name_fv(derivs);
#                         push!(needed_derivatives, derivs);
                        
#                         push!(needed_coef_deriv, [v, derivs]);
                        
#                         push!(coef_mods, others);
                        
#                         push!(needed_coef_name, make_coef_name_fv(v, others, derivs, ind));
                        
#                     else
#                         push!(needed_coef_deriv, [v, []]);
#                         push!(coef_mods, []);
#                         push!(needed_coef_name, make_coef_name_fv(v, [], [], ind));
#                     end
#                     push!(needed_coef, v);
#                     push!(needed_coef_ind, ind);
                    
#                     push!(coef_derivs, needed_coef_deriv[end]);
                    
                    
#                 else
#                     push!(coef_derivs, [nothing, []]);
#                 end
                
#                 push!(coef_facs, v);
#                 push!(coef_inds, ind);
#             end
#         end
        
#     end # factors loop
    
#     # If there's no var part, need to do this
#     if var_part === nothing
#         var_part = :(Iamnotsure);
#     end
    
#     # build coefficient parts
#     name_ind = 0;
#     multbydt = false;
#     if length(coef_facs) > 0
#         for j=1:length(coef_facs)
#             if typeof(coef_facs[j]) <: Number
#                 # Numbers are put in the constant part
#                 const_part = const_part * coef_facs[j];
#             elseif coef_facs[j] ===:dt
#                 multbydt = true;
#             else
#                 tmp = coef_facs[j];
#                 #println("coef_facs: "*string(tmp)*" : "*string(typeof(tmp)));
#                 if typeof(tmp) == Symbol && !(tmp ===:dt)
#                     name_ind += 1;
#                     tmp = Symbol(needed_coef_name[name_ind]);
                    
#                 end
#                 if !(coef_part === nothing)
#                     coef_part = :($coef_part .* $tmp);
#                 else
#                     coef_part = tmp;
#                 end
#             end
#         end
#     end
    
#     if length(coef_expr_facs) > 0
#         for j=1:length(coef_expr_facs)
#             tmp = coef_expr_facs[j];
#             if !(coef_part === nothing)
#                 coef_part = :($coef_part .* $tmp);
#             else
#                 coef_part = tmp;
#             end
#         end
#     end
    
#     # Now the const_part, coef_part and var_part are complete.
#     # Build an expression for the term.
#     # 2*D2__c_1*D1__u_1    ->  (Q .* 2 .* (D1y * coef_0_1)) * D1x               on LHS
#     #                      ->  2 * (Q * ((D1y * coef_0_1) .* (???????????))) on RHS
#     qc_part = :Q;
#     if multbydt
#         const_part = :(($const_part * dt));
#     end
    
#     if lorr == LHS
#         if !(const_part == 1)
#             qc_part = :($qc_part .* $const_part);
#         end
#         if !(coef_part === nothing)
#             qc_part = :($qc_part .* coef_part);
#         end
#         if qc_part === :Q && var_part == :Iamnotsure
#             term = 1;
#         else
#             term = :(($qc_part) * $var_part);
#         end
#     else #RHS
#         if const_part == 1
#             term = :($coef_part);
#         elseif !(coef_part === nothing)
#             term = :($const_part * ($coef_part));
#         else#only constant
#             term = :($const_part);
#         end
#     end
    
#     # Apply the negative if needed
#     if neg
#         negex = :(-a);
#         negex.args[2] = term;
#         term = negex;
#     end
    
#     return (term, needed_coef, needed_coef_ind, needed_coef_name, needed_coef_deriv, var_component);
# end

# #=
# Changes the symbolic layer flux term into a code layer term.
# An example term with known variable u, coefficient c
# Q is a quadrature vector like: (weights*detJ)' * refel.surf_Q[face]

# 2*_c_1*DGSIDE1__u_1     ->  2 * (Q * (coef_0_1 .* coef_DGSIDE1_u_1))

# 2*D2__c_1*DGSIDE1_D1_D1__u_1    ->  2 * (Q * ((D1y * coef_0_1) .* (D2x * coef_DGSIDE1_u_1)))

# Note: These result in a single value for the cell.

# Also records the needed derivative matrices like: J.rx * refel.Ddr
# D1 -> D1x
# D2 -> D1y
# D1_D1 -> D2x
# etc.

# Returns:
# - expression for term
# - needed coefficients
# - coefficient indices(vector components)
# - assigned names for coefficients (_c_1 -> coef_n_1 when c is the nth coefficient in coefficients array)
# - needed derivative matrices

# =#
# function process_flux_term_fv_julia(sterm, var, lorr, offset_ind=0)
#     if typeof(sterm) == Symbol
#         term = sterm;
#     else
#         term = copy(sterm);
#     end
    
#     need_derivative = false;
#     need_derivative_for_coefficient = false;
#     needed_coef = [];
#     needed_coef_name = [];
#     needed_coef_ind = [];
#     needed_coef_deriv = [];
#     needed_derivatives = [];
    
#     coef_part = nothing;
#     var_part = nothing;
#     weight_part = :face_wdetj;
#     const_part = 1;
#     var_component = 0;
    
#     # extract each of the factors.
#     factors = separate_factors(term);
    
#     # strip off all negatives, combine and reattach at the end
#     neg = false;
#     for i=1:length(factors)
#         if typeof(factors[i]) == Expr && factors[i].args[1] === :- && length(factors[i].args) == 2
#             neg = !neg;
#             factors[i] = factors[i].args[2];
#         end
#     end
    
#     # Separate factors into var/coefficient parts
#     coef_facs = [];
#     coef_inds = [];
#     coef_names = [];
#     coef_derivs = [];
#     coef_mods = [];
#     coef_expr_facs = [];
#     for i=1:length(factors)
#         if typeof(factors[i]) <: Number
#             push!(coef_facs, factors[i]);
#             push!(coef_inds, -1);
#             push!(coef_names, "");
#             push!(coef_derivs, [nothing, "", ""]);
            
#         elseif typeof(factors[i]) == Expr && factors[i].head === :call
#             # These should both be purely coefficient/known expressions. 
#             if factors[i].args[1] === :./
#                 # The second arg should be 1, the third should not contain an unknown
#                 # The denominator expression needs to be processed completely
#                 (piece, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[3]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[3] = piece;
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
                
#             elseif factors[i].args[1] === :.^
#                 # The second arg is the thing raised
#                 (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[2]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[2] = piece1;
#                 # Do the same for the power just in case
#                 (piece2, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[3]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[3] = piece2;
                
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
#             elseif factors[i].args[1] === :sqrt
#                 factors[i].args[1] = :.^
#                 # The second arg is the thing sqrted
#                 (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[2]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[2] = piece1;
#                 # add a 1/2 power argument
#                 push!(factors[i].args, 1/2);
                
#                 push!(coef_expr_facs, factors[i]);
#                 #push!(coef_inds, 0);
#             elseif factors[i].args[1] === :abs
#                 # The second arg is the thing in abs
#                 (piece1, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(factors[i].args[2]);
#                 need_derivative = need_derivative || nd;
#                 append!(needed_coef, nc);
#                 append!(needed_coef_name, ncm);
#                 append!(needed_coef_ind, nci);
#                 append!(needed_coef_deriv, ncd);
#                 factors[i].args[2] = piece1;
                
#                 push!(coef_expr_facs, factors[i]);
#             end
            
#         else
#             (index, v, mods) = extract_symbols(factors[i]);
            
#             if is_unknown_var(v, var) && lorr == LHS # If rhs, treat as a coefficient
#                 if !(var_part === nothing)
#                     # Two unknowns multiplied in this term. Nonlinear LHS. abort.
#                     printerr("Nonlinear term found on LHS. Code layer incomplete.");
#                     return (-1, -1, -1, -1, -1, -1);
#                 end
#                 var_component = index;
#                 #offset component for multivar
#                 mvar_var = var;
#                 if typeof(var) <:Array
#                     for vi=1:length(var)
#                         if v === var[vi].symbol
#                             var_component = var_component .+ offset_ind[vi];
#                             mvar_var = var[vi];
#                         end
#                     end
#                 end
#                 if length(mods) > 0
#                     need_derivative = true;
#                     (derivs, others) = get_derivative_mods_fv(mods);
#                     dmatname = make_deriv_matrix_name_fv(derivs);
#                     push!(needed_derivatives, derivs);
#                     var_part = dmatname; # LHS part is just the matrix
#                 else
#                     # no derivative mods
#                     var_part = :(Iamnotsure);
#                 end
#             else # coefficients or RHS vars
#                 if length(index) == 1 # This should always be the case
#                     ind = index[1];
#                 end
#                 # Check for derivative mods
#                 if typeof(v) == Symbol && !(v ===:dt)
#                     if length(mods) > 0
#                         need_derivative = true;
#                         need_derivative_for_coefficient = true;
#                         (derivs, others) = get_derivative_mods_fv(mods);
#                         dmatname = make_deriv_matrix_name_fv(derivs);
#                         push!(needed_derivatives, derivs);
#                         push!(needed_coef_deriv, [v, derivs]);
#                         push!(coef_mods, others);
#                         push!(needed_coef_name, make_coef_name_fv(v, others, derivs, ind));
#                         push!(coef_names, make_coef_name_fv(v, others, derivs, ind));
                        
#                     else
#                         push!(needed_coef_deriv, [v, []]);
#                         push!(coef_mods, []);
#                         push!(needed_coef_name, make_coef_name_fv(v, [], [], ind));
#                         push!(coef_names, make_coef_name_fv(v, [], [], ind));
#                     end
#                     push!(needed_coef, v);
#                     push!(needed_coef_ind, ind);
                    
#                     push!(coef_derivs, needed_coef_deriv[end]);
                    
                    
#                 else
#                     push!(coef_derivs, [nothing, []]);
#                     push!(coef_names, "");
#                 end
                
#                 push!(coef_facs, v);
#                 push!(coef_inds, ind);
#             end
#         end
        
#     end # factors loop
    
#     # build coefficient parts
#     name_ind = 0;
#     multbydt = false;
#     if length(coef_facs) > 0
#         for j=1:length(coef_facs)
#             tmp = coef_facs[j];
#             #println("coef_facs: "*string(tmp)*" : "*string(typeof(tmp)));
#             if typeof(coef_facs[j]) <: Number
#                 # Numbers are put in the constant part
#                 const_part = const_part * coef_facs[j];
#             elseif tmp ===:dt
#                 multbydt = true;
#             elseif typeof(tmp) == Symbol && !(tmp ===:dt)
#                 name_ind += 1;
#                 tmp = Symbol(coef_names[j]);
                
#                 if occursin("FACENORMAL1", coef_names[j])
#                     ind = coef_inds[j];
#                     tmp = :(normal[$ind]);
#                 elseif occursin("FACENORMAL2", coef_names[j])
#                     ind = coef_inds[j];
#                     tmp = :(-normal[$ind]);
#                 end
                
#                 if !(coef_part === nothing)
#                     coef_part = :($coef_part .* $tmp);
#                 else
#                     coef_part = tmp;
#                 end
#             end
#         end
#     end
    
#     # These are factors coming from subexpressions
#     if length(coef_expr_facs) > 0
#         for j=1:length(coef_expr_facs)
#             tmpc = coef_expr_facs[j];
#             if !(coef_part === nothing)
#                 cp = coef_part;
#                 coef_part = :($cp .* $tmpc);
#             else
#                 coef_part = tmpc;
#             end
#         end
#     end
    
#     # Now the const_part, coef_part and var_part are complete.
#     # Build an expression for the term.
#     # 2*D2__c_1*D1_D1__u_1    ->  (Q .* 2 .* (D1y * coef_0_1)) * D2x               on LHS
#     #                         ->  2 * (Q * ((D1y * coef_0_1) .* (D2x * coef_u_1))) on RHS
#     qc_part = :(surf_Q);
#     if multbydt
#         const_part = :(($const_part * dt));
#     end
    
#     if lorr == LHS
#         if !(const_part == 1)
#             qc_part = :($qc_part .* $const_part);
#         end
#         if !(coef_part === nothing)
#             qc_part = :($qc_part .* coef_part);
#         end
#         if qc_part === :(surf_Q) && var_part == :Iamnotsure
#             term = 1;
#         else
#             term = :(($qc_part) * $var_part);
#         end
        
#     else #RHS
#         if const_part == 1
#             term = :($coef_part);
#         elseif !(coef_part === nothing)
#             term = :($const_part * ($coef_part));
#         else#only constant
#             term = :($const_part);
#         end
#     end
    
#     # Apply the negative if needed
#     if neg
#         negex = :(-a);
#         negex.args[2] = term;
#         term = negex;
#     end
    
#     return (term, needed_coef, needed_coef_ind, needed_coef_name, needed_coef_deriv, var_component);
# end

# # Special processing for sub expressions of known things.(denominators, sqrt, etc.)
# # It should only contain knowns/coeficients, so just replace symbols (f -> coef_0_1)
# # And return the needed coefficient info
# function process_known_expr_fv_julia(ex)
#     need_derivative = false;
#     needed_coef = [];
#     needed_coef_name = [];
#     needed_coef_ind = [];
#     needed_coef_deriv = [];
    
#     # Work recursively through the expression
#     if typeof(ex) <: Number
#         return (ex, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
        
#     elseif typeof(ex) == Symbol
#         # turn arithmetic ops into dotted versions
#         if ex === :+ || ex === :.+
#             return (:.+ , false, [], [], [], []);
#         elseif ex === :- || ex === :.-
#             return (:.- , false, [], [], [], []);
#         elseif ex === :* || ex === :.*
#             return (:.* , false, [], [], [], []);
#         elseif ex === :/ || ex === :./
#             return (:./ , false, [], [], [], []);
#         elseif ex === :^ || ex === :.^
#             return (:.^ , false, [], [], [], []);
#         elseif ex === :abs
#             return (:abs , false, [], [], [], []);
#         end
        
#         (index, v, mods) = extract_symbols(ex);
#         if length(index) == 1
#             ind = index[1];
#         end
        
#         if !(v ===:dt)
#             if v === :FACENORMAL
#                 tmp = :(normal[$ind]);
#             elseif v === :FACENORMAL1
#                 tmp = :(normal[$ind]);
#             elseif v === :FACENORMAL2
#                 tmp = :(-normal[$ind]);
#             else
#                 # Check for derivative mods
#                 (derivs, others) = get_derivative_mods_fv(mods);
#                 push!(needed_coef_deriv, [v, derivs]);
#                 push!(needed_coef, v);
#                 push!(needed_coef_ind, ind);
                
#                 tmps = make_coef_name_fv(v, others, derivs, ind);
#                 push!(needed_coef_name, tmps);
#                 tmp = Symbol(tmps); # The symbol to return
#             end
#         else
#             tmp = ex;
            
#         end
        
#         return (tmp, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
        
#     elseif typeof(ex) == Expr
#         newex = copy(ex);
#         for i=1:length(ex.args)
#             (piece, nd, nc, ncm, nci, ncd) = process_known_expr_fv_julia(ex.args[i]);
#             newex.args[i] = piece;
#             need_derivative = need_derivative || nd;
#             append!(needed_coef, nc);
#             append!(needed_coef_name, ncm);
#             append!(needed_coef_ind, nci);
#             append!(needed_coef_deriv, ncd);
#         end
        
#         return (newex, need_derivative, needed_coef, needed_coef_name, needed_coef_ind, needed_coef_deriv);
#     end
    
# end

# ###############################################################################################################
# # utils
# ###############################################################################################################

# # # Separates expressions into terms.
# # # Assumes the expressions are expanded as they should be coming from the parser.
# # # 2*a + b*c*5 - w*2 -> [2*a, b*c*5, w*2]
# # function separate_terms(ex)
# #     terms = [ex];
# #     if typeof(ex) == Expr && ex.head === :call
# #         if ex.args[1] === :+ || (ex.args[1] === :- && length(ex.args) > 2)
# #             terms = [];
# #             for i=2:length(ex.args)
# #                 append!(terms, separate_terms(ex.args[i]));
# #             end
# #         end
# #     end
    
# #     return terms;
# # end

# # # Parses symengine terms into julia expressions
# # function terms_to_expr(symex)
# #     terms = [];
# #     sz = size(symex);
# #     if length(sz) == 1 # scalar or vector
# #         for i=1:length(symex)
# #             for ti=1:length(symex[i])
# #                 push!(terms, Meta.parse(string(symex[i][ti])));
# #             end
# #         end
# #     elseif length(sz) == 2 # matrix
# #         #TODO
# #         printerr("sorry. still need to implement code layer for tensors.")
# #     end
    
# #     # for i=1:length(terms)
# #     #     convert_sqrt_to_power(terms[i]);
# #     # end
    
# #     return terms;
# # end

# # # Separates terms into factors.
# # # Assumes the term is only multiplied symbols or numbers
# # # 2*a*thing -> [2, a, thing]
# # function separate_factors(ex)
# #     factors::Array{Any,1} = [ex];
# #     if typeof(ex) == Expr && ex.head === :call
# #         if ex.args[1] === :* || ex.args[1] === :.*
# #             factors = [];
# #             for i=2:length(ex.args)
# #                 append!(factors, separate_factors(ex.args[i]));
# #             end
            
# #         elseif ex.args[1] === :^ || ex.args[1] === :.^
# #             factors = [];
# #             power = ex.args[3];
# #             if power == 1
# #                 factors = [ex.args[2]]; #a^1 = a
# #             else
# #                 ex.args[1] = :.^ ;
# #                 factors = [ex]; # the power is handled later
# #             end
            
# #         elseif ex.args[1] === :/ || ex.args[1] === :./
# #             factors = [];
# #             append!(factors, separate_factors(ex.args[2])); # a/b = a * 1/b
# #             divex = :(1 ./ a);
# #             divex.args[3] = ex.args[3];
# #             push!(factors, divex);
            
# #         elseif ex.args[1] === :sqrt
# #             ex.args[1] = :.^ ;
# #             push!(ex.args, 0.5);
# #             factors = [ex];
            
# #         elseif ex.args[1] === :- && length(ex.args) == 2
# #             # strip off negetive and place on first factor
# #             subex = ex.args[2];
# #             factors = separate_factors(subex);
# #             negex = :(-a);
# #             negex.args[2] = factors[1];
# #             factors[1] = negex;
# #         end
# #     end
    
# #     return factors;
# # end

# # # Checks if the coefficient has constant value.
# # # If so, also returns the value.
# # function is_constant_coef(c)
# #     isit = false;
# #     val = 0;
# #     for i=1:length(coefficients)
# #         if c === coefficients[i].symbol
# #             isit = (typeof(coefficients[i].value[1]) <: Number);
# #             if isit
# #                 val = coefficients[i].value[1];
# #             else
# #                 val = coefficients[i].value[1].name;
# #             end
# #         end
# #     end
    
# #     return (isit, val);
# # end

# # # Checks the type of coefficient: constant, genfunction, or variable
# # # Returns: (type, val)
# # # constant: type=1, val=number
# # # genfunction: type=2, val= index in genfunctions array
# # # variable: type=3, val=index in variables array
# # function get_coef_val(c, comp)
# #     type = 0;
# #     val = 0;
# #     for i=1:length(coefficients)
# #         if c === coefficients[i].symbol
# #             isit = (typeof(coefficients[i].value[comp]) <: Number);
# #             if isit
# #                 type = 1;
# #                 val = coefficients[i].value[comp];
# #             else
# #                 type = 2;
# #                 name = coefficients[i].value[comp].name;
# #                 for j=1:length(genfunctions)
# #                     if name == genfunctions[j].name
# #                         val = j;
# #                     end
# #                 end
# #             end
# #         end
# #     end
# #     if type == 0
# #         for i=1:length(variables)
# #             if c === variables[i].symbol
# #                 type = 3;
# #                 val = variables[i].index;
# #             end
# #         end
# #     end
    
# #     return (type, val);
# # end

# # function get_coef_index(c)
# #     ind = -1;
# #     for i=1:length(coefficients)
# #         if c === coefficients[i].symbol
# #             ind = coefficients[i].index;
# #         end
# #     end
    
# #     return ind;
# # end

# function make_coef_name_fv(c, othermods, derivs, cind)
#     ind = get_coef_index(c);
#     if ind >= 0
#         name = string(ind);
#     else
#         name = string(c);
#     end
    
#     if occursin("FACENORMAL1", name)
#         return "FACENORMAL1_"*string(cind);
#     elseif occursin("FACENORMAL2", name)
#         return "FACENORMAL2_"*string(cind);
#     end
    
#     name = make_deriv_matrix_name_fv(derivs) * name;
#     for i=length(othermods):-1:1
#         name = string(othermods[i]) * "_" * name;
#     end
#     name = "coef_"*name*"_"*string(cind);
#     return name;
# end

# # function is_unknown_var(v, vars)
# #     if typeof(vars) == Variable
# #         return v===vars.symbol;
# #     end
# #     for t in vars
# #         if t.symbol === v
# #             return true;
# #         end
# #     end
# #     return false;
# # end

# # # Extract meaning from the symbolic object name
# # # The format of the input symbol should look like this
# # #   MOD1_MOD2_..._var_n
# # # There could be any number of mods, _var_n is the symvar symbol (n is the vector index, or 1 for scalar)
# # # Returns ([n], var, [MOD1,MOD2,...]) all as strings
# # function extract_symbols(ex)
# #     str = string(ex);
# #     #println("extracting from: "*str);
# #     index = [];
# #     var = nothing;
# #     mods = [];
# #     l = lastindex(str);
# #     e = l; # end of variable name
# #     b = l; # beginning of variable name
    
# #     # dt is a special symbol that will be assigned a number value in the generated function.
# #     if str == "dt"
# #         return([0], ex, []);
# #     end
    
# #     for i=l:-1:0
# #         if e==l
# #             if str[i] == '_'
# #                 e = i-1;
# #             else
# #                 # These are the indices at the end. Parse one digit at a time.
# #                 try
# #                     index = [parse(Int, str[i]); index] # The indices on the variable
# #                 catch
# #                     return ([0],ex,[]);
# #                 end
                
# #             end
# #         elseif b==l && i>0
# #             if str[i] == '_'
# #                 b = i+1;
# #             end
            
# #         else
# #             # At this point we know b and e
# #             if var === nothing
# #                 var = Symbol(SubString(str, b, e));
# #                 b = b-2;
# #                 e = b;
# #             end
# #         end
# #     end
    
# #     # Change index from [a, b] to [a*d + b]
# #     newindex = index[end];
# #     for i=1:(length(index)-1)
# #         newindex = newindex + (index[i]-1)*config.dimension^(length(index)-i);
# #     end
# #     index = [newindex];
    
# #     # extract the modifiers like D1_ 
# #     if b>1
# #         e = b-1;
# #         for i=e:-1:1
# #             if str[i] == '_'
# #                 push!(mods, SubString(str, b, e));
                
# #                 e = i-1;
                
# #             elseif i == 1
# #                 push!(mods, SubString(str, 1, e));
# #             end
# #             b = i;
# #         end
# #     end
    
# #     #println("got: "*string(mods)*" "*string(var)*" "*string(index));
    
# #     return (index, var, mods);
# # end

# # Separate the derivative mods(Dn_) from the other mods and return an array 
# # of derivative indices and the remaining mods.
# function get_derivative_mods_fv(mods)
#     derivs = [];
#     others = [];
    
#     for i=1:length(mods)
#         if length(mods[i]) == 2 && mods[i][1] == 'D'
#             try
#                 index = parse(Int, mods[i][2])
#                 push!(derivs, index);
#             catch
#                 printerr("Unexpected modifier: "*mods[i]*" see get_derivative_mods_fv() in generate_code_layer")
#                 push!(others, mods[i]);
#             end
#         else
#             push!(others, mods[i]);
#         end
#     end
    
#     return (derivs, others);
# end

# # Makes a name for a derivative matrix
# # [1,2,3] -> D3zD2yD1x
# function make_deriv_matrix_name_fv(deriv)
#     dname = "";
#     dn = [0;0;0];
    
#     for i=1:length(deriv)
#         dn[deriv[i]] += 1;
#     end
    
#     if dn[1] > 0
#         dname = "D"*string(dn[1])*"x" * dname;
#     end
#     if dn[2] > 0
#         dname = "D"*string(dn[2])*"y" * dname;
#     end
#     if dn[3] > 0
#         dname = "D"*string(dn[3])*"z" * dname;
#     end
    
#     return dname;
# end