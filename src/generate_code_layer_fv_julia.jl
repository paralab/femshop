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
vol_J =     args[13]; # jacobian for both elements
time =      args[14]; # time for time dependent coefficients
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
                        code *= cname * " = " * cname * "[refel.face2local[frefelind[1]]]; # extract face values only.";
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
                cellside = 0; # 0 means no side flag
                for flagi=1:length(entities[i].flags)
                    if occursin("DGSIDE1", entities[i].flags[flagi]) || occursin("CELL1", entities[i].flags[flagi])
                        cellside = 1;
                    elseif occursin("DGSIDE2", entities[i].flags[flagi]) || occursin("CELL2", entities[i].flags[flagi])
                        cellside = 2;
                    end
                end
                if vors == "surface"
                    l2gsymbol = "els[1]"
                else
                    l2gsymbol = "el"
                end
                if cellside == 1
                    l2gsymbol = "els[1]"
                elseif cellside == 2
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
                        if cellside == 0
                            # No side was specified, so use the average
                            code *= cname * " = 0.5 * (Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", els[1]] + Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", els[2]]);\n";
                        else
                            code *= cname * " = Femshop.variables["*string(cval)*"].values["*string(entities[i].index)*", "*l2gsymbol*"];\n";
                        end
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
    # 0.1*D1__u_1*_FACENORMAL1_1    ->  ???               on LHS
    #                             ->  0.1 * (coef_D1xu_1 .* normal[1])     on RHS
    if lorr == LHS
        (test_part, var_part, coef_part, test_ind, var_ind) = separate_factors(term, var);
        # # LHS: ??
        # TODO
        #
        # Seriously, need to do
        
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
