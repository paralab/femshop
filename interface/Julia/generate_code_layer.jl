#=
Use the symbolic layer expressions to generate the FEM code
=#
function generate_code_layer(ex, rhsvar=nothing)
    if language == 0 || language == JULIA
        return generate_code_layer_julia(ex, rhsvar);
    elseif language == CPP
        return generate_code_layer_dendro(ex, rhsvar);
    elseif language == MATLAB
        # not ready
        printerr("Matlab code gen not yet supported")
    end
end

# Julia version returns an expression for the generated function for linear or bilinear term
function generate_code_layer_julia(ex, rhsvar=nothing)
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
    push!(code.args, :((detJ, J) = geometric_factors(refel, x)));
    
    need_derivative = false;
    needed_coef = [];
    
    # Separate terms
    terms = separate_terms(ex);
    
    # Process the terms turning them into the code layer
    cindex = 0;
    for i=1:length(terms)
        (codeterm, der, coe) = process_term_julia(terms[i], cindex, rhsvar);
        need_derivative = need_derivative || der;
        append!(needed_coef, coe);
        cindex = length(needed_coef);
        terms[i] = codeterm;
    end
    
    # If derivatives are needed, prepare the appropriate matrices
    if need_derivative
        if language == 0 # No external generation
            if config.dimension == 1
                push!(code.args, Expr(:(=), :Rxmatrix, :(diagm(J.rx))));
                push!(code.args, Expr(:(=), :Qrmatrix, :(refel.Qr)));
            elseif config.dimension == 2
                push!(code.args, Expr(:(=), :Rxmatrix, :([diagm(J.rx) diagm(J.sx); diagm(J.ry) diagm(J.sy)])));
                push!(code.args, Expr(:(=), :Qrmatrix, :([refel.Qr ; refel.Qs])));
            elseif config.dimension == 3
                push!(code.args, Expr(:(=), :Rxmatrix, :([diagm(J.rx) diagm(J.sx) diagm(J.tx); diagm(J.ry) diagm(J.sy) diagm(J.ty); diagm(J.rz) diagm(J.sz) diagm(J.tz)])));
                push!(code.args, Expr(:(=), :Qrmatrix, :([refel.Qr ; refel.Qs ; refel.Qt])));
            end
            
        end
    end
    
    # If coefficients need to be computed, do so
    # # First remove duplicates
    # unique_coef = [];
    # for i=1:length(needed_coef)
    #     already = false;
    #     for j=1:length(unique_coef)
    #         if unique_coef[j] === needed_coef[i]
    #             already = true;
    #         end
    #     end
    #     if !already
    #         push!(unique_coef, needed_coef[i]);
    #     end
    # end
    # needed_coef = unique_coef;
    
    # For constant coefficients, this generates something like:
    ######################################
    # coef_n = a.value[1];
    ######################################
    
    # For variable coefficients, this generates something like:
    ######################################
    # coef_n = zeros(refel.Np);
    # for coefi = 1:refel.Np
    #     coef_n[coefi] = a.value[1].func(x[coefi,1], x[coefi,2],x[coefi,3],time);
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
                tmps = "coef_"*string(i);
                tmpc = Symbol(tmps);
                (ctype, cval) = get_coef_val(needed_coef[i]);
                if ctype == 1
                    # constant coefficient -> coef_n = cval
                    tmpn = cval;
                    push!(code.args, Expr(:(=), tmpc, tmpn)); # allocate coef_n
                elseif ctype == 2
                    # genfunction coefficients -> coef_n = coef.value[1].func(cargs)
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc;
                    tmpn = :(a.value[1]);
                    tmpn.args[1].args[1] = needed_coef[i];
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                elseif ctype == 3
                    # variable values -> coef_n = variable.values
                    tmpb = :(Femshop.variables[$cval].values[gbl]);
                    push!(code.args, Expr(:(=), tmpc, tmpb));
                end
            end
        end
        
        if length(cloopin.args) > 0
            cloop.args[2] = cloopin;
            push!(code.args, cloop); # add loop to code
        end
        
    end
    
    # finally add the code expression
    if length(terms) > 1
        tmp = :(a+b);
        tmp.args = [:+];
        for i=1:length(terms)
            push!(tmp.args, terms[i]);
        end
        terms[1] = tmp;
    end
    push!(code.args, Expr(:return, terms[1]));
    return code;
end

# Dendro version returns a string that will be included in either the feMatrix or feVector
function generate_code_layer_dendro(ex, rhsvar=nothing)
    code = ""; # the code will be in this string
    
    need_derivative = false;
    term_derivative = [];
    needed_coef = [];
    
    to_delete = []; # arrays allocated with new that need deletion
    
    # Separate terms
    terms = separate_terms(ex);
    code_terms = [];
    
    # Process the terms turning them into the code layer
    cindex = 0;
    for i=1:length(terms)
        (codeterm, der, coe) = process_term_dendro(terms[i], cindex, length(terms), i, rhsvar);
        need_derivative = need_derivative || der;
        push!(term_derivative, der);
        append!(needed_coef, coe);
        cindex = length(needed_coef);
        push!(code_terms, codeterm);
    end
    
    # If coefficients need to be computed, do so
    
    # For constant coefficients, this generates something like:
    ######################################
    # double coef_n = a;
    ######################################
    
    # For variable coefficients, this generates something like:
    ######################################
    # double* coef_n = new double[nPe];
    # double val;
    # for(unsigned int i=0; i<nPe; i++){
    #     a(coords[i*3+0], coords[i*3+0], coords[i*3+0], &val);
    #     coef_n[coefi] = val;
    # }
    ######################################
    coef_alloc = "";
    coef_loop = "";
    if length(needed_coef) > 0
        for i=1:length(needed_coef)
            if !(typeof(needed_coef[i]) <: Number || needed_coef[i] === :dt)
                (isconst, cval) = is_constant_coef(needed_coef[i]);
                if isconst
                    # constant coefficient
                    coef_alloc *= "double coef_"*string(i)*" = "*string(cval)*";\n";
                else
                    # variable coefficients
                    coef_alloc *= "double* coef_"*string(i)*" = new double[nPe];\n";
                    push!(to_delete, "coef_"*string(i));
                    coef_loop *= "    "*cval*"(coords[i*3+0], coords[i*3+1], coords[i*3+2], &cval);\n";
                    coef_loop *= "    coef_"*string(i)*"[i] = cval;\n";
                end
            end
        end
        # add the coef loop if needed
        if length(coef_loop) > 1
            coef_loop = "double cval;\nfor(unsigned int i=0; i<nPe; i++){\n"*coef_loop*"}\n";
        end
    end
    
    # If temp storage was needed for multiple terms, allocate
    temp_alloc = "";
    if length(terms) > 1
        for i=1:length(terms)
            temp_alloc *= "double* out_"*string(i)*" = new double[nPe];\n";
            push!(to_delete, "out_"*string(i));
        end
    end
    
    # delete any allocated arrays
    dealloc = "";
    for i=1:length(to_delete)
        dealloc *= "delete [] "*to_delete[i]*";\n";
    end
    
    # Put the pieces together
    code = temp_alloc * coef_alloc * coef_loop;
    for i=1:length(terms)
        code *= code_terms[i] * "\n";
    end
    if length(terms) > 1 # combine the terms' temp arrays
        combine_loop = "out_1[i]";
        for i=2:length(terms)
            combine_loop *= "+out_"*string(i)*"[i]";
        end
        code *= 
"\nfor(unsigned int i=0;i<nPe;i++){
    out[i]="*combine_loop*";
}\n";
    end
    code *= dealloc;
    
    return code;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term_julia(sterm, coef_index, rhsvar)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    weight_part = :(refel.wg .* detJ);
    
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
    
    # check for derivatives
    for j=1:length(factors)
        if factors[j] === :GRADTEST || factors[j] === :GRADTRIAL
            need_derivative = true;
        end
    end
    
    # Separate factors into test/trial/coefficient parts
    coef_facs = [];
    for i=1:length(factors)
        if factors[i] == :TEST || factors[i] == :GRADTEST
            test_part = factors[i];
        elseif factors[i] == :TRIAL || factors[i] == :GRADTRIAL
            trial_part = factors[i];
        else
            push!(coef_facs, factors[i]);
        end
    end
    
    # If rhsvar is a variable symbol, change TRIAL into var.values and treat as a coefficient
    if !(rhsvar === nothing) && trial_part === :TRIAL
        trial_part = nothing;
        tmpv = :(a.values[gbl]);
        tmpv.args[1].args[1] = rhsvar;
        push!(coef_facs, tmpv);
    end
    
    if test_part === :TEST
        test_part = :(refel.Q');
    elseif test_part === :GRADTEST
        test_part = :(transpose(Rxmatrix * Qrmatrix));
    end
    
    if trial_part === :TRIAL
        trial_part = :(refel.Q);
    elseif trial_part === :GRADTRIAL
        trial_part = :(Rxmatrix * Qrmatrix);
    elseif trial_part === nothing
        trial_part = :(refel.Q);
    end
    
    # build coefficient parts
    if length(coef_facs) > 0
        for j=1:length(coef_facs)
            tmp = coef_facs[j];
            if typeof(tmp) == Symbol && !(tmp ===:dt)
                push!(needed_coef, tmp);
                tmps = "coef_"*string(coef_index + length(needed_coef));
                tmp = Symbol(tmps);
            end
            if j>1
                coef_part = :($coef_part .* $tmp);
            else
                coef_part = tmp;
            end
        end
    end
    
    # put together weight/detJ and coefficient parts
    if need_derivative && config.dimension == 2
        weight_part = :(vcat($weight_part, $weight_part));
        if rhsvar === nothing && length(needed_coef) > 0
            coef_part = :(vcat($coef_part, $coef_part));
        end
    elseif need_derivative && config.dimension == 3
        weight_part = :(vcat($weight_part, vcat($weight_part, $weight_part)));
        if rhsvar === nothing && length(needed_coef) > 0
            coef_part = :(vcat($coef_part, vcat($coef_part, $coef_part)));
        end
    end
    
    term = test_part;
    if !(trial_part === nothing)
        term = :($test_part * diagm($weight_part) * $trial_part);
        
        if !(coef_part === nothing)
            if rhsvar === nothing
                term = :($test_part * (diagm($weight_part .* $coef_part) * $trial_part));
            else
                term = :($test_part * (diagm($weight_part) * ($trial_part * $coef_part)));
            end
            
        end
    else
        if !(coef_part === nothing)
            term = :($test_part * $coef_part);
        end
    end
    
    if neg
        negex = :(-a);
        negex.args[2] = copy(term);
        term = negex;
    end
    
    return (term, need_derivative, needed_coef);
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term_dendro(sterm, coef_index, termcount, thisterm, rhsvar)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    
    test_part = "";
    trial_part = "";
    coef_part = "";
    
    # Will the output of this term be stored in "out" or a temp value "out_n"
    out_name = "out";
    if termcount > 1
        out_name = "out_"*string(thisterm);
    end
    
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
    
    # check for derivatives
    for j=1:length(factors)
        if factors[j] === :GRADTEST || factors[j] === :GRADTRIAL
            need_derivative = true;
        end
    end
    
    # Separate factors into test/trial/coefficient parts
    coef_facs = [];
    for i=1:length(factors)
        if factors[i] == :TEST || factors[i] == :GRADTEST
            test_part = factors[i];
        elseif factors[i] == :TRIAL || factors[i] == :GRADTRIAL
            trial_part = factors[i];
        else
            push!(coef_facs, factors[i]);
        end
    end
    
    if test_part === :TEST
        # transpose(Q)
        test_part = 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,"*out_name*");\n";
    elseif test_part === :GRADTEST
        # transpose(Rx*Qr)
        test_part = 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,Qx,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qx);

DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qy,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qy);

DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qz,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,Qz);\n";
    end
    
    if trial_part === :TRIAL
        # Q
        trial_part = 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
    elseif trial_part === :GRADTRIAL
        # Rx*Qr
        trial_part = 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qx);

DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qy);

DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,Qz);\n";
    elseif trial_part == ""
        # Q
        trial_part = 
"DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,"*out_name*");\n";
    end
    
    # build weight/coefficient parts
    idxstr = "[i]";
    weight_part = "W1d[i]";
    if config.dimension == 2
        idxstr = "[j*nrp+i]";
        weight_part = "W1d[i]*W1d[j]";
    elseif config.dimension == 3
        idxstr = "[(k*nrp+j)*nrp+i]";
        weight_part = "W1d[i]*W1d[j]*W1d[k]";
    end
    
    # If term is negative, apply it here
    if neg
        weight_part = "-"*weight_part;
    end
    
    # Multiply by coefficients if needed
    if length(coef_facs) > 0 && rhsvar === nothing
        for j=1:length(coef_facs)
            tmp = coef_facs[j];
            if typeof(tmp) == Symbol && !(tmp ===:dt)
                push!(needed_coef, tmp);
                tmp = "coef_"*string(coef_index + length(needed_coef))*idxstr;
            else
                tmp = string(tmp);
            end
            if j>1
                coef_part = :($coef_part * $tmp);
            else
                coef_part = tmp;
            end
        end
        weight_part = weight_part*" * "*coef_part;
    end
    
    # build the inner weight/coef loop
    body = out_name*"[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=(Jx*Jy*Jz*"*weight_part*");";
    if need_derivative
        body = "Qx[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jy*Jz)/Jx)*"*weight_part*");
            Qy[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jz)/Jy)*"*weight_part*");
            Qz[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jy)/Jz)*"*weight_part*");";
    end
    
    wcloop =
"for(unsigned int k=0;k<(eleOrder+1);k++){
    for(unsigned int j=0;j<(eleOrder+1);j++){
        for(unsigned int i=0;i<(eleOrder+1);i++){
            "*body*"
        }
    }
}\n";
    
    # Put the pieces together
    term = trial_part * "\n" * wcloop * "\n" * test_part;
    if need_derivative
        term *= 
"\n"*"for(unsigned int i=0;i<nPe;i++){
    "*out_name*"[i]=Qx[i]+Qy[i]+Qz[i];
}\n";
    end
    
    return (term, need_derivative, needed_coef);
end

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
# Returns: value, name, or array
function get_coef_val(c)
    type = 0;
    val = 0;
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            isit = (typeof(coefficients[i].value[1]) <: Number);
            if isit
                type = 1;
                val = coefficients[i].value[1];
            else
                type = 2;
                val = coefficients[i].value[1].name;
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

