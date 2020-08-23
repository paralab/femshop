#=
Use the symbolic layer expressions to generate the FEM code
=#
function generate_code_layer(ex, var, lorr)
    if language == 0 || language == JULIA
        return generate_code_layer_julia(ex, var, lorr);
    elseif language == CPP
        return generate_code_layer_dendro(ex, var, lorr);
    elseif language == MATLAB
        return generate_code_layer_homg(ex, var, lorr);
    end
end

# External code gen in these similar files
include("generate_code_layer_dendro.jl");
include("generate_code_layer_homg.jl");

# Julia and utils are in this file

###############################################################################################################
# julia
###############################################################################################################

# Julia version returns an expression for the generated function for linear or bilinear term
function generate_code_layer_julia(symex, var, lorr)
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
    for i=1:length(terms)
        (codeterm, der, coe, coeind, testi, trialj) = process_term_julia(terms[i], var, lorr);
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
        terms[i] = codeterm;
    end
    
    # If derivatives are needed, prepare the appropriate matrices
    if need_derivative
        if config.dimension == 1
            push!(code.args, Expr(:(=), :R1matrix, :(diagm(J.rx))));
            push!(code.args, Expr(:(=), :Q1matrix, :(refel.Qr)));
        elseif config.dimension == 2
            push!(code.args, Expr(:(=), :R1matrix, :([diagm(J.rx) diagm(J.sx)])));
            push!(code.args, Expr(:(=), :Q1matrix, :([refel.Qr ; refel.Qs])));
            push!(code.args, Expr(:(=), :R2matrix, :([diagm(J.ry) diagm(J.sy)])));
            push!(code.args, Expr(:(=), :Q2matrix, :([refel.Qr ; refel.Qs])));
        elseif config.dimension == 3
            push!(code.args, Expr(:(=), :R1matrix, :([diagm(J.rx) diagm(J.sx) diagm(J.tx)])));
            push!(code.args, Expr(:(=), :Q1matrix, :([refel.Qr ; refel.Qs ; refel.Qt])));
            push!(code.args, Expr(:(=), :R2matrix, :([diagm(J.ry) diagm(J.sy) diagm(J.ty)])));
            push!(code.args, Expr(:(=), :Q2matrix, :([refel.Qr ; refel.Qs ; refel.Qt])));
            push!(code.args, Expr(:(=), :R3matrix, :([diagm(J.rz) diagm(J.sz) diagm(J.tz)])));
            push!(code.args, Expr(:(=), :Q3matrix, :([refel.Qr ; refel.Qs ; refel.Qt])));
        end
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
                if cind < 0
                    # probably a variable
                    cind = string(needed_coef[i]);
                end
                tmps = "coef_"*string(cind)*"_"*string(needed_coef_ind[i]);
                tmpc = Symbol(tmps);
                (ctype, cval) = get_coef_val(needed_coef[i], needed_coef_ind[i]);
                if ctype == 1
                    # constant coefficient -> coef_n = cval
                    tmpn = cval;
                    push!(code.args, Expr(:(=), tmpc, tmpn)); # allocate coef_n
                elseif ctype == 2
                    # genfunction coefficients -> coef_n_i = coef.value[i].func(cargs)
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc;
                    tmpn = :(a.value[1]);
                    tmpn.args[2] = needed_coef_ind[i];
                    tmpn.args[1].args[1] = needed_coef[i];
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                elseif ctype == 3
                    # variable values -> coef_n = variable.values
                    tmpb = :(Femshop.variables[$cval].values[gbl]); #TODO this only works for scalars!
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
    dofsper = 1;
    if typeof(var) <: Array
        for vi=1:length(var)
            dofsper += length(var[vi].symvar.vals); # The number of components for this variable
        end
    elseif !(var.type == SCALAR)
        dofsper = length(var.symvar.vals);
    end
    
    if dofsper > 1
        if lorr == RHS
            push!(code.args, Expr(:(=), :full_vector, :(zeros(refel.Np*$dofsper)))); # allocate vector
        else
            push!(code.args, Expr(:(=), :full_matrix, :(zeros(refel.Np*$dofsper, refel.Np*$dofsper)))); # allocate matrix
        end
    end
    
    if length(terms) > 1
        if typeof(var)==Variable # Only one variable
            if var.type == SCALAR # Only one component
                tmp = :(a+b);
                tmp.args = [:+];
                for i=1:length(terms)
                    push!(tmp.args, terms[i]);
                end
                terms[1] = tmp;
            else # More than one component
                # add terms into full matrix according to testind/trialind
                tmp = :(a+b);
                tmp.args = [:+];
                for i=1:length(terms)
                    ti = test_ind[i][1]-1;
                    #rows = :(($ti * refel.Np + 1):(($ti + 1)*refel.Np));
                    tj = trial_ind[i][1]-1;
                    #cols = :(($tj * refel.Np + 1):(($tj + 1)*refel.Np));
                    
                    if lorr == LHS
                        push!(code.args, Expr(:(+=), :(full_matrix[($ti * refel.Np + 1):(($ti + 1)*refel.Np),($tj * refel.Np + 1):(($tj + 1)*refel.Np)]), terms[i]));
                    else
                        push!(code.args, Expr(:(+=), :(full_vector[($ti * refel.Np + 1):(($ti + 1)*refel.Np)]), terms[i]));
                    end
                    
                end
                if lorr == LHS
                    terms[1] = :full_matrix;
                else
                    terms[1] = :full_vector;
                end
                
            end
        else #more than one variable
            #TODO
        end
    end
    
    # At this point everything is packed into terms[1]
    push!(code.args, Expr(:return, terms[1]));
    return code;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term_julia(sterm, var, lorr)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    needed_coef_ind = [];
    
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    weight_part = :(refel.wg .* detJ);
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
        (index, v, mods) = extract_symbols(factors[i]);
        
        if is_test_func(v)
            test_component = index; # the vector index
            if length(mods) > 0
                # TODO more than one derivative mod
                need_derivative = true;
                dmatr = Symbol("R"*mods[1][2]*"matrix");
                dmatq = Symbol("Q"*mods[1][2]*"matrix");
                test_part = :(transpose($dmatr * $dmatq));
            else
                # no derivative mods
                test_part = :(refel.Q');
            end
        elseif is_unknown_var(v, var)
            if !(trial_part === nothing)
                # Two unknowns multiplied in this term. Nonlinear. abort.
                printerr("Nonlinear term. Code layer incomplete.");
                return (-1, -1, -1, -1, -1, -1);
            end
            trial_component = index;
            if length(mods) > 0
                # TODO more than one derivative mod
                need_derivative = true;
                dmatr = Symbol("R"*mods[1][2]*"matrix");
                dmatq = Symbol("Q"*mods[1][2]*"matrix");
                trial_part = :($dmatr * $dmatq);
            else
                # no derivative mods
                trial_part = :(refel.Q);
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
    if lorr == RHS && trial_part === :(refel.Q)
        tmpv = :(a.values[gbl]);
        tmpv.args[1].args[1] = var.symbol; #TODO, will not work for var=array
        push!(coef_facs, tmpv);
        push!(coef_inds, trial_component);
    end
    
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
                    tmps = "coef_"*string(ind)*"_"*string(coef_inds[j]);
                    tmp = Symbol(tmps);
                else
                    tmps = "coef_"*string(tmp)*"_"*string(coef_inds[j]);
                    tmp = Symbol(tmps);
                end
            end
            if j>1
                coef_part = :($coef_part .* $tmp);
            else
                coef_part = tmp;
            end
        end
    end
    
    term = test_part;
    if !(coef_part === nothing)
        if lorr == LHS
            term = :($test_part * (diagm($weight_part .* $coef_part) * $trial_part));
        else # RHS
            term = :($test_part * (diagm($weight_part) * ($trial_part * $coef_part)));
        end
        
    else
        term = :($test_part * diagm($weight_part) * $trial_part);
    end
    
    if neg
        negex = :(-a);
        negex.args[2] = copy(term);
        term = negex;
    end
    
    return (term, need_derivative, needed_coef, needed_coef_ind, test_component, trial_component);
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
                val = coefficients[i].value[comp].name;
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
                    index = [parse(Int, str[i]); index] # The indices on the variable
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
