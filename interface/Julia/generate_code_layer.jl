#=
Use the symbolic layer expressions to generate the FEM code
=#

function generate_code_layer(ex, rhsvar=nothing)
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
        (codeterm, der, coe) = process_term(terms[i], cindex, rhsvar);
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
                if is_constant_coef(needed_coef[i])
                    # constant coefficient
                    tmpn = :(a.value[1]);
                    tmpn.args[1].args[1] = needed_coef[i];
                    push!(code.args, Expr(:(=), tmpc, tmpn)); # allocate coef_n
                else
                    # variable coefficients
                    tmpv = :(a[coefi]);
                    tmpv.args[1] = tmpc;
                    tmpn = :(a.value[1]);
                    tmpn.args[1].args[1] = needed_coef[i];
                    tmpb = :(a.func());
                    tmpb.args[1].args[1]= tmpn;
                    append!(tmpb.args, cargs);
                    push!(code.args, Expr(:(=), tmpc, :(zeros(refel.Np)))); # allocate coef_n
                    push!(cloopin.args, Expr(:(=), tmpv, tmpb)); # add it to the loop
                end
            end
        end
        cloop.args[2] = cloopin;
        push!(code.args, cloop); # add loop to code
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

# checks if the coefficient has constant value
function is_constant_coef(c)
    isit = false;
    for i=1:length(coefficients)
        if c === coefficients[i].symbol
            isit = (typeof(coefficients[i].value[1]) <: Number);
        end
    end
    
    return isit;
end

# Changes the symbolic layer term into a code layer term
# also records derivative and coefficient needs
function process_term(sterm, coef_index, rhsvar)
    term = copy(sterm);
    need_derivative = false;
    needed_coef = [];
    
    test_part = nothing;
    trial_part = nothing;
    weight_part = :(refel.wg .* detJ);
    coef_part = nothing;
    
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
    
    #TEST -> transpose(Q)
    #GRADTEST -> transpose(RxQr)
    #TRIAL -> Q
    #GRADTRIAL -> RxQr
    #coef -> diagm(coef.*W.*detJ)
    
    if language == 0 # No external generation
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
        end
        
    elseif language == CPP
        
    elseif language == MATLAB
        
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