#=
# Various utilities for the SymbolicParser module
=#

# Is s a test function or operator applied to test function?
function is_test_op(s, t)
    if s === t
        return true;
    end
    if typeof(s) == Expr
        if s.args[1] === :-
            return is_test_op(s.args[2], t);
        end
        for p in test_ops
            if s.args[1] === p && s.args[2] === t
                return true;
            end
        end
    end
    return false;
end

function is_var_op(s, v)
    if s === v
        return true;
    end
    if typeof(s) == Expr
        if s.args[1] === :-
            return is_var_op(s.args[2], v);
        end
        for p in var_ops
            if s.args[1] === p && s.args[2] === v
                return true;
            end
        end
        # Dt is special
        if s.args[1] === :Dt
            return true;
        end
    end
    return false;
end

# example: (1-3*a) + b  ->  [(1-3*a), b]
# a + b + c + (d + e)  ->  [a, b, c, (d+e)]
# assumes minus has been converted to plus(negative)
function get_top_terms(ex)
    if !(typeof(ex) == Expr)
        return [ex];
    end
    terms = [ex];
    if ex.head === :call && ex.args[1] === :+
        terms = ex.args[2:end];
    end
    return terms;
end

# example: (1 + -3*a) + b  ->  [1, -3*a, b]
# assumes minus has been converted to plus(negative)
function get_all_terms(ex)::Array{Any,1}
    if !(typeof(ex) == Expr)
        return [ex];
    end
    terms = [];
    if ex.head === :call && ex.args[1] === :+
        for i=2:length(ex.args)
            terms = [terms; get_all_terms(ex.args[i])];
        end
    else
        terms = [ex];
    end
    return terms;
end

# example: (3*a)*(-b)*c ->  [3, a, -b, c]
# assumes expanded term with only negative and * calls
function get_all_factors(ex, neg=false)
    if ex === nothing
        return [];
    end
    if typeof(ex) <: Number
        if neg
            return -ex;
        else
            return ex;
        end
    end
    if typeof(ex) == Symbol
        if neg
            negf = :(-a);
            negf.args[2] = ex;
            return [negf];
        else
            return [ex];
        end
    end
    facs = [];
    if ex.args[1] === :- && length(ex.args) == 2
        if neg
            facs = get_all_factors(ex.args[2], false);
        else
            facs = get_all_factors(ex.args[2], true);
        end
    elseif ex.args[1] === :*
        for i=2:length(ex.args)
            facs = [facs; get_all_factors(ex.args[i])];
        end
        # remove double negatives
        if neg 
            if typeof(facs[1]) == Expr && facs[1].args[1] === :-
                facs[1] = facs[1].args[2];
            else
                negf = :(-a);
                negf.args[2] = facs[1];
                facs[1] = negf;
            end
        end
        
        # collect all the negatives and apply to the first factor
        n = 0;
        for i=1:length(facs)
            if typeof(facs[i]) == Expr && facs[i].args[1] === :-
                n += 1;
                facs[i] = facs[i].args[2];
            end
        end
        if n%2 == 1
            negf = :(-a);
            negf.args[2] = facs[1];
            facs[1] = negf;
        end
    else
        if neg
            facs = [distribute_negative(ex)];
        else
            facs = [ex];
        end
    end
    return facs;
end

# multiplies all factors
# [3, a, u]  ->  3*a*u
function assemble_term(facs)
    if length(facs) == 1
        return facs[1];
    elseif length(facs) > 1
        term = :(a*b);
        term.args = [:*]
        for i=1:length(facs)
            if !(facs[i] === nothing)
                if facs[i] == 0
                    return 0;
                end
                term.args = [term.args ; facs[i]];
            end
        end
        if length(term.args) > 2
            return term;
        elseif length(term.args) == 2
            return term.args[2];
        else
            return 0;
        end
    else
        return 0;
    end
end

# (3*a*7*u)  ->  (21, a*u)
function extract_constant(ex)
    c = 1;
    v = [];
    facs = get_all_factors(ex);
    for i=1:length(facs)
        if typeof(facs[i]) <: Number
            c = c*facs[1];
        else
            v = [v ; facs[i]];
        end
    end
    v = assemble_term(v);
    if c == 1
        c = nothing;
    end
    
    return (c, v);
end

# (a+b-c)  ->  -a + -b + c
# (a*b)  ->  -a*b
function distribute_negative(ex)
    negex = :(-a);
    if typeof(ex) <: Number
        return -ex;
    elseif typeof(ex) == Symbol
        negex.args[2] = ex;
        return negex;
    elseif typeof(ex) == Expr
        if ex.args[1] === :-
            if length(ex.args) > 2 # -(a-b)  ->  -a + b
                newex = copy(ex);
                newex.args[1] = :+;
                newex.args[2] = distribute_negative(newex.args[2]);
                return newex;
            else # -(-a) ->  a
                return ex.args[2];
            end
        elseif ex.args[1] === :+ # -(a+b)  ->  -a + -b
            newex = copy(ex);
            for i=2:length(ex.args)
                newex.args[i] = distribute_negative(newex.args[i]);
            end
            return newex;
        elseif ex.args[1] === :* #-(a*b)  ->  -a*b
            newex = copy(ex);
            newex.args[2] = distribute_negative(newex.args[2]);
            return newex;
        else # -(?expr) -> -(?expr)
            negex.args[2] = ex;
            return negex;
        end
    end
end

# dot(u+v, w+x)  ->  dot(u,w) + dot(v,w) + ...
function expand_dot(ex)
    left = expand(ex.args[2]);
    right = expand(ex.args[3]);
    if typeof(left) <: Number  || typeof(right) <: Number # shouldn't happen. 
        return 0;
    end
    # 
    lterms = get_all_terms(left);
    rterms = get_all_terms(right);
    plusex = :(a+b);
    dotex = :(dot(a,b));
    if length(lterms) < 1 || length(rterms) < 1 # one of the arguments disappeared(became 0?)
        return 0;
    end
    if length(lterms) == 1 && length(rterms) == 1 # dot(a,b) = dot(a,b)
        dotex.args[2] = lterms[1];
        dotex.args[3] = rterms[1];
        result = dotex;
    else
        plusex.args = [:+];
        for i=1:length(lterms)
            for j=1:length(rterms);
                dotex.args[2] = lterms[i];
                dotex.args[3] = rterms[j];
                plusex.args = [plusex.args ; copy(dotex)];
            end
        end
        result = plusex;
    end
    
    return result;
end

# Takes an expression like grad(ex)
# grad(au - bv) = grad(a)u + a*grad(u) + -grad(b)v + -b*grad(v)
function expand_grad(ex)
    inside = expand(ex.args[2]);
    newex = :(grad(a));
    if typeof(inside) <: Number # grad(4) = 0
        return 0;
    elseif typeof(inside) == Symbol # grad(a) = grad(a)
        return ex;
    elseif typeof(inside) == Expr
        if inside.args[1] === :- && length(inside.args) == 2
            inside = inside.args[2];
            newex.args[2] = inside;
            newex = expand_grad(newex);
            if newex == 0
                return 0;
            end
            newex = distribute_negative(newex);
            return newex;
        elseif inside.args[1] === :+
            plusex = copy(inside);
            for i=2:length(inside.args)
                gradex = :(grad(a));
                gradex.args[2] = inside.args[i];
                plusex.args[i] = expand_grad(gradex);
            end
            return plusex;
        elseif inside.args[1] === :*
            plusex = copy(inside);
            plusex.args[1] = :+;
            mulex = copy(inside);
            gradex = :(grad(a));
            for i=2:length(inside.args);
                gradex.args[2] = inside.args[i];
                mulex.args[i] = expand_grad(gradex);
                plusex.args[i] = copy(mulex);
                mulex.args[i] = inside.args[i];
            end
            return plusex;
        end
    end
    return ex;
end

# div(5*a*u + b)  ->  (dot(5*grad(a),u) + a*div(u)) + div(b)
# assume the last factor is the vector (above u and b are vectors)
function expand_div(ex)
    inside = ex.args[2];
    if typeof(inside) <: Number # div(4) = 0
        return 0;
    elseif typeof(inside) == Symbol # div(a) = div(a)
        return ex;
    elseif typeof(inside) == Expr
        if inside.args[1] === :- && length(inside.args) == 2 # div(-u) = -div(u)
            inside = inside.args[2];
            newex = :(div(a));
            newex.args[2] = inside;
            newex = expand_div(newex);
            if newex == 0
                return 0;
            end
            newex = distribute_negative(newex);
            return newex;
        elseif inside.args[1] === :+ # div(a+b) = div(a) + div(b)
            plusex = copy(inside);
            for i=2:length(inside.args)
                divex = :(div(a));
                divex.args[2] = inside.args[i];
                plusex.args[i] = expand_div(divex);
            end
            return plusex;
        elseif inside.args[1] === :* # div(a*b*c*u) = dot(grad(a*b*c),u) + a*b*c*div(u)
            mulex = copy(inside);
            gradex = :(grad(a));
            divex = :(div(a));
            dotex = :(dot(a,b));
            plusex = :(a+b);
            
            vect = mulex.args[end];
            if length(mulex.args) > 3
                scal = copy(mulex);
                scal.args = scal.args[1:end-1];
            else
                scal = mulex.args[2];
            end
            
            gradex.args[2] = scal;
            gradex = expand_grad(gradex);
            divex.args[2] = vect;
            divex = expand_div(divex);
            dotex.args[2] = gradex;
            dotex.args[3] = vect;
            dotex = expand_dot(dotex);
            mulex.args[end] = divex;
            plusex.args[2] = dotex;
            plusex.args[3] = mulex;
            
            return plusex;
        end
    end
    
    return ex;
end

# (Dt(2*a*u*v))  ->  2*a*Dt(u)*v
function handle_dt(ex, var)
    inside = ex.args[2];
    facs = get_all_factors(inside);
    out = [];
    in = [];
    for i=1:length(facs)
        if is_var_op(facs[i], var)
            in = [in ; facs[i]];
        else
            out = [out ; facs[i]];
        end
    end
    if length(in) == 0 # this shouldn't happen
        return assemble_term(out);
    end
    if length(out) == 0
        term = :(Dt(a));
        term.args[2] = assemble_term(in);
        return term;
    else
        term = :(a*Dt(b));
        term.args[2] = assemble_term(out);
        term.args[3].args[2] = assemble_term(in);
        return term;
    end
end

# expand a polynomial expression and replace minus with +(-1*) and divide with *(arg^(-1))
# ((1-3*r)*u + grad(u))*v  ->  1*u*v - 3*r*u*v + grad(u)*v
function expand(ex)
    if !(typeof(ex) == Expr)
        return ex;
    end
    newex = copy(ex);
    if newex.head === :call
        if newex.args[1] === :+     # recursively expand args
            for i=2:length(newex.args)
                newex.args[i] = expand(newex.args[i]);
            end
        elseif newex.args[1] === :-
            if length(newex.args) > 2 # swap minus with plus(-) and expand args
                minusex = :(-a);
                for i=3:length(newex.args)
                    minusex.args[2] = newex.args[i];
                    newex.args[i] = copy(minusex);
                end
                newex.args[1] = :+;
                newex = expand(newex);
            else # -(ex)  ->  -(expand(ex))
                newex = distribute_negative(expand(newex.args[2]));
            end
        elseif newex.args[1] === :* # expand each arg and distribute terms (a+b)*(c+d)  ->  ac+bc+ad+bd
            for i=2:length(newex.args) # first expand each factor
                newex.args[i] = expand(newex.args[i]);
            end
            newargs = get_top_terms(newex.args[2]);
            grownargs = [];
            grownargs = [grownargs; copy(newargs)]; 
            for i=3:length(newex.args)
                tempargs = get_top_terms(newex.args[i]);
                mulex = :(a*b);
                leftN = length(newargs);
                
                for j=2:length(tempargs) # grow array to match needed number of terms
                    grownargs = vcat(grownargs, copy(newargs));
                end
                newargs = grownargs;
                
                for j=1:leftN # multiply left terms by right terms
                    mulex.args[2] = newargs[j];
                    for k=1:length(tempargs)
                        mulex.args[3] = tempargs[k];
                        newargs[(k-1)*leftN+j] = copy(mulex);
                    end
                end
            end
            if length(newargs) > 1
                newex.args = [:+ ; newargs]; # replace the expression with the expanded one
            else
                newex = newargs[1];
            end
            
        elseif newex.args[1] === :/ # replace (/) with(*()^(-1)) and expand
            powex = :(a^(-1));
            powex.args[2] = newex.args[3];
            newex.args[1] = :*;
            newex.args[3] = powex;
            newex = expand(newex);
        elseif newex.args[1] === :grad # use expand_grad to expand this 
            newex = expand_grad(newex);
        elseif newex.args[1] === :div # use expand_div to expand this 
            newex = expand_div(newex);
        end
        # TODO consider powers: (a+b)^2  ->  aa+ab+ba+bb
        
    end
    
    # Restructure the expression into one long +() call
    # Merge terms into single *() or -(*()) calls
    if typeof(newex) == Expr && length(newex.args) > 2
        terms = get_all_terms(newex);
        mulex = :(a*b);
        for i=1:length(terms)
            facs = get_all_factors(terms[i]);
            # If any of the factors is 0, set this whole term to zero.
            for j=1:length(facs)
                if facs[j] == 0
                    facs = [:(0)];
                    break;
                end
            end
            if length(facs) > 1
                mulex.args = [:* ; facs];
                terms[i] = copy(mulex);
            else
                terms[i] = facs[1];
            end
        end
        # remove any zero terms
        nzterms = [];
        for i=1:length(terms)
            if !(terms[i] == 0 || terms[i] === :(0) || terms[i] === :(-0))
                nzterms = [nzterms; terms[i]];
            end
        end
        if length(nzterms) > 1
            newex = :(a+b);
            newex.args = [:+ ; nzterms];
        elseif length(nzterms) == 1
            newex = nzterms[1];
        else
            newex = 0;
        end
    end
    
    return newex;
end