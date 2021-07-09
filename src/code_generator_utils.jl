########################################################
# Utility functions used by the various code generators
########################################################

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
# special: type=-1, val=0
# number: type=0, val=number
# constant coef: type=1, val=number
# genfunction: type=2, val= index in genfunctions array
# variable: type=3, val=index in variables array
function get_coef_val(c)
    if c.index == -1
        # It is a special symbol like dt
        return (-1, 0);
    elseif c.index == 0
        # It is a number
        return(0, c.name);
    end
    
    type = 0;
    val = 0;
    for i=1:length(coefficients)
        if c.name == string(coefficients[i].symbol)
            isit = (typeof(coefficients[i].value[c.index]) <: Number);
            if isit
                type = 1; # a constant wrapped in a coefficient ... not ideal
                val = coefficients[i].value[c.index];
            else
                type = 2; # a function
                name = coefficients[i].value[c.index].name;
                for j=1:length(genfunctions)
                    if name == genfunctions[j].name
                        val = j;
                    end
                end
            end
        end
    end
    if type == 0 # a variable
        for i=1:length(variables)
            if c.name == string(variables[i].symbol)
                type = 3;
                val = variables[i].index;
            end
        end
    end
    
    return (type, val);
end

# Returns the index in femshop's coefficient array
# or -1 if it is not there.
function get_coef_index(c)
    ind = -1;
    for i=1:length(coefficients)
        if c.name == string(coefficients[i].symbol)
            ind = coefficients[i].index;
        end
    end
    
    return ind;
end

# Makes a name for this value based on coefficient or variable name, derivative and other modifiers, and component index.
function make_coef_name(c)
    tag = "";
    for i=1:length(c.flags)
        tag *= c.flags[i];
    end
    for i=1:length(c.derivs)
        tag *= "D"*string(c.derivs[i]);
    end
    
    str = "coef_"*tag*"_"*string(c.name)*"_"*string(c.index);
    return str;
end

function is_test_function(ent)
    for t in test_functions
        if string(t.symbol) == ent.name
            return true;
        end
    end
    return false;
end

function is_unknown_var(ent, vars)
    if typeof(vars) == Variable
        return ent.name == string(vars.symbol);
    elseif typeof(vars) <:Array
        for v in vars
            if string(v.symbol) == ent.name
                return true;
            end
        end
    end
    return false;
end

# They are the same if all parts of them are the same.
function is_same_entity(a, b)
    return a.name == b.name && a.index == b.index && a.derivs == b.derivs && a.flags == b.flags;
end

# symex could be an array. If so, make a similar array containing an array of terms.
# If symex is a SymExpression, make an array of terms.
# A term is an Expr or Symbol representing tn in the top level t1+t2+t3...
# The SymEntities in the Expr have been replaced with their corresponding make_coef_name symbols
function process_terms(symex)
    if typeof(symex) <: Array
        terms = [];
        for i=1:length(symex)
            push!(terms, process_terms(symex[i]));
        end
        
    else
        if typeof(symex) == SymExpression
            newex = copy(symex.tree);
        elseif typeof(symex) == Expr
            newex = copy(symex);
        else
            newex = symex;
        end
        
        #println(string(newex) * " : " * string(typeof(newex)))
        
        if typeof(newex) == Expr && newex.head === :call
            if (newex.args[1] === :+ || newex.args[1] === :.+ || newex.args[1] === :- || newex.args[1] === :.-)
                # each arg is a term
                # Do this recursively to handle t1+(t2-(t3+...))
                terms = [];
                for i=2:length(newex.args)
                    append!(terms, process_terms(newex.args[i]));
                end
            else
                if newex.args[1] === :/ || newex.args[1] === :./
                    # change a/b to a*1/b
                    mulex = :(a.*b);
                    mulex.args[2] = newex.args[2];
                    newex.args[2] = 1;
                    mulex.args[3] = newex;
                    newex = mulex;
                end
                terms = newex;
            end
        else
            # There is just one term
            terms = newex;
        end
    end
    
    return terms;
end

# Assume the term looks like a*b*c or a*(b*c) or something similar.
function separate_factors(ex, var=nothing)
    test_part = nothing;
    trial_part = nothing;
    coef_part = nothing;
    
    if typeof(ex) == Expr
        if (ex.args[1] === :.- || ex.args[1] === :-) && length(ex.args)==2
            # a negative sign. Stick it on one of the factors
            (test_part, trial_part, coef_part) = separate_factors(ex.args[2], var);
            negex = :(-a);
            if !(coef_part === nothing)
                negex.args[2] = coef_part;
                coef_part = negex;
            elseif !(test_part === nothing)
                negex.args[2] = test_part;
                test_part = negex;
            elseif !(trial_part === nothing)
                negex.args[2] = trial_part;
                trial_part = negex;
            end
            
        elseif ex.args[1] === :.* || ex.args[1] === :*
            tmpcoef = [];
            # Recursively separate each arg to handle a*(b*(c*...))
            for i=2:length(ex.args)
                (testi, triali, coefi) = separate_factors(ex.args[i], var);
                if !(testi === nothing)
                    test_part = testi;
                end
                if !(triali === nothing)
                    trial_part = triali;
                end
                if !(coefi === nothing)
                    push!(tmpcoef, coefi);
                end
            end
            if length(tmpcoef) > 1
                coef_part = :(a.*b);
                coef_part.args = [:.*];
                append!(coef_part.args, tmpcoef);
            elseif length(tmpcoef) == 1
                coef_part = tmpcoef[1];
            end
            
        end
    elseif typeof(ex) == SymEntity
        if is_test_function(ex)
            test_part = ex;
        elseif !(var === nothing) && is_unknown_var(ex, var)
            trial_part = ex;
        else
            coef_part = ex;
        end
    else # a number?
        coef_part = ex;
    end
    
    # println(ex)
    # println(test_part)
    # println(trial_part)
    # println(coef_part)
    # println("")
    
    return (test_part, trial_part, coef_part);
end

function replace_entities_with_symbols(ex)
    if typeof(ex) <: Array
        for i=1:length(ex)
            ex[i] = replace_entities_with_symbols(ex[i]);
        end
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = replace_entities_with_symbols(ex.args[i]);
        end
    elseif typeof(ex) == SymEntity
        return Symbol(make_coef_name(ex));
    
    end
    
    return ex;
end