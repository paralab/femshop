#=
# The final symbolic form will be an Expr expression tree.
# The leaf nodes are all SymEntitys as described below.
# The parents will be math operators. 
# All used entities will also be kept in a list.
#
# A SymEntity has
# - name, the name of this variable, coefficient, etc. (or number for constants)
# - index, the component of this, such as a vector component (-1 for constants)
# - derivs
# - flags
=#
export SymExpression, SymEntity, build_symexpressions

struct SymExpression
    tree                # The expression tree
    entities::Array     # The SymEntitys corresponding to the leaf nodes of the tree
end

struct SymEntity
    name::Union{Float64, String}    # The string or number 
    index::Int64                    # vector component, 1 for scalars, -1 for numbers
    derivs::Array                   # any derivatives applied
    flags::Array                    # any non-derivative modifiers or flags
end

# For printing, write the entity in a more readable format
Base.show(io::IO, x::SymEntity) = print(io, symentity_string(x));
function symentity_string(a::SymEntity)
    result = "";
    for i=1:length(a.derivs)
        result *= "D"*string(a.derivs[i]) * "(";
    end
    
    for i=1:length(a.flags)
        result *= a.flags[i];
    end
    
    result *= a.name;
    if a.index > 0
        result *= "_"*string(a.index);
    end
    
    for i=1:length(a.derivs)
        result *= ")";
    end
    
    return result;
end

function build_symexpressions(var, lhsvol, rhsvol, lhssurf=nothing, rhssurf=nothing)
    lhsvolsym = build_symexpression(var, lhsvol, "L", "V");
    rhsvolsym = build_symexpression(var, rhsvol, "R", "V");
    lhssurfsym = build_symexpression(var, lhssurf, "L", "S");
    rhssurfsym = build_symexpression(var, rhssurf, "R", "S");
    
    if lhssurf === nothing && rhssurf === nothing
        return (lhsvolsym, rhsvolsym);
    else
        return (lhsvolsym, rhsvolsym, lhssurfsym, rhssurfsym);
    end
end

function build_symexpression(var, ex, lorr, vors)
    if ex === nothing
        return nothing
    end
    
    if typeof(ex) <: Array
        exarray = [];
        for i=1:length(ex)
            push!(exarray, build_symexpression(var, ex[i], lorr, vors));
        end
        return exarray;
        
    elseif typeof(ex) <: SymEngine.Basic
        # ex is a SymEngine.Basic object
        # First parse it into an Expr
        jex = Meta.parse(string(ex));
        
        # Then traverse it to find SymTerm leaves and build a SymExpression
        (newex, sents) = extract_symentities(var, jex);
        
        # Change all math operators to broadcast(dot) versions
        newex = broadcast_all_ops(newex);
        
        # If using FEM, reorder factors for easier generation
        if config.solver_type == CG || config.solver_type == DG
            newex = order_expression_for_fem(var, test_functions, newex);
        end
        
        symex = SymExpression(newex, sents);
        
        return symex;
        
    end
end

# Traverse the tree and replace appropriate symbols with SymEntitys.
# Also, keep a list of symentities used.
function extract_symentities(var, ex)
    sents = [];
    if typeof(ex) == Expr
        # Recursively traverse the tree
        newex = copy(ex);
        if ex.head === :call # an operation like +(a,b)
            for i=2:length(ex.args)
                (subex, subsents) = extract_symentities(var, ex.args[i]);
                newex.args[i] = subex;
                append!(sents, subsents);
            end
        elseif ex.head === :if # a conditional like if a<b c else d end
            # TODO
        else
            # What else could it be?
        end
        
    elseif typeof(ex) == Symbol
        sen = build_symentity(var, ex);
        push!(sents, sen);
        newex = sen;
        
    else # probably a number
        newex = ex;
    end
    
    return (newex, sents);
end

# Builds a SymEntity struct from a symbolic name such as:
# DGSIDE1_D2__u_1 -> name=u, index=1, derivs=D2, flags=DGSIDE1
function build_symentity(var, ex)
    if typeof(ex) <: Number
        return SymEntity(ex, 0, [], []);
    elseif typeof(ex) == Symbol
        # Extract the index, derivs, flags
        (ind, symb, mods) = extract_entity_parts(ex);
        # Separate the derivative mods from the other flags
        (derivs, others) = get_derivative_mods(mods);
        
        return SymEntity(symb, ind, derivs, others);
        
    else
        printerr("unexpected type in build_symentity: "*string(typeof(ex))*" : "*string(ex));
    end
end

# Extract meaning from the symbolic object name
# The format of the input symbol should look like this
#   MOD1_MOD2_..._var_n
# There could be any number of mods, _var_n is the symvar symbol (n is the vector index, or 1 for scalar)
# Returns (n, var, [MOD1,MOD2,...]) all as strings
function extract_entity_parts(ex)
    str = string(ex);
    index = [];
    var = nothing;
    mods = [];
    l = lastindex(str);
    e = l; # end of variable name
    b = l; # beginning of variable name
    
    # dt is a special symbol that will be passes as a number value in the generated function.
    if str == "dt"
        return(-1, str, []);
    end
    
    # Extract the name and index
    for i=l:-1:0
        if e==l
            if str[i] == '_'
                e = i-1;
            else
                # These are the indices at the end. Parse one digit at a time.
                try
                    index = [parse(Int, str[i]); index] # The indices on the variable
                catch
                    return ([0],ex,[]);
                end
                
            end
        elseif b==l && i>0
            if str[i] == '_'
                b = i+1;
            end
            
        else
            # At this point we know b and e
            if var === nothing
                var = string((SubString(str, b, e)));
                b = b-2;
                e = b;
            end
        end
    end
    
    # Change index from [a, b] to a*d + b
    newindex = index[end];
    for i=1:(length(index)-1)
        newindex = newindex + (index[i]-1)*config.dimension^(length(index)-i);
    end
    index = newindex;
    
    # extract the modifiers like D1_ etc. separated by underscores
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
    
    return (index, var, mods);
end

# Separate the derivative mods(Dn_) from the other mods and return an array 
# of derivative indices and the remaining mods.
function get_derivative_mods(mods)
    derivs = [];
    others = [];
    
    for i=1:length(mods)
        if length(mods[i]) == 2 && mods[i][1] == 'D'
            try
                index = parse(Int, mods[i][2])
                push!(derivs, index);
            catch
                printerr("Unexpected modifier: "*mods[i]*" see get_derivative_mods() in symexpression.jl")
                push!(others, mods[i]);
            end
        else
            push!(others, mods[i]);
        end
    end
    
    return (derivs, others);
end

# changes all ops to broadcast versions (:+ -> :.+, sqrt(a) -> Expr(:., :sqrt, :((a,))))
function broadcast_all_ops(ex)
    if typeof(ex) == Expr
        newex = copy(ex);
        
        
    else
        newex = ex;
    end
    
    return newex;
end

# Order the factors in each term so that the test function is first and the unknown(when present) is second
function order_expression_for_fem(var, test, ex)
    # Assuming the upper levels of the tree look like t1+t2-t3+...
    # also possible: t1+(t2-t3)+... or just t1
    # Handle recursively in case there is something like t1+(t2+(t3+(t4...)))
    if typeof(ex) == Expr
        newex = copy(ex);
        if newex.head === :call
            if newex.args[1] === :+ || newex.args[1] === :.+ || newex.args[1] === :- || newex.args[1] === :.-
                # It is a + or - operation. Recursively order each term
                for i=2:length(newex.args)
                    newex.args[i] = order_expression_for_fem(var, test, newex.args[i])
                end
            elseif newex.args[1] === :* || newex.args[1] === :.*
                # search for the test function
                test_part = 0;
                var_part = 0;
                for i=2:length(newex.args)
                    if has_specific_part(newex.args[i], test)
                        test_part = i;
                    elseif has_specific_part(newex.args[i], var)
                        var_part = i;
                    end
                end
                
                # Place test first, var second
                if test_part > 2
                    tmp = newex.args[2];
                    newex.args[2] = newex.args[test_part];
                    newex.args[test_part] = tmp;
                    if var_part == 2
                        var_part = test_part;
                    end
                    if var_part > 3
                        tmp = newex.args[3];
                        newex.args[3] = newex.args[var_part];
                        newex.args[var_part] = tmp;
                    end
                    
                elseif var_part > 2
                    tmp = newex.args[2];
                    newex.args[2] = newex.args[var_part];
                    newex.args[var_part] = tmp;
                end
            end
            
        end
        
        return newex;
        
    else # ex is not an Expr
        # It is a single symentity. Do nothing.
        return ex;
    end
    
end

# Returns true if a.name matches any of the parts
function has_specific_part(a, parts)
    if typeof(a) == SymEntity
        if typeof(parts) <: Array
            for i=1:length(parts)
                if has_specific_part(a, parts[i])
                    return true;
                end
            end
        elseif typeof(parts) == Variable
            return a.name == string(parts.symbol);
        elseif typeof(parts) == Coefficient
            return a.name == string(parts.symbol);
        end
    end
    return false;
end