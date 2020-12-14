#= 
# A set of tools for parsing the variational forms into symEngine expressions.
=#
module SymbolicParser
export sp_parse, add_custom_op, add_custom_op_file

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, NEUMANN, ROBIN, NO_BC,
        MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR, SYM_TENSOR,
        LHS, RHS,
        LINEMESH, QUADMESH, HEXMESH
import ..Femshop: Femshop_config, Femshop_prob, Variable, Coefficient
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, coefficients, parameters, test_functions

using SymEngine, LinearAlgebra

#### globals ########################
# Basic symbolic operators that are included automatically and have precedence
op_names = [];
ops = [];

# Custom operators that can be added by a user
custom_op_names = [];
custom_ops = [];
#####################################

include("symtype.jl");
include("symoperator.jl");
include("basic_ops.jl");

# a special operator for dealing with scalar multiplication when the scalar is in an array
import Base.*
function *(a::Array{Basic,1}, b::Array{Basic,1})
    if size(a) == (1,)
        return b.*a[1];
    elseif size(b) == (1,)
        return a.*b[1];
    elseif size(a) == size(b)
        # This should be an error, but I will treat it as a dot product for lazy people
        return [transpose(a) * b];
    else
        # un oh... How do I redirect this to the symengine method?
        # (sigh) This will be an error until I figure it out.
    end
end
# import Base.sqrt
# function sqrt(a::Array{Basic,1})
#     return a .^ (1//2);
# end

# Adds a single custom operator
function add_custom_op(s, handle)
    global custom_op_names;
    global custom_ops;
    push!(custom_op_names, s); # add the Symbol s to the list
    push!(custom_ops, SymOperator(s, handle));
end

# Includes a set of custom operators
# The file will include an array of symbols and an array of function handles
function add_custom_op_file(file)
    include(file);
    for i=1:length(_names)
        add_custom_op(_names[i], _handles[i]);
        log_entry("Added custom operator: "*string(_names[i]));
    end
end

# Parses a variational form expression into a SymEngine Basic expression 
# Takes an Expr, variable symbol, and test function symbol
# The symbols are only for determining lhs vs. rhs
# Returns a tuple of SymEngine expressions (lhs, rhs) for the equation lhs = rhs
# lhs contains terms including the unknown variable
# rhs contains terms without it
function sp_parse(ex, var)
    debug = false;
    lhs = nothing;
    rhs = nothing;
    varcount = 1;
    timederiv = false;
    if debug println("expr = "*string(ex)); end
    
    # Check that there are as many vars as exs
    if typeof(var) <: Array
        if !(typeof(ex) <:Array && length(var) == length(ex))
            printerr("Error: Need same # of unknowns and equations");
            return (nothing,nothing);
        end
        varcount = length(var);
    end
    
    # Work with a copy of ex
    symex = copy(ex);
    
    # Insert parameters
    symex = insert_parameters(symex);
    if debug println("insert parameters -> "*string(symex)); end
    
    # Replace symbols for variables, coefficients, test functions, and special operators
    symex = replace_symbols(symex);
    if debug println("replace symbols -> "*string(symex)); end
    
    # Evaluate the expression to apply symbolic operators
    symex = apply_ops(symex);
    if debug println("apply ops -> "*string(symex)); end
    
    # Expand the expression and separate terms
    sterms = get_sym_terms(symex);
    if debug println("sterms = "*string(sterms)); end
    
    # Check for time derivatives
    timederiv = check_for_dt(sterms);
    
    # Check for surface integrals
    has_surface = check_for_surface(sterms);
    
    # Each element has an lhs and rhs
    lhs = copy(sterms); # set up the container right
    rhs = copy(sterms);
    if varcount > 1
        for i=1:length(symex)
            sz = size(symex[i]);
            (lhs[i],rhs[i]) = split_left_right(sterms[i],sz,var);
            if length(rhs[i]) == 0
                rhs[i] = [Basic(0)];
            end
        end
    else
        sz = size(symex);
        (lhs,rhs) = split_left_right(sterms,sz,var);
        if length(rhs[1]) == 0
            rhs[1] = [Basic(0)];
        end
    end
    
    # If there was a time derivative, separate those terms as well
    if timederiv
        dtlhs = copy(lhs);
        dtrhs = copy(rhs);
        if varcount > 1
            for i=1:length(symex)
                sz = size(symex[i]);
                (dtlhs[i],lhs[i]) = split_dt(lhs[i],sz);
                (dtrhs[i],rhs[i]) = split_dt(rhs[i],sz);
            end
        else
            sz = size(symex);
            (dtlhs,lhs) = split_dt(lhs,sz);
            (dtrhs,rhs) = split_dt(rhs,sz);
        end
    end
    
    # If there was a surface integral, separate those terms as well
    if has_surface
        surflhs = copy(lhs);
        surfrhs = copy(rhs);
        if varcount > 1
            for i=1:length(symex)
                sz = size(symex[i]);
                (surflhs[i],lhs[i]) = split_surf(lhs[i],sz);
                (surfrhs[i],rhs[i]) = split_surf(rhs[i],sz);
            end
        else
            sz = size(symex);
            (surflhs,lhs) = split_surf(lhs,sz);
            (surfrhs,rhs) = split_surf(rhs,sz);
        end
    end
    
    if debug println("volLHS = "*string(lhs)); end
    if debug println("volRHS = "*string(rhs)); end
    if timederiv
        if debug println("dtLHS = "*string(dtlhs)); end
        if debug println("dtRHS = "*string(dtrhs)); end
        if has_surface
            if debug println("surfLHS = "*string(surflhs)); end
            if debug println("surfRHS = "*string(surfrhs)); end
            return ((dtlhs,lhs), rhs, surflhs, surfrhs);
        else
            return ((dtlhs,lhs), rhs);
        end
        
    else
        if has_surface
            if debug println("surfLHS = "*string(surflhs)); end
            if debug println("surfRHS = "*string(surfrhs)); end
            return (lhs, rhs, surflhs, surfrhs);
        else
            return (lhs, rhs);
        end
    end
end

# Replaces parameter symbols with their values
function insert_parameters(ex)
    if typeof(ex) == Symbol
        # parameter?
        for p in parameters
            if ex === p.symbol
                if p.type == SCALAR
                    return insert_parameters(p.value[1]);
                else
                    return insert_parameters(p.value);
                end
            end
        end
        # nope
        return ex;
        
    elseif typeof(ex) == Expr && length(ex.args) > 0
        for i=1:length(ex.args)
            ex.args[i] = insert_parameters(ex.args[i]); # recursively replace on args if ex
        end
        return ex;
    elseif typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = insert_parameters(ex[i]);
        end
        return result;
    else
        return ex;
    end
end

# Replaces variable, coefficient and operator symbols in the expression
function replace_symbols(ex)
    if typeof(ex) == Symbol
        # variable?
        for v in variables
            if ex === v.symbol
                return v.symvar.vals;
            end
        end
        # coefficient?
        for c in coefficients
            if ex === c.symbol
                return c.symvar.vals;
            end
        end
        # operator?
        for i=1:length(ops)
            if ex === ops[i].symbol
                #return Symbol(string(ops[i].op));
                return :(ops[$i].op);
            end
        end
        for i=1:length(custom_ops)
            if ex === custom_ops[i].symbol
                #return Symbol(string(custom_ops[i].op));
                return :(custom_ops[$i].op);
            end
        end
        # test function?
        for c in test_functions
            if ex === c.symbol
                return c.symvar.vals;
            end
        end
        # none of them?
        return ex;
    elseif typeof(ex) == Expr && length(ex.args) > 0
        for i=1:length(ex.args)
            ex.args[i] = replace_symbols(ex.args[i]); # recursively replace on args if ex
        end
        return ex;
    elseif typeof(ex) <:Array
        result = copy(ex);
        for i=1:length(ex)
            result[i] = replace_symbols(ex[i]);
        end
        return result;
    else
        return ex;
    end
end

# Eval to apply the sym_*_op ops to create a SymEngine expression
function apply_ops(ex)
    try
        if typeof(ex) <:Array
            result = [];
            for i=1:length(ex)
                push!(result, eval(ex[i]));
            end
            return result;
        else
            return eval(ex);
        end
    catch e
        printerr("Problem evaluating the symbolic expression: "*string(ex));
        println(string(e));
        return 0;
    end
end

function has_unknown(ex, var)
    str = string(ex);
    result = false;
    if typeof(var) <: Array
        for i=1:length(var)
            vs = "_"*string(var[i])*"_";
            if occursin(vs, str)
                result = true;
            end
        end
    else
        vs = "_"*string(var)*"_";
        result = occursin(vs, str);
    end
    
    return result;
end

# Separate terms for each element of ex.
# Elements of the returned array are arrays of terms
function get_sym_terms(ex)
    # Recursively work on each term of the array
    if typeof(ex) <: Array
        sz = size(ex);
        result = Array{Array,length(sz)}(undef, sz);
        
        if length(sz) == 1 # vector or scalar
            for i=1:sz[1]
                result[i] = get_sym_terms(ex[i]);
            end
        elseif length(sz) == 2 # matrix
            for j=1:sz[2]
                for i=1:sz[1]
                    result[i,j] = get_sym_terms(ex[i,j]);
                end
            end
        elseif length(sz) == 3 # rank 3
            for k=1:sz[3]
                for j=1:sz[2]
                    for i=1:sz[1]
                        result[i,j,k] = get_sym_terms(ex[i,j,k]);
                    end
                end
            end
        end
        
        return result;
    end
    
    # ex is a symbolic expression(not array of them)
    # First expand it
    newex = expand(ex);
    #println("expanded = "*string(newex));
    # Then separate the terms into an array
    return get_all_terms(newex);
end

function get_all_terms(ex)
    if typeof(ex) == Basic
        # convert to Expr, separate, convert to Basic
        expr = Meta.parse(string(ex));
        #println("Expr = "*string(expr));
        #dump(expr);
        terms = get_all_terms(expr);
        #println("Exprterms = "*string(terms));
        bterms = Array{Basic,1}(undef,0);
        for i=1:length(terms)
            if !(terms[i] == 0)
                push!(bterms, Basic(terms[i]));
            end
        end
        return bterms;
    end
    
    # At this point ex must be an Expr or symbol or number
    terms = [];
    if !(typeof(ex) == Expr)
        push!(terms, ex);
        return terms;
    end
    if ex.head === :call
        if ex.args[1] === :+
            for i=2:length(ex.args)
                terms = append!(terms, get_all_terms(ex.args[i]));
            end
        elseif ex.args[1] === :-
            # apply -() to minused terms
            # Remember that the result of this will only be added terms, no minuses
            terms = get_all_terms(ex.args[2]);
            if length(ex.args) < 3
                for j=1:length(terms)
                    terms[j] = apply_negative(terms[j]);
                end
            else
                for i=3:length(ex.args)
                    tmp = get_all_terms(ex.args[i]);
                    for j=1:length(tmp)
                        tmp[j] = apply_negative(tmp[j]);
                    end
                    append!(terms, tmp);
                end
            end
        else
            push!(terms, ex);
        end
    else
        push!(terms, ex);
    end
    return terms;
end

function check_for_dt(terms)
    result = false;
    if typeof(terms) <: Array
        for i=1:length(terms)
            result = result || check_for_dt(terms[i]);
        end
    else
        result = occursin("TIMEDERIV", string(terms));
    end
    return result;
end

function check_for_surface(terms)
    result = false;
    if typeof(terms) <: Array
        for i=1:length(terms)
            result = result || check_for_surface(terms[i]);
        end
    else
        result = occursin("SURFACEINTEGRAL", string(terms));
    end
    return result;
end

function split_left_right(sterms,sz,var)
    lhs = copy(sterms); # set up the container right
    rhs = copy(sterms);
    if length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            lhs[i] = Array{Basic,1}(undef,0);
            rhs[i] = Array{Basic,1}(undef,0);
            for ti=1:length(sterms[i])
                if has_unknown(sterms[i][ti], var)
                    #println("lhs: "*string(sterms[i][ti]));
                    push!(lhs[i], sterms[i][ti]);
                else
                    #println("rhs: "*string(sterms[i][ti]));
                    # switch sign to put on RHS
                    push!(rhs[i], -sterms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                lhs[i,j] = Basic(0);
                rhs[i,j] = Basic(0);
                for ti=1:length(sterms[i,j])
                    if has_unknown(sterms[i,j][ti], var)
                        #println("lhs: "*string(sterms[i,j][ti]));
                        push!(lhs[i,j], sterms[i,j][ti]);
                    else
                        #println("rhs: "*string(sterms[i,j][ti]));
                        push!(rhs[i,j], -sterms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    lhs[i,j,k] = Basic(0);
                    rhs[i,j,k] = Basic(0);
                    for ti=1:length(sterms[i,j,k])
                        if has_unknown(sterms[i,j,k][ti], var)
                            #println("lhs: "*string(sterms[i,j,k][ti]));
                            push!(lhs[i,j,k], sterms[i,j,k][ti]);
                        else
                            #println("rhs: "*string(sterms[i,j,k][ti]));
                            push!(rhs[i,j,k], -sterms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    return (lhs,rhs);
end

function split_dt(terms, sz)
    hasdt = copy(terms); # set up the container right
    nodt = copy(terms);
    TIMEDERIV = symbols("TIMEDERIV"); # will be removed from terms
    if length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            hasdt[i] = Array{Basic,1}(undef,0);
            nodt[i] = Array{Basic,1}(undef,0);
            for ti=1:length(terms[i])
                if check_for_dt(terms[i][ti])
                    #println("hasdt: "*string(terms[i][ti]));
                    terms[i][ti] = subs(terms[i][ti], TIMEDERIV=>1);
                    push!(hasdt[i], terms[i][ti]);
                else
                    #println("nodt: "*string(terms[i][ti]));
                    push!(nodt[i], terms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                hasdt[i,j] = Basic(0);
                nodt[i,j] = Basic(0);
                for ti=1:length(terms[i,j])
                    if check_for_dt(terms[i,j][ti])
                        #println("hasdt: "*string(terms[i,j][ti]));
                        terms[i,j][ti] = subs(terms[i,j][ti], TIMEDERIV=>1);
                        push!(hasdt[i,j], terms[i,j][ti]);
                    else
                        #println("nodt: "*string(terms[i,j][ti]));
                        push!(nodt[i,j], terms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    hasdt[i,j,k] = Basic(0);
                    nodt[i,j,k] = Basic(0);
                    for ti=1:length(terms[i,j,k])
                        if check_for_dt(terms[i,j,k][ti])
                            #println("hasdt: "*string(terms[i,j,k][ti]));
                            terms[i,j,k][ti] = subs(terms[i,j,k][ti], TIMEDERIV=>1);
                            push!(hasdt[i,j,k], terms[i,j,k][ti]);
                        else
                            #println("nodt: "*string(terms[i,j,k][ti]));
                            push!(nodt[i,j,k], terms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    
    return (hasdt, nodt);
end

function split_surf(terms, sz)
    hassurf = copy(terms); # set up the container right
    nosurf = copy(terms);
    SURFACEINTEGRAL = symbols("SURFACEINTEGRAL"); # will be removed from terms
    if length(sz) == 1 # vector or scalar
        for i=1:sz[1]
            hassurf[i] = Array{Basic,1}(undef,0);
            nosurf[i] = Array{Basic,1}(undef,0);
            for ti=1:length(terms[i])
                if check_for_surface(terms[i][ti])
                    terms[i][ti] = subs(terms[i][ti], SURFACEINTEGRAL=>1);
                    push!(hassurf[i], terms[i][ti]);
                else
                    push!(nosurf[i], terms[i][ti]);
                end
            end
        end
    elseif length(sz) == 2 # matrix
        for j=1:sz[2]
            for i=1:sz[1]
                hassurf[i,j] = Basic(0);
                nosurf[i,j] = Basic(0);
                for ti=1:length(terms[i,j])
                    if check_for_surface(terms[i,j][ti])
                        terms[i,j][ti] = subs(terms[i,j][ti], SURFACEINTEGRAL=>1);
                        push!(hassurf[i,j], terms[i,j][ti]);
                    else
                        push!(nosurf[i,j], terms[i,j][ti]);
                    end
                end
            end
        end
    elseif length(sz) == 3 # rank 3
        for k=1:sz[3]
            for j=1:sz[2]
                for i=1:sz[1]
                    hassurf[i,j,k] = Basic(0);
                    nosurf[i,j,k] = Basic(0);
                    for ti=1:length(terms[i,j,k])
                        if check_for_surface(terms[i,j,k][ti])
                            terms[i,j,k][ti] = subs(terms[i,j,k][ti], SURFACEINTEGRAL=>1);
                            push!(hassurf[i,j,k], terms[i,j,k][ti]);
                        else
                            push!(nosurf[i,j,k], terms[i,j,k][ti]);
                        end
                    end
                end
            end
        end
    end
    
    return (hassurf, nosurf);
end

function apply_negative(ex)
    negex = :(-a);
    if typeof(ex) == Symbol
        negex.args[2] = ex;
        return negex;
    elseif typeof(ex) <: Number
        return -ex;
    elseif typeof(ex) == Expr
        if ex.args[1] == :- && length(ex.args) == 2
            return ex.args[2];
        else
            negex.args[2] = ex;
            return negex;
        end
    end
end

end # module