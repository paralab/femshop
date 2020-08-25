#= 
# A set of tools for parsing the variational forms into symEngine expressions.
=#
module SymbolicParser
export sp_parse, add_custom_op

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, NEUMANN, ROBIN,
        MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR, SYM_TENSOR,
        LHS, RHS,
        LINEMESH, QUADMESH, HEXMESH
import ..Femshop: Femshop_config, Femshop_prob, Variable, Coefficient
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, coefficients, test_functions

using SymEngine

include("symtype.jl");
include("symoperator.jl");

# Predefined operator symbols that will be replaced with sym_*_op 
op_names = [:dot, :inner, :cross, :grad, :div, :curl, :laplacian];
# Build SymOperator objects for those
ops = init_ops();

# Custom operators that can be added by a user
custom_op_names = [];
custom_ops = [];

function add_custom_op(s, handle)
    push!(custom_op_names, s); # add the Symbol s to the list
    push!(custom_ops, SymOperator(s, handle));
end

# a special operator for dealing with scalar multiplication when the scalar is in an array
import Base.*
function *(a::Array{Basic,1}, b::Array{Basic,1})
    if size(a) == (1,)
        return b.*a[1];
    elseif size(b) == (1,)
        return a.*b[1];
    else
        # This should be an error, but I will treat it as a dot product for lazy people
        # Note, it will still be an error if dimensions are not right.
        return [transpose(a) * b];
    end
end

# Parses a variational form expression into a SymEngine Basic expression 
# Takes an Expr, variable symbol, and test function symbol
# The symbols are only for determining lhs vs. rhs
# Returns a tuple of SymEngine expressions (lhs, rhs) for the equation lhs = rhs
# lhs contains terms including the unknown variable
# rhs contains terms without it
function sp_parse(ex, var)
    lhs = nothing;
    rhs = nothing;
    varcount = 1;
    #println("expr = "*string(ex));
    
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
    
    # Replace symbols for variables, coefficients, test functions, and special operators
    symex = replace_symbols(symex);
    #println("symex1 = "*string(symex));
    
    # Evaluate the expression to apply symbolic operators
    symex = apply_ops(symex);
    #println("symex2 = "*string(symex));
    
    # Expand the expression and separate terms
    sterms = get_sym_terms(symex);
    #println("sterms = "*string(sterms));
    
    # Each element has an lhs and rhs
    lhs = copy(sterms); # set up the container right
    rhs = copy(sterms);
    if varcount > 1
        for i=1:length(symex)
            sz = size(symex[i]);
            #println("sz="*string(sz));
            #println("ln="*string(length(sterms[i][1])))
            (lhs[i],rhs[i]) = split_left_right(sterms[i],sz,var);
            #println("lhs"*string(i)*" = "*string(lhs[i]));
            #println("rhs"*string(i)*" = "*string(rhs[i]));
        end
    else
        sz = size(symex);
        (lhs,rhs) = split_left_right(sterms,sz,var);
    end
    
    #println("LHS = "*string(lhs));
    #println("RHS = "*string(rhs));
    return (lhs, rhs);
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
        for p in ops
            if ex === p.symbol
                return Symbol("sym_"*string(ex)*"_op");
            end
        end
        for p in custom_ops
            if ex === p.symbol
                return Symbol("sym_"*string(ex)*"_op");
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
        result = [];
        for i=1:length(ex)
            push!(result, replace_symbols(ex[i]));
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
        bterms = Array{Basic,1}(undef,size(terms));
        for i=1:length(terms)
            bterms[i] = Basic(terms[i]);
        end
        return bterms;
    end
    
    # At this point ex must be an Expr or symbol or number
    if !(typeof(ex) == Expr)
        return [ex];
    end
    terms = [];
    if ex.head === :call
        if ex.args[1] === :+
            for i=2:length(ex.args)
                terms = [terms; get_all_terms(ex.args[i])];
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
            terms = [ex];
        end
    else
        terms = [ex];
    end
    return terms;
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