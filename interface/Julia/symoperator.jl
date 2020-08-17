#=
Operators that work on SymType objects
=#

# An abstract operator type.
# Use sym_op() to build.
struct SymOperator
    input_count::Int    # The number of input arguments
    input_rank::Array   # The possible input SymType ranks (an array of tuples for possible rank combos)
    op                  # Function handle for the operator
end

function sym_dot_op(a,b)
    if (typeof(a) == Array{Basic,1} || typeof(a) == Array{Float64,1}) && (typeof(b) == Array{Basic,1} || typeof(b) == Array{Float64,1})
        return [transpose(a)*b];
    else
        printerr("Unsupported operation: dot("*string(typeof(a))*", "*string(typeof(b))*")");
    end
end

function sym_grad_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) <: Array
        d = config.dimension;
        sz = size(u);
        rank = 0;
        if length(sz) == 1 && sz[1] > 1
            rank = 1;
        elseif length(sz) == 2 
            rank = 2;
        end
        
        if rank == 0
            # result is a vector
            for i=1:d
                push!(result, sym_deriv(u[1], i));
            end
        elseif rank == 1
            # result is a tensor
            for i=1:d
                for j=1:d
                    push!(result, sym_deriv(u[i], j));
                end
            end
            reshape(result, d,d);
        elseif rank == 2
            # not yet ready
            printerr("unsupported operator, grad(tensor)");
            return nothing;
        end
    elseif typeof(u) == Basic
        # result is a vector
        d = config.dimension;
        for i=1:d
            push!(result, sym_deriv(u, i));
        end
    elseif typeof(u) <: Number
        return zeros(config.dimension);
    end
    
    # wrap result in a SymType
    return result;
end

function sym_div_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) == SymType
        d = config.dimension;
        if u.rank == 0
            # Not allowed
            printerr("unsupported operator, div(scalar)");
            return nothing;
        elseif u.rank == 1
            # result is a scalar
            if d==1
                result = [sym_deriv(u.vals[1], 1)];
            else
                ex = :(a+b);
                ex.args = [:+];
                for i=1:d
                    push!(ex.args, sym_deriv(u.vals[i], i))
                end
                result = [Basic(ex)];
            end
        elseif u.rank == 2
            # not yet ready
            printerr("unsupported operator, div(tensor)");
            return nothing;
        end
    elseif typeof(u) <: Number
        # Not allowed
        printerr("unsupported operator, div(number)");
        return nothing;
    end
    
    # wrap result in a SymType
    return SymType(u.rank-1, u.dim, result);
end

function sym_curl_op(u)
    result = Array{Basic,1}(undef,0);
    if typeof(u) == SymType
        d = config.dimension;
        if u.rank == 0
            # Not allowed
            printerr("unsupported operator, curl(scalar)");
            return nothing;
        elseif u.rank == 1
            # result is a vector
            if d==1
                result = [sym_deriv(u.vals[1], 1)];
            else
                #TODO
                printerr("curl not ready");
                return nothing;
            end
        elseif u.rank == 2
            # not yet ready
            printerr("unsupported operator, curl(tensor)");
            return nothing;
        end
    elseif typeof(u) <: Number
        # Not allowed
        printerr("unsupported operator, curl(number)");
        return nothing;
    end
    
    # wrap result in a SymType
    return SymType(u.rank, u.dim, result);
end

function sym_laplacian_op(u)
    # simply use the above ops
    return sym_div_op(sym_grad_op(u));
end