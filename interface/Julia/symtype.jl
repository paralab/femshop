#=
Types: scalar, vector, tensor(r=2), symmetric tensor(r=2)
Higher rank tensors not yet ready.
=#

# An abstract value type for variables, coefficients, etc.
# Use sym_var() to build.
struct SymType
    rank::Int   # tensor rank
    dim::Int    # dimension
    vals::Array # Basic symbols for values (scalar->[u_1], vector->[u_1, u_2,...], matrix->[u_11 u_12 ...])
end

# Builds a SymType
function sym_var(name, type, dim)
    symvar = [];
    rank = 0;
    if type == SCALAR
        symvar = [symbols("_"*name*"_1")];
    elseif type == VECTOR
        symvar = [symbols("_"*name*"_$i") for i=1:dim];
        rank = 1;
    elseif type == TENSOR
        symvar = [symbols("_"*name*"_$i$j") for i=1:dim, j=1:dim];
        rank = 2;
    elseif type == SYM_TENSOR
        if dim == 1
            symvar = [symbols("_"*name*"_1")];
        elseif dim == 2
            symvar = [symbols("_"*name*"_11") symbols("_"*name*"_12"); symbols("_"*name*"_12") symbols("_"*name*"_22")];
        elseif dim == 3
            symvar = [symbols("_"*name*"_11") symbols("_"*name*"_12") symbols("_"*name*"_13"); symbols("_"*name*"_12") symbols("_"*name*"_22") symbols("_"*name*"_23"); symbols("_"*name*"_13") symbols("_"*name*"_23") symbols("_"*name*"_33")];
        elseif dim == 4
            symvar = [symbols("_"*name*"_11") symbols("_"*name*"_12") symbols("_"*name*"_13") symbols("_"*name*"_14"); symbols("_"*name*"_12") symbols("_"*name*"_22") symbols("_"*name*"_23") symbols("_"*name*"_24"); symbols("_"*name*"_13") symbols("_"*name*"_23") symbols("_"*name*"_33") symbols("_"*name*"_34"); symbols("_"*name*"_14") symbols("_"*name*"_24") symbols("_"*name*"_34") symbols("_"*name*"_44")];
        end
        rank = 2;
    else
        # unknown type
        printerr("unknown type for: "*name);
    end
    
    sym_type = SymType(rank, dim, symvar);
    
    return sym_type;
end

# Applies a derivative prefix. wrt is the axis index
# sym_deriv(u_12, 1) -> D1_u_12
function sym_deriv(var, wrt)
    prefix = "D"*string(wrt)*"_";
    
    return symbols(prefix*string(var));
end
