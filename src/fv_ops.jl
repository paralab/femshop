# The finite volume operators that are included automatically if FV is used.

function sym_upwind_op(v, a, alpha=[Basic(0)])
    side1 = "DGSIDE1_";
    side2 = "DGSIDE2_";
    if typeof(a) <: Array
        result = copy(a);
        for i=1:length(result)
            result[i] = sym_upwind_op(v, a[i], alpha);
        end
    elseif typeof(a) == Basic
        side1 = symbols(side1*string(a));
        side2 = symbols(side2*string(a));
        result = Basic(0.5) .* ( (sym_dot_op(v, sym_normal_op())) .* (side1 + side2) .+ abs.(sym_dot_op(v, sym_normal_op())) .* ([Basic(1.0)]-alpha) .* (side1 - side2) );
        result = result[1];
    elseif typeof(a) <: Number
        result = Basic(a);
    end
    return result;
end

_names = [:upwind];
_handles = [sym_upwind_op];
