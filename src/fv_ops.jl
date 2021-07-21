# The finite volume operators that are included automatically if FV is used.

function sym_upwind_op(v, a, alpha=Basic(0))
    # Input will be an array of Basic. Output should be in a similar array.
    if typeof(a) <: Array
        result = copy(a);
        for i=1:length(result)
            result[i] = sym_upwind_op(v, a[i], alpha);
        end
        
    elseif typeof(a) == Basic
        # Apply a DGSIDEn tag to signal to the code generator which side of the face the value is from.
        # TODO: add alias like CELL1/CELL2 more appropriate for FV
        side1 = symbols("DGSIDE1_"*string(a));
        side2 = symbols("DGSIDE2_"*string(a));
        # The formula for an adjustable first order upwind scheme.
        # F(u) = 0.5*(side1+side2)*(v.normal) + 0.5*(side1-side2)*abs(v.normal)*(1-alpha)
        result = Basic(0.5) .* ( (sym_dot_op(v, sym_normal_op())) .* (side1 + side2) .+ abs.(sym_dot_op(v, sym_normal_op())[1]) .* (Basic(1.0)-alpha) .* (side1 - side2) );
        result = result[1]; # Since dot and normal returned arrays, result was placed in an array, but this is done above redundantly.
        
    elseif typeof(a) <: Number
        # If the input was just a constant, the result will just be that constant.
        result = Basic(a);
    end
    return result;
end

_names = [:upwind];
_handles = [sym_upwind_op];
