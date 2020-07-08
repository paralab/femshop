#=
# Construct the elemental linear terms and bilinear operators
# For each of these "args" contains:
#   1. variable being solved for
#   2. node coordinates for evaluating coefficient functions
#   3. global indices of nodes for extracting values from variable arrays
#   4. reference element
#   5. bilinear or linear? bilinear builds a matrix operator, linear builds a vector
#   6. time for time dependent coefficients
=#
export mass_operator, stiffness_operator

function zero_operator(ex, args)
    ref = args[4];
    return zeros(ref.Np);
end

function one_operator(ex,args)
    ref = args[4];
    return ones(ref.Np);
end

function mass_operator(ex, args)
    var = args[1];  # 
    x = args[2];    # global coords of element's nodes
    gbl = args[3];  # global indices of the nodes
    refel = args[4];# 
    borl = args[5]; # bilinear or linear? lhs or rhs?
    time = args[6]; # time for time dependent coefficients
    
    # get geometric factors for this element
    (detJ, J) = geometric_factors(refel, x);
    
    if borl == RHS # linear
        f = zeros(refel.Np); # function values
        # ex could be a coefficient, other variable, or...
        if typeof(ex) == Coefficient && typeof(ex.value[1]) == GenFunction
            if refel.dim == 1
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i],0,0,time);
                end
            elseif refel.dim == 2
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i,1],x[i,2],0,time);
                end
            elseif refel.dim == 3
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i,1],x[i,2],x[i,3],time);
                end
            end
        elseif typeof(ex) == Variable && ex.type == SCALAR
            f = ex.values[gbl];
        end
        
        result = elemental_mass(refel, detJ) * f;
        
    else # bilinear
        result = elemental_mass(refel, detJ);
    end
    
    return result;
end

function stiffness_operator(ex, args)
    var = args[1];  # 
    x = args[2];    # global coords of element's nodes
    gbl = args[3];  # global indices of the nodes
    refel = args[4];# 
    borl = args[5]; # bilinear or linear? lhs or rhs?
    time = args[6]; # time for time dependent coefficients
    
    (detJ, J) = geometric_factors(refel, x);
    
    if borl == RHS # linear
        f = zeros(refel.Np); # function values
        # ex could be a coefficient, other variable, or...
        if typeof(ex) == Coefficient && typeof(ex.value[1]) == GenFunction
            if refel.dim == 1
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i],0,0,time);
                end
            elseif refel.dim == 2
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i,1],x[i,2],0,time);
                end
            elseif refel.dim == 3
                for i=1:refel.Np
                    f[i] = ex.value[1].func(x[i,1],x[i,2],x[i,3],time);
                end
            end
        elseif typeof(ex) == Variable && ex.type == SCALAR
            f = ex.values[gbl];
        end
        
        result = elemental_stiffness(refel, detJ, J, x) * f;
        
    else # bilinear
        result = elemental_stiffness(refel, detJ, J, x);
    end
    
    return result;
end
