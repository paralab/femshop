function FV_flux_bc_rhs_only(val, facex, Qvec, t=0, dofind = 1, totaldofs = 1)
    if typeof(val) <: Number
        return val;
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        bvals = zeros(size(facex,2));
        if config.dimension == 1
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],facex[2,i],0,t);
            end
        else
            for i=1:length(bvals)
                bvals[i]=val.value[1].func(facex[1,i],facex[2,i],facex[3,i],t);
            end
        end
        
    elseif typeof(val) == GenFunction
        bvals = zeros(size(facex,2));
        if config.dimension == 1
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],facex[2,i],0,t);
            end
        else
            for i=1:length(bvals)
                bvals[i]=val.func(facex[1,i],facex[2,i],facex[3,i],t);
            end
        end
    end
    
    # Do quadrature over the face
    b = Qvec * bvals;
    
    return b[1];
end

function FV_flux_bc_rhs_only_simple(val, facex, t=0)
    if typeof(val) <: Number
        return val;
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        if config.dimension == 1
            bval = val.value[1].func(facex[1],0,0,t);
        elseif config.dimension == 2
            bval = val.value[1].func(facex[1],facex[2],0,t);
        else
            bval = val.value[1].func(facex[1],facex[2],facex[3],t);
        end
        
    elseif typeof(val) == GenFunction
        if config.dimension == 1
            bval = val.func(facex[1],0,0,t);
        elseif config.dimension == 2
            bval = val.func(facex[1],facex[2],0,t);
        else
            bval = val.func(facex[1],facex[2],facex[3],t);
        end
    end
    
    return bval;
end