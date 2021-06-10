function FV_dirichlet_bc_rhs_only(val, facex, t=0, dofind = 1, totaldofs = 1)
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
    # TODO, for now just average
    b = sum(bvals)/length(bvals);
    
    return b;
end