# Apply boundary conditions to the system

# zeros the row and puts a 1 on the diagonal
function identity_rows(A, rows, N)
    if issparse(A)
        (I, J, V) = findnz(A);
        topush = [];
        for i=1:length(rows)
            ind = rows[i];
            gotit = false;
            for k in 1:length(I)
                if I[k] == ind
                    if J[k] == ind
                        V[k] = 1;
                        gotit = true;
                    else
                        V[k] = 0;
                    end
                end
            end
            if !gotit
                push!(topush, ind);
            end
        end
        append!(I, topush);
        append!(J, topush);
        append!(V, ones(length(topush)));
        return sparse(I,J,V);
    else
        for i=1:length(rows)
            A[rows[i],:] = zeros(N);
            A[rows[i],rows[i]] = 1;
        end
        return A;
    end
end

function dirichlet_bc(A, b, val, bdry, t=0)
    N = length(b);
    A = identity_rows(A, bdry, N);
    
    if config.dimension == 1
        for i=1:length(bdry)
            if typeof(val) <: Number
                b[bdry[i]]=val;
            elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
                b[bdry[i]]=val.value[1].func(grid_data.allnodes[i,1],0,0,t);
            end
        end
        
    elseif config.dimension == 2
        for i=1:length(bdry)
            if typeof(val) <: Number
                b[bdry[i]]=val;
            elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
                b[bdry[i]]=val.value[1].func(grid_data.allnodes[i,1],grid_data.allnodes[i,2],0,t);
            end
        end
        
    else
        
    end
    
    return (A, b);
end

function dirichlet_bc_rhs_only(b, val, bdry, t=0)
    N = length(b);
    
    if config.dimension == 1
        for i=1:length(bdry)
            if typeof(val) <: Number
                b[bdry[i]]=val;
            elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
                b[bdry[i]]=val.value[1].func(grid_data.allnodes[i,1],0,0,t);
            end
        end
        
    elseif config.dimension == 2
        for i=1:length(bdry)
            if typeof(val) <: Number
                b[bdry[i]]=val;
            elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
                b[bdry[i]]=val.value[1].func(grid_data.allnodes[i,1],grid_data.allnodes[i,2],0,t);
            end
        end
        
    else
        
    end
    
    return b;
end