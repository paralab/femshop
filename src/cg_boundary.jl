# Apply boundary conditions to the system

# This is just to speed things up (huge improvement)
function is_in_rows(k, rows)
    # rows are sorted
    b = 1;
    e = length(rows);
    c = Int(round((e+1)/2));
    while e-b > 1
        if k == rows[c]
            return true;
        elseif k > rows[c]
            b=c;
        else
            e=c;
        end
        c = Int(round((b+e)/2));
    end
    return k==rows[b] || k==rows[e];
end

# zeros the rows and puts a 1 on the diagonal
function identity_rows(A, rows, N)
    if issparse(A)
        (I, J, V) = findnz(A);
        eN = length(I);
        newI = zeros(Int,eN);
        newJ = zeros(Int,eN);
        newV = zeros(eN);
        newind = 0;
        (M,N) = size(A);
        
        sort!(rows);
        # Remove bdry row elements
        for k in 1:length(I)
            isbdry = is_in_rows(I[k], rows);
            if !isbdry
                newind = newind+1;
                newI[newind] = I[k];
                newJ[newind] = J[k];
                newV[newind] = V[k];
            end
        end
        
        # Add diogonal 1s
        rn = length(rows);
        if newind + rn <= eN
            for i=1:rn
                newind = newind+1;
                newI[newind] = rows[i];
                newJ[newind] = rows[i];
                newV[newind] = 1.0;
            end
            newI = newI[1:newind];
            newJ = newJ[1:newind];
            newV = newV[1:newind];
        else
            toadd = rn - (eN-newind);
            for i=1:(eN-newind)
                newind = newind+1;
                newI[newind] = rows[i];
                newJ[newind] = rows[i];
                newV[newind] = 1.0;
            end
            append!(newI, rows[toadd:rn]);
            append!(newJ, rows[toadd:rn]);
            append!(newV, ones(toadd));
        end
        
        return sparse(newI,newJ,newV, M,N);
        
    else
        for i=1:length(rows)
            A[rows[i],:] = zeros(N);
            A[rows[i],rows[i]] = 1;
        end
        return A;
    end
end

# Inserts the rows of S in A
function insert_rows(A, S, rows, N)
    if issparse(A)
        (I, J, V) = findnz(A);
        (M,N) = size(A);
        # determine which elements to remove
        toskip = zeros(Int, length(I));
        includen = 0;
        for k = 1:length(I)
            for i=1:length(rows)
                if I[k] == rows[i]
                    toskip[k] = 1;
                end
            end
            if toskip[k] == 0
                includen = includen + 1;
            end
        end
        
        # put the remaining elements in a new matrix
        newI = zeros(Int, includen);
        newJ = zeros(Int, includen);
        newV = zeros(includen);
        ind = 1;
        for k = 1:length(I)
            if toskip[k] == 0
                newI[ind] = I[k];
                newJ[ind] = J[k];
                newV[ind] = V[k];
                ind = ind + 1;
            end
        end
        
        # Do something similar to extract elements of S
        (I, J, V) = findnz(S);
        
        toskip = ones(Int, length(I));
        includen = 0;
        for k = 1:length(I)
            for i=1:length(rows)
                if I[k] == rows[i]
                    toskip[k] = 0;
                end
            end
            if toskip[k] == 0
                includen = includen + 1;
            end
        end
        newI2 = zeros(Int, includen);
        newJ2 = zeros(Int, includen);
        newV2 = zeros(includen);
        ind = 1;
        for k = 1:length(I)
            if toskip[k] == 0
                newI2[ind] = I[k];
                newJ2[ind] = J[k];
                newV2[ind] = V[k];
                ind = ind + 1;
            end
        end
        
        append!(newI, newI2);
        append!(newJ, newJ2);
        append!(newV, newV2);
        return sparse(newI,newJ,newV, M, N);
    else
        for i=1:length(rows)
            A[rows[i],:] = S[rows[i],:];
        end
        return A;
    end
end

function insert_sparse_rows(A, SI, SJ, SV)
    (I, J, V) = findnz(A);
    (M,N) = size(A);
    
    # Zero existing elements in SI rows
    for k in 1:length(I)
        if is_in_rows(I[k], SI)
            V[k] = 0;
        end
    end
    
    # append S values
    append!(I, SI);
    append!(J, SJ);
    append!(V, SV);
    
    # Form sparse matrix
    return sparse(I, J, V, M, N);
end

function dirichlet_bc(A, b, val, bdryind, t=0, dofind = 1, totaldofs = 1)
    N = length(b);
    bdry = copy(bdryind);
    if totaldofs > 1
        bdry = (bdry .- 1) .* totaldofs .+ dofind;
    end
    #A = identity_rows(A, bdry, N);
    
    b = dirichlet_bc_rhs_only(b, val, bdryind, t, dofind, totaldofs); # how convenient
    
    #return (A, b);
    return (bdry, b);
end

function dirichlet_bc_rhs_only(b, val, bdryind, t=0, dofind = 1, totaldofs = 1)
    N = length(b);
    bdry = copy(bdryind);
    vecbdry = copy(bdryind);
    if totaldofs > 1
        vecbdry = (vecbdry .- 1) .* totaldofs .+ dofind;
    end
    
    if typeof(val) <: Number
        for i=1:length(bdry)
            b[vecbdry[i]]=val;
        end
        
    elseif typeof(val) == Coefficient && typeof(val.value[1]) == GenFunction
        if config.dimension == 1
            for i=1:length(bdry)
                b[vecbdry[i]]=val.value[1].func(grid_data.allnodes[1,bdry[i]],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bdry)
                b[vecbdry[i]]=val.value[1].func(grid_data.allnodes[1,bdry[i]],grid_data.allnodes[2,bdry[i]],0,t);
            end
        else
            for i=1:length(bdry)
                b[vecbdry[i]]=val.value[1].func(grid_data.allnodes[1,bdry[i]],grid_data.allnodes[2,bdry[i]],grid_data.allnodes[3,bdry[i]],t);
            end
        end
        
    elseif typeof(val) == GenFunction
        if config.dimension == 1
            for i=1:length(bdry)
                b[vecbdry[i]]=val.func(grid_data.allnodes[1,bdry[i]],0,0,t);
            end
        elseif config.dimension == 2
            for i=1:length(bdry)
                b[vecbdry[i]]=val.func(grid_data.allnodes[1,bdry[i]],grid_data.allnodes[2,bdry[i]],0,t);
            end
        else
            for i=1:length(bdry)
                b[vecbdry[i]]=val.func(grid_data.allnodes[1,bdry[i]],grid_data.allnodes[2,bdry[i]],grid_data.allnodes[3,bdry[i]],t);
            end
        end
    end
    
    return b;
end

function neumann_bc(A, b, val, bdryind, bid, t=0, dofind = 1, totaldofs = 1)
    N = length(b);
    bdry = copy(bdryind);
    if totaldofs > 1
        bdry = (bdry .- 1) .* totaldofs .+ dofind;
    end
    
    # Assemble the differentiation matrix and swap rows in A
    # (na, nb) = size(A);
    # S1 = spzeros(na, nb);
    # S2 = spzeros(na, nb);
    # S3 = spzeros(na, nb);
    # S = spzeros(na, nb);
    S_I = zeros(Int, 0);
    S_J = zeros(Int, 0);
    S1_V = zeros(0);
    S2_V = zeros(0);
    S3_V = zeros(0);
    S_V = zeros(0);
    
    bfc = grid_data.bdryface[bid];
    
    # This can be precomputed for uniform grid meshes
    if config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
        precomputed = true;
        glb = grid_data.loc2glb[:,1];
        xe = grid_data.allnodes[:,glb[:]];
        (detJ, J) = geometric_factors(refel, xe);
        if config.dimension == 1
            R1D1 = diagm(J.rx) * refel.Dr; #R1matrix * D1matrix;
            
        elseif config.dimension == 2
            R1D1 = [diagm(J.rx) diagm(J.sx)] * [refel.Ddr ; refel.Dds]; #R1matrix * D1matrix;
            R2D2 = [diagm(J.ry) diagm(J.sy)] * [refel.Ddr ; refel.Dds]; #R2matrix * D2matrix;
            
        elseif config.dimension == 3
            R1D1 = [diagm(J.rx) diagm(J.sx) diagm(J.tx)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R1matrix * D1matrix;
            R2D2 = [diagm(J.ry) diagm(J.sy) diagm(J.ty)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R2matrix * D2matrix;
            R3D3 = [diagm(J.rz) diagm(J.sz) diagm(J.tz)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R3matrix * D3matrix;
        end
    else
        precomputed = false;
    end
    
    for fi = 1:length(bfc)
        # find the relevant element
        e = grid_data.face2element[1,bfc[fi]] # >0 ? grid_data.face2element[1,bfc[fi]] : grid_data.face2element[2,bfc[fi]];
        glb = grid_data.loc2glb[:,e];       # global indices of this element's nodes
        xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes
        
        # Local indices for the nodes on the face
        #flocal = get_face_local_index(grid_data.face2glb[:,1,bfc[fi]], glb);
        flocal = refel.face2local[grid_data.faceRefelInd[1,bfc[fi]]];
        
        # offset for multi dof
        if totaldofs > 1
            sglb = (glb.-1) .* totaldofs .+ dofind;
        else
            sglb = glb;
        end
        
        # Build diff matrices if needed
        if !precomputed
            (detJ, J) = geometric_factors(refel, xe);
            if config.dimension == 1
                R1D1 = diagm(J.rx) * refel.Dr; #R1matrix * D1matrix;
                
            elseif config.dimension == 2
                R1D1 = [diagm(J.rx) diagm(J.sx)] * [refel.Ddr ; refel.Dds]; #R1matrix * D1matrix;
                R2D2 = [diagm(J.ry) diagm(J.sy)] * [refel.Ddr ; refel.Dds]; #R2matrix * D2matrix;
                
            elseif config.dimension == 3
                R1D1 = [diagm(J.rx) diagm(J.sx) diagm(J.tx)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R1matrix * D1matrix;
                R2D2 = [diagm(J.ry) diagm(J.sy) diagm(J.ty)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R2matrix * D2matrix;
                R3D3 = [diagm(J.rz) diagm(J.sz) diagm(J.tz)] * [refel.Ddr ; refel.Dds ; refel.Ddt]; #R3matrix * D3matrix;
            end
        end
        
        # Add to S matrix
        Imat = zeros(Int,length(flocal),length(sglb));
        Jmat = zeros(Int,size(Imat));
        for i=1:length(sglb)
            Imat[:,i] = sglb[flocal];
            Jmat[:,i] = fill(sglb[i], length(flocal));
        end
        append!(S_I, Imat[:]);
        append!(S_J, Jmat[:]);
        if config.dimension == 1
            append!(S1_V, R1D1[flocal,:][:]);
            
        elseif config.dimension == 2
            append!(S1_V, R1D1[flocal,:][:]);
            append!(S2_V, R2D2[flocal,:][:]);
            
        elseif config.dimension == 3
            append!(S1_V, R1D1[flocal,:][:]);
            append!(S2_V, R2D2[flocal,:][:]);
            append!(S3_V, R3D3[flocal,:][:]);
        end
    end
    
    # Add the right components of S1,S2,S3 according to normal vector
    #S_V = zeros(length(S1_V));
    for i=1:length(bdry)
        norm = grid_data.bdrynorm[bid][:,i];
        if config.dimension == 1
            S_V = norm[1] .* S1_V;
        elseif config.dimension == 2
            S_V = norm[1] .* S1_V + norm[2] .* S2_V;
        elseif config.dimension == 3
            S_V = norm[1] .* S1_V + norm[2] .* S2_V + norm[3] .* S3_V;
        end
    end
    
    #A = insert_rows(A, S, bdry, N);
    
    b = dirichlet_bc_rhs_only(b, val, bdryind, t, dofind, totaldofs); # how convenient
    
    #return (A, b);
    return (bdry, S_I, S_J, S_V, b);
end

function neumann_bc_rhs_only(b, val, bdryind, bid, t=0, dofind = 1, totaldofs = 1)
    return dirichlet_bc_rhs_only(b, val, bdryind, t, dofind, totaldofs); # how convenient
end

# # Returns the local indices for a face given the face and element global indices
# function get_face_local_index(f2glb, e2glb)
#     ind = zeros(Int, length(f2glb));
#     for fi=1:length(f2glb)
#         for ei=1:length(e2glb)
#             if f2glb[fi] == e2glb[ei]
#                 ind[fi] = ei;
#                 break;
#             end
#         end
#     end
    
#     return ind;
# end