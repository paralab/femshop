#=
# Contains info about all nodes on the domain
=#
struct DG_Grid
    allnodes::Array{Float64}        # All node coordinates
    # boundaries
    bdry::Array{Array{Int,1},1}     # Indices of boundary nodes for each BID (bdry[bid][nodes])
    bdryface::Array{Array{Int,1},1} # Indices of faces touching each BID 
    bdrynorm::Array{Array{Float64,2},1} # Normal vector for boundary nodes for each BID
    bids::Array{Int,1}              # BID corresponding to rows of bdrynodes
    # elements
    loc2glb::Array{Int,2}           # local to global map for each element's nodes
    glbvertex::Array{Int,2}         # global indices of each elements' vertices
    # faces
    face2glb::Array{Int,3}          # global map for faces
    faceVertex2glb::Array{Int,3}    # global indices of face vertices
    face2local::Array{Int,3}        # local map for faces
    facenorm::Array{Float64,2}      # Normal vector for each face
    faceRefelInd::Array{Int,2}      # Index for face matching the one used by refel for each side
end

function cg_grid_to_dg_grid(cggrid, mesh)
    dim = size(cggrid.allnodes, 1);
    nel = size(cggrid.loc2glb, 2);
    np = size(cggrid.loc2glb, 1);
    nf = size(mesh.element2face, 1);
    nface = size(cggrid.face2glb, 2);
    nfp = size(cggrid.face2glb, 1);
    nfv = size(cggrid.faceVertex2glb, 1);
    
    dgnodes = zeros(dim, nel*np);
    
    dgbids = copy(cggrid.bids);
    dgbdryface = copy(cggrid.bdryface);
    dgbdry = Array{Array{Int,1},1}(undef,length(dgbids));
    dgbdrynorm = Array{Array{Float64,2},1}(undef,length(dgbids));
    for i=1:length(dgbids)
        dgbdry[i] = zeros(Int, length(dgbdryface[i])*nfp);
        dgbdrynorm[i] = zeros(dim, length(dgbdryface[i])*nfp);
    end
    
    dgloc2glb = zeros(Int, size(cggrid.loc2glb));
    dgglbvertex = zeros(Int, size(cggrid.glbvertex));
    
    dgface2glb = zeros(Int, nfp, 2, nface);
    dgfaceVertex2glb = zeros(Int, size(cggrid.faceVertex2glb,1), 2, nface);
    dgface2local = zeros(Int, nfp, 2, nface);
    dgfacenorm = copy(mesh.normals);
    dgfaceRefelInd = zeros(Int, 2, nface);
    
    # Loop over elements.
    #   Add each element's nodes to dgnodes.
    #   Set loc2glb and glbvertex.
    # Loop over that element's faces.
    #   Set face2glb and faceVertex2glb.
    #   Set normals.
    # Loop over boundary faces.
    #   Set bdry and bdrynorm.
    
    # element loop
    for ei=1:nel
        offset = (ei-1)*np+1;
        dgnodes[:, offset:offset+np-1] = cggrid.allnodes[:, cggrid.loc2glb[:,ei]];
        dgloc2glb[:, ei] = offset:offset+np-1;
        vi = 1;
        fnodeind = ones(Int, nf);
        for ni=1:np 
            for nj=1:length(cggrid.glbvertex[:,ei])
                if cggrid.glbvertex[nj,ei] == cggrid.loc2glb[ni,ei]
                    dgglbvertex[vi,ei] = dgloc2glb[ni,ei];
                    vi = vi+1;
                    break;
                end
            end
        end
    end
    
    # face loop
    for fi=1:nface
        e1 = mesh.face2element[1,fi];
        e2 = mesh.face2element[2,fi];
        if e1==0 # boundary face
            e1=e2;
        elseif e2==0
            e2=e1;
        end
        
        for ni=1:nfp
            for nj=1:np
                #e1
                if cggrid.loc2glb[nj,e1] == cggrid.face2glb[ni,fi]
                    dgface2glb[ni,1,fi] = dgloc2glb[nj,e1];
                    dgface2local[ni,1,fi] = nj;
                end
                #e2
                if cggrid.loc2glb[nj,e2] == cggrid.face2glb[ni,fi]
                    dgface2glb[ni,2,fi] = dgloc2glb[nj,e2];
                    dgface2local[ni,2,fi] = nj;
                end
            end
        end
        for ni=1:nfv
            for nj=1:np
                #e1
                if cggrid.loc2glb[nj,e1] == cggrid.faceVertex2glb[ni,fi]
                    dgfaceVertex2glb[ni,1,fi] = dgloc2glb[nj,e1];
                end
                #e2
                if cggrid.loc2glb[nj,e2] == cggrid.faceVertex2glb[ni,fi]
                    dgfaceVertex2glb[ni,2,fi] = dgloc2glb[nj,e2];
                end
            end
        end
        
        dgfaceRefelInd[1,fi] = which_face(dgface2glb[:,1,fi], dgloc2glb[:,e1], dim, nfp);
        dgfaceRefelInd[2,fi] = which_face(dgface2glb[:,2,fi], dgloc2glb[:,e2], dim, nfp);
        
    end
    
    # bdry loop
    for bid=1:length(dgbids)
        for bi=1:length(dgbdryface[bid])
            dgbdry[bid][((bi-1)*nfp+1):(bi*nfp)] = dgface2glb[:,1,dgbdryface[bid][bi]];
            for fcind=1:nfp
                dgbdrynorm[bid][:, (bi-1)*nfp+fcind] = mesh.normals[:, dgbdryface[bid][bi]];
            end
        end
        # set duplicates to zero
        nbp = length(dgbdry[bid])
        nzeros = 0;
        for bi=1:nbp
            ind = dgbdry[bid][bi];
            if ind > 0
                for bj=(bi+1):nbp
                    if dgbdry[bid][bj] == ind
                        dgbdry[bid][bj] = 0;
                        nzeros = nzeros + 1;
                    end
                end
            end
        end
        # remove zeros from the list
        # println(nbp);
        # println(nzeros);
        newdgbdry = zeros(Int, nbp - nzeros);
        newdgbdrynorm = zeros(dim, nbp - nzeros);
        ind = 0;
        for bi=1:nbp
            if dgbdry[bid][bi] > 0
                ind = ind + 1;
                newdgbdry[ind] = dgbdry[bid][bi];
                newdgbdrynorm[:,ind] = dgbdrynorm[bid][:,bi];
            end
        end
        dgbdry[bid] = newdgbdry;
        dgbdrynorm[bid] = newdgbdrynorm;
    end
    
    return DG_Grid(dgnodes, dgbdry, dgbdryface, dgbdrynorm, dgbids, dgloc2glb, dgglbvertex, dgface2glb, dgfaceVertex2glb, dgface2local, dgfacenorm, dgfaceRefelInd);
    
end

# Determine which face of the element this is.
# Note: currently only works for line/quad/hex
function which_face(f2g, l2g, dim, nfp)
    # To determine face refel index, look at first two nodes
    ind = 0;
    if dim == 1
        if f2g[1] == l2g[1]
            # left face
            ind = 1;
        else
            # right face
            ind = 2;
        end
        
    elseif dim == 2
        if f2g[1] == l2g[1]
            if f2g[2] == l2g[2]
                # bottom
                ind = 2;
            else
                # left
                ind = 1;
            end
        else
            if f2g[1] == l2g[nfp]
                # right
                ind = 3;
            else
                # top
                ind = 4;
            end
        end
        
    elseif dim == 3
        if f2g[1] == l2g[1]
            if f2g[2] == l2g[2]
                if f2g[nfp] == l2g[nfp]
                    # front
                    ind = 3;
                else
                    # bottom
                    ind = 2;
                end
            else
                # left
                ind = 1;
            end
        else
            if f2g[nfp-1] == l2g[end-1]
                if f2g[1] == l2g[end-nfp+1]
                    # back
                    ind = 6;
                else
                    # top
                    ind = 5;
                end
            else
                # right
                ind = 4;
            end
        end
    end
    
    return ind;
end