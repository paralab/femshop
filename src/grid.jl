#=
# Contains info about all nodes on the domain
# Unlike MeshData struct, this accounts for interior nodes and corresponds to nodal DOFs.
# This is a CG grid. There is a separate DGGrid struct.
=#
struct Grid
    allnodes::Array{Float64}        # All node coordinates size = (dim, nnodes)
    # boundaries
    bdry::Array{Array{Int,1},1}     # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
    bdryface::Array{Array{Int,1},1} # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
    bdrynorm::Array{Array{Float64,2},1} # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, faces])*note:array of arrays
    bids::Array{Int,1}              # BID corresponding to rows of bdrynodes
    # elements
    loc2glb::Array{Int,2}           # local to global map for each element's nodes (size is (Np, nel))
    glbvertex::Array{Int,2}         # global indices of each elements' vertices (size if (Nvertex, nel))
    # faces
    face2glb::Array{Int,2}          # local to global map for faces (size is (Nfp, Nfaces))
    faceVertex2glb::Array{Int,2}    # global indices of face vertices (size is (Nfvertex, Nfaces))
    #facenormals::Array{Array{Float64,2},1}    # normal vector for each face
end

etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types

# Build a grid from a mesh
function grid_from_mesh(mesh)
    if config.dimension == 1
        return grid_from_mesh_1d(mesh);
    elseif config.dimension == 2
        if mesh.etypes[1] == 2 # triangles
            return grid_from_mesh_2d_triangle(mesh);
        elseif mesh.etypes[1] == 3 # quads
            return grid_from_mesh_2d_quad(mesh);
        end
    elseif config.dimension == 3
        return grid_from_mesh_3d(mesh); # hex and tet
    end
end

function grid_from_mesh_1d(mesh)
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]];
    Nf = size(mesh.face2vertex, 2);
    nx = mesh.nx;
    nel = mesh.nel;
    
    refel = build_refel(1, ord, nfaces, config.elemental_nodes);
    
    N = (nx-1)*ord + 1;         # number of total nodes
    Np = refel.Np;              # number of nodes per element
    x = zeros(1,N);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryfc = [];                # index of elements touching each BID
    bdrynorm = [];              # normal at boundary nodes
    bids = collectBIDs(mesh);
    loc2glb = zeros(Int, Np, nel)# local to global index map for each element's nodes
    glbvertex = zeros(Int, 2, nel);# local to global for vertices
    f2glb = zeros(Int, 1, Nf);# face local to global
    fvtx2glb = zeros(Int, 1, Nf);# face vertex local to global

    # Elements/nodes
    for ei=1:mesh.nel
        elem = mesh.elements[:, ei];
        vx = mesh.nodes[:,elem];
        x1 = vx[1]; # left vertex
        h = vx[2]-vx[1]; # size of element
        glbvertex[1, ei] = (ei-1)*(Np-1) + 1;
        glbvertex[2, ei] = ei*(Np-1) + 1;
        
        for ni=1:Np
            gi = (ei-1)*(Np-1) + ni; # global index of this node
            x[1,gi] = x1 .+ h*0.5 .* (refel.r[ni] + 1); # coordinates of this node
            loc2glb[ni,ei] = gi; # local to global map
        end
        
        # Do the faces while we're here
        f1 = mesh.element2face[1,ei];
        f2 = mesh.element2face[2,ei];
        f2glb[1,f1] = glbvertex[1,ei];
        f2glb[1,f2] = glbvertex[2,ei];
        fvtx2glb[1,f1] = glbvertex[1,ei];
        fvtx2glb[1,f2] = glbvertex[2,ei];
    end
    
    # boundary
    # 
    bdry = Array{Array{Int,1},1}(undef,length(bids));
    bdryfc = Array{Array{Int,1},1}(undef,length(bids));
    bdrynorm = Array{Array{Float64,2},1}(undef,length(bids));
    for i=1:length(bids)
        bdry[i] = Array{Int,1}(undef,0);
        bdryfc[i] = Array{Int,1}(undef,0);
        bdrynorm[i] = Array{Float64,2}(undef,1,0);
    end
    for i=1:Nf
        if mesh.bdryID[i] > 0
            # face i is a boundary face
            # add its nodes to the bdry list
            push!(bdry[mesh.bdryID[i]], f2glb[1,i]);
            push!(bdryfc[mesh.bdryID[i]], i);
            #push!(bdrynorm[mesh.bdryID[i]], [mesh.normals[i]]);
            bdrynorm[mesh.bdryID[i]] = [bdrynorm[mesh.bdryID[i]] mesh.normals[:,i]];
        end
    end
    
    return (refel, Grid(x, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb));
end

function grid_from_mesh_2d_quad(mesh)
    ord = config.basis_order_min;
    nfaces = 4;
    totalfaces = nfaces*mesh.nel;
    nx = mesh.nx;
    nel = mesh.nel;
    
    refel = build_refel(2, ord, 4, config.elemental_nodes);
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(Int, 0));
        push!(bdryfc, zeros(Int,0));
        push!(bdrynorm, zeros(config.dimension,0));
    end
    loc2glb = zeros(Int, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(Int, 4, nel);     # local to global for vertices
    f2glb = zeros(Int, refel.Nfp[1], totalfaces);  # face node local to global
    fvtx2glb = zeros(Int, 2, totalfaces);   # face vertex local to global
    
    tmpallnodes = zeros(2, mesh.nel*refel.Np);
    for ei=1:nel
        # Find this element's nodes
        e_vert = mesh.nodes[1:2, mesh.elements[1:4, ei]];
        (e_x, e_y) = quad_refel_to_xy(refel.r[:,1], refel.r[:,2], e_vert);
        
        # Add them to the tmp global nodes
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] = e_x;
        tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] = e_y;
        
        # temporary mapping
        loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
    end
    
    # Go back and remove duplicate nodes. Adjust loc2glb.
    to_remove = [];
    remove_count = 0;
    tol = 1e-9;
    found = false;
    next_ind = Np+1;
    allnodes = zeros(size(tmpallnodes));
    allnodes[:,1:Np] = tmpallnodes[:,1:Np];
    for ei=2:nel
        for ni=1:Np
            found = false;
            for ej=1:ei-1
                for nj=1:Np
                    if is_same_node(tmpallnodes[:,loc2glb[ni,ei]], allnodes[:,loc2glb[nj,ej]], tol)
                        # duplicates: keep the ej one, remove ei
                        push!(to_remove, loc2glb[ni,ei]);
                        loc2glb[ni,ei] = loc2glb[nj,ej];
                        remove_count += 1;
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found
                allnodes[:,next_ind] = tmpallnodes[:,loc2glb[ni,ei]];
                loc2glb[ni,ei] = next_ind;
                next_ind += 1;
            end
        end
    end
    N = next_ind-1;
    allnodes = allnodes[:,1:N];
    
    # face to global
    for ei=1:nel
        for fi=1:nfaces
            fid = (ei-1)*nfaces + fi;
            f2glb[:,fid] = loc2glb[refel.face2local[fi], ei];
        end
    end
    
    # vertices and boundary
    for ei=1:nel
        mfids = mesh.element2face[:,ei];
        normals = mesh.normals[:,mfids];
        
        # vertices
        for ni=1:Np
            for vi=1:4
                if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol)
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        
        # mesh and grid face indices may be different, so we need to map them
        gfids = ((ei-1)*nfaces+1):(ei*nfaces);
        g2mfids = [0,0,0,0];
        meshfaces = mesh.element2face[:,ei];
        for gfi=1:4
            for mfi=1:4
                if is_same_line(mesh.nodes[:,mesh.face2vertex[:,meshfaces[mfi]]], allnodes[:, f2glb[:,gfids[gfi]]], tol)
                    g2mfids[gfi] = mfi;
                end
            end
        end
        # println("grid2mesh faces for element "*string(ei)*":")
        # println(g2mfids)
        
        # Now that we have a mapping to mesh faces, copy boundary info
        # bdry, bdryface, bdrynorm
        for fi=1:4
            b = mesh.bdryID[mesh.element2face[g2mfids[fi],ei]];
            bind = indexin([b], bids)[1];
            if !(bind === nothing)
                append!(bdry[bind], f2glb[:,gfids[fi]]);
                push!(bdryfc[bind], gfids[fi]);
                bdrynorm[bind] = hcat(bdrynorm[bind], normals[:,indexin([fi],g2mfids)]);
            end
        end
        
    end # element loop
    
    # There are duplicates in the bdry info. Remove them
    newbdry = similar(bdry);
    for i=1:length(bdry)
        newbdry[i] = [];
    end
    
    for bidi=1:length(bids)
        for bi=1:length(bdry[bidi])
            found = false;
            for bidj=1:bidi
                for bj=1:length(newbdry[bidj])
                    if bdry[bidi][bi] == newbdry[bidj][bj]
                        # bi already exists in newbdry
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found # it's a new bdry node
                push!(newbdry[bidi], bdry[bidi][bi]);
            end
        end
    end
    bdry = newbdry;
    
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb));
end

function grid_from_mesh_2d_triangle(mesh)
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]];
    totalfaces = nfaces*mesh.nel;
    nx = mesh.nx;
    nel = mesh.nel;
    
    refel = build_refel(2, ord, nfaces, config.elemental_nodes);
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(Int, 0));
        push!(bdryfc, zeros(Int,0));
        push!(bdrynorm, zeros(config.dimension,0));
    end
    loc2glb = zeros(Int, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(Int, 3, nel);     # local to global for vertices
    f2glb = zeros(Int, refel.Nfp[1], totalfaces);  # face node local to global
    fvtx2glb = zeros(Int, 2, totalfaces);   # face vertex local to global
    
    
    tmpallnodes = zeros(2, mesh.nel*refel.Np);
    for ei=1:nel
        # Find this element's nodes
        e_vert = mesh.nodes[1:2, mesh.elements[1:3, ei]];
        (e_x, e_y) = triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], e_vert);
        
        # Add them to the tmp global nodes
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] = e_x;
        tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] = e_y;
        
        # temporary mapping
        loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
    end
    
    # Go back and remove duplicate nodes. Adjust loc2glb.
    to_remove = [];
    remove_count = 0;
    tol = 1e-9;
    found = false;
    next_ind = Np+1;
    allnodes = zeros(size(tmpallnodes));
    allnodes[:,1:Np] = tmpallnodes[:,1:Np];
    for ei=2:nel
        for ni=1:Np
            found = false;
            for ej=1:ei-1
                for nj=1:Np
                    if is_same_node(tmpallnodes[:,loc2glb[ni,ei]], allnodes[:,loc2glb[nj,ej]], tol)
                        # duplicates: keep the ej one, remove ei
                        push!(to_remove, loc2glb[ni,ei]);
                        loc2glb[ni,ei] = loc2glb[nj,ej];
                        remove_count += 1;
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found
                allnodes[:,next_ind] = tmpallnodes[:,loc2glb[ni,ei]];
                loc2glb[ni,ei] = next_ind;
                next_ind += 1;
            end
        end
    end
    N = next_ind-1;
    allnodes = allnodes[:,1:N];

    # face to global
    for ei=1:nel
        for fi=1:nfaces
            fid = (ei-1)*nfaces + fi;
            f2glb[:,fid] = loc2glb[refel.face2local[fi], ei];
        end
    end
    
    # vertices and boundary
    for ei=1:nel
        mfids = mesh.element2face[:,ei];
        normals = mesh.normals[:,mfids];
        
        # vertices
        for ni=1:Np
            for vi=1:3
                if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol)
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        
        # mesh and grid face indices may be different, so we need to map them
        gfids = ((ei-1)*nfaces+1):(ei*nfaces);
        g2mfids = [0,0,0];
        meshfaces = mesh.element2face[:,ei];
        for gfi=1:3
            for mfi=1:3
                if is_same_line(mesh.nodes[:,mesh.face2vertex[:,meshfaces[mfi]]], allnodes[:, f2glb[:,gfids[gfi]]], tol)
                    g2mfids[gfi] = mfi;
                end
            end
        end
        # println("grid2mesh faces for element "*string(ei)*":")
        # println(g2mfids)
        
        # Now that we have a mapping to mesh faces, copy boundary info
        # bdry, bdryface, bdrynorm
        for fi=1:3
            b = mesh.bdryID[mesh.element2face[g2mfids[fi],ei]];
            bind = indexin([b], bids)[1];
            if !(bind === nothing)
                append!(bdry[bind], f2glb[:,gfids[fi]]);
                push!(bdryfc[bind], gfids[fi]);
                bdrynorm[bind] = hcat(bdrynorm[bind], normals[:,indexin([fi],g2mfids)]);
            end
        end
        
    end # element loop 
    # There are duplicates in the bdry info. Remove them
    newbdry = similar(bdry);
    for i=1:length(bdry)
        newbdry[i] = [];
    end
    
    for bidi=1:length(bids)
        for bi=1:length(bdry[bidi])
            found = false;
            for bidj=1:bidi
                for bj=1:length(newbdry[bidj])
                    if bdry[bidi][bi] == newbdry[bidj][bj]
                        # bi already exists in newbdry
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found # it's a new bdry node
                push!(newbdry[bidi], bdry[bidi][bi]);
            end
        end
    end
    bdry = newbdry;
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb));
end

function grid_from_mesh_3d(mesh)
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]];
    totalfaces = nfaces*mesh.nel;
    nx = mesh.nx;
    nel = mesh.nel;
    if nfaces == 6
        nvtx = 8;
        facenvtx = 4;
    else
        nvtx = 4;
        facenvtx = 3;
    end
    
    refel = build_refel(3, ord, nfaces, config.elemental_nodes);
    
    Np = refel.Np;                      # number of nodes per element
    bdry = [];                          # index(in x) of boundary nodes for each BID
    bdryfc = [];                        # index of faces touching each BID
    bdrynorm = [];                      # normal at boundary nodes
    bids = collectBIDs(mesh);           # BID list
    nbids = length(bids);
    for i=1:nbids
        push!(bdry, zeros(Int, 0));
        push!(bdryfc, zeros(Int,0));
        push!(bdrynorm, zeros(config.dimension,0));
    end
    loc2glb = zeros(Int, Np, nel)       # local to global index map for each element's nodes
    glbvertex = zeros(Int, nvtx, nel);     # local to global for vertices
    f2glb = zeros(Int, refel.Nfp[1], totalfaces);  # face node local to global
    fvtx2glb = zeros(Int, facenvtx, totalfaces);   # face vertex local to global
    
    tmpallnodes = zeros(3, mesh.nel*refel.Np);
    for ei=1:nel
        # Find this element's nodes
        e_vert = mesh.nodes[1:3, mesh.elements[1:8, ei]];
        if nvtx == 8
            (e_x, e_y, e_z) = hex_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
        else
            (e_x, e_y, e_z) = tet_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
        end
        
        # Add them to the tmp global nodes
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] = e_x;
        tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] = e_y;
        tmpallnodes[3, ((ei-1)*Np+1):(ei*Np)] = e_z;
        
        # temporary mapping
        loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
    end
    
    # Go back and remove duplicate nodes. Adjust loc2glb.
    to_remove = [];
    remove_count = 0;
    tol = 1e-9;
    found = false;
    next_ind = Np+1;
    allnodes = zeros(size(tmpallnodes));
    allnodes[:,1:Np] = tmpallnodes[:,1:Np];
    for ei=2:nel
        for ni=1:Np
            found = false;
            for ej=1:ei-1
                for nj=1:Np
                    if is_same_node(tmpallnodes[:,loc2glb[ni,ei]], allnodes[:,loc2glb[nj,ej]], tol)
                        # duplicates: keep the ej one, remove ei
                        push!(to_remove, loc2glb[ni,ei]);
                        loc2glb[ni,ei] = loc2glb[nj,ej];
                        remove_count += 1;
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found
                allnodes[:,next_ind] = tmpallnodes[:,loc2glb[ni,ei]];
                loc2glb[ni,ei] = next_ind;
                next_ind += 1;
            end
        end
    end
    N = next_ind-1;
    allnodes = allnodes[:,1:N];
    
    # face to global
    for ei=1:nel
        for fi=1:nfaces
            fid = (ei-1)*nfaces + fi;
            f2glb[:,fid] = loc2glb[refel.face2local[fi], ei];
        end
    end
    
    # vertices and boundary
    for ei=1:nel
        mfids = mesh.element2face[:,ei];
        normals = mesh.normals[:,mfids];
        
        # vertices
        for ni=1:Np
            for vi=1:8
                if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol)
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        
        # mesh and grid face indices may be different, so we need to map them
        gfids = ((ei-1)*nfaces+1):(ei*nfaces);
        g2mfids = zeros(Int, nfaces);
        meshfaces = mesh.element2face[:,ei];
        for gfi=1:nfaces
            for mfi=1:nfaces
                if is_same_plane(mesh.nodes[:,mesh.face2vertex[:,meshfaces[mfi]]], allnodes[:, f2glb[:,gfids[gfi]]], tol)
                    g2mfids[gfi] = mfi;
                end
            end
        end
        # println("grid2mesh faces for element "*string(ei)*":")
        # println(g2mfids)
        
        # Now that we have a mapping to mesh faces, copy boundary info
        # bdry, bdryface, bdrynorm
        for fi=1:nfaces
            b = mesh.bdryID[mesh.element2face[g2mfids[fi],ei]];
            bind = indexin([b], bids)[1];
            if !(bind === nothing)
                append!(bdry[bind], f2glb[:,gfids[fi]]);
                push!(bdryfc[bind], gfids[fi]);
                bdrynorm[bind] = hcat(bdrynorm[bind], normals[:,indexin([fi],g2mfids)]);
            end
        end
        
    end # element loop
    
    # There are duplicates in the bdry info. Remove them
    newbdry = similar(bdry);
    for i=1:length(bdry)
        newbdry[i] = [];
    end
    
    for bidi=1:length(bids)
        for bi=1:length(bdry[bidi])
            found = false;
            for bidj=1:bidi
                for bj=1:length(newbdry[bidj])
                    if bdry[bidi][bi] == newbdry[bidj][bj]
                        # bi already exists in newbdry
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if !found # it's a new bdry node
                push!(newbdry[bidi], bdry[bidi][bi]);
            end
        end
    end
    bdry = newbdry;
    
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb));
end

function collectBIDs(mesh)
    bids = [];
    for i=1:length(mesh.bdryID)
        if mesh.bdryID[i] > 0
            already = false;
            for j=1:length(bids)
                if mesh.bdryID[i] == bids[j]
                    already = true;
                end
            end
            if !already
                push!(bids, mesh.bdryID[i]);
            end
        end
    end
    return bids;
end

#Extra remove later 

function triangle_element_nodes_(refel, v)
    return  triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], v);
end

function triangle_refel_to_xy_(r, s, v)
    x = 0.5 * (-(r .+ s) * v[1,1] .+ (1 .+ r) * v[1,2] .+ (1 .+ s) * v[1,3]);
    y = 0.5 * (-(r .+ s) * v[2,1] .+ (1 .+ r) * v[2,2] .+ (1 .+ s) * v[2,3]);
    
    return (x, y);
end

function quad_refel_to_xy(r, s, v)
    vx = v[1,:];
    vy = v[2,:];
    
    rp = 1 .+ r;
    rm = 1 .- r;
    sp = 1 .+ s;
    sm = 1 .- s;
    x = 0.25 .* (rm .* sm .* vx[1] .+ rp .* sm .* vx[2] .+ rp .* sp .* vx[3] .+ rm .* sp .* vx[4]); 
    y = 0.25 .* (rm .* sm .* vy[1] .+ rp .* sm .* vy[2] .+ rp .* sp .* vy[3] .+ rm .* sp .* vy[4]); 
    
    return (x, y);
end

function hex_refel_to_xyz(r, s, t, v)
    vx = v[1,:];
    vy = v[2,:];
    vz = v[3,:];
    
    rp = 1 .+ r;
    rm = 1 .- r;
    sp = 1 .+ s;
    sm = 1 .- s;
    tp = 1 .+ t;
    tm = 1 .- t;
    x = 0.125 .* (rm .* sm .* tm .* vx[1] .+ rp .* sm .* tm .* vx[2] .+ rp .* sp .* tm .* vx[3] .+ rm .* sp .* tm .* vx[4]
                    .+ rm .* sm .* tp .* vx[5] .+ rp .* sm .* tp .* vx[6] .+ rp .* sp .* tp .* vx[7] .+ rm .* sp .* tp .* vx[8]);
    y = 0.125 .* (rm .* sm .* tm .* vy[1] .+ rp .* sm .* tm .* vy[2] .+ rp .* sp .* tm .* vy[3] .+ rm .* sp .* tm .* vy[4]
                    .+ rm .* sm .* tp .* vy[5] .+ rp .* sm .* tp .* vy[6] .+ rp .* sp .* tp .* vy[7] .+ rm .* sp .* tp .* vy[8]);
    z = 0.125 .* (rm .* sm .* tm .* vz[1] .+ rp .* sm .* tm .* vz[2] .+ rp .* sp .* tm .* vz[3] .+ rm .* sp .* tm .* vz[4]
                    .+ rm .* sm .* tp .* vz[5] .+ rp .* sm .* tp .* vz[6] .+ rp .* sp .* tp .* vz[7] .+ rm .* sp .* tp .* vz[8]);
    
    return (x, y, z);
end

function tetrahedron_refel_to_xyz(r, s, t, v)
    x = 0.5 .* (-(r.+s.+t) .* 2/3*v[1,1] .+ (1 .+r) .* v[1,2] .+ (1 .+s) .* v[1,3] .+ (1 .+t) .* v[1,4]);
    y = 0.5 .* (-(r.+s.+t) .* 2/3*v[2,1] .+ (1 .+r) .* v[2,2] .+ (1 .+s) .* v[2,3] .+ (1 .+t) .* v[2,4]);
    z = 0.5 .* (-(r.+s.+t) .* 2/3*v[3,1] .+ (1 .+r) .* v[3,2] .+ (1 .+s) .* v[3,3] .+ (1 .+t) .* v[3,4]);
    
    return (x, y, z);
end

# Returns true if the nodes are within tol of each other.
function is_same_node(x1, x2, tol)
    return sum(abs.(x1 - x2)) < tol
end

# Returns true if the two node lists have at least two of the same nodes.
function is_same_line(l1, l2, tol)
    found = 0;
    n1 = size(l1,2);
    n2 = size(l2,2);
    for i=1:n1
        for j=1:n2
            if is_same_node(l1[:,i], l2[:,j], tol)
                found += 1;
            end
            if found >= 2
                return true;
            end
        end
    end
    
    return false;
end

# Returns true if the two node lists have at least three of the same nodes.
function is_same_plane(p1, p2, tol)
    found = 0;
    n1 = size(p1,2);
    n2 = size(p2,2);
    for i=1:n1
        for j=1:n2
            if is_same_node(p1[:,i], p2[:,j], tol)
                found += 1;
            end
            if found >= 3
                return true;
            end
        end
    end
    
    return false;
end