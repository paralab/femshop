#=
# Contains info about all nodes on the domain
# Unlike MeshData struct, this accounts for interior nodes and corresponds to nodal DOFs.
=#

struct Grid
    allnodes::Array{Float64}        # All node coordinates size = (dim, nnodes)
    # boundaries
    bdry::Array{Array{Int,1},1}     # Indices of boundary nodes for each BID (bdry[bid][nodes])*note:array of arrays
    bdryface::Array{Array{Int,1},1} # Indices of faces touching each BID (bdryface[bid][faces])*note:array of arrays
    bdrynorm::Array{Array{Float64,2},1} # Normal vector for boundary nodes for each BID (bdrynorm[bid][dim, nodes])*note:array of arrays
    bids::Array{Int,1}              # BID corresponding to rows of bdrynodes
    # elements
    loc2glb::Array{Int,2}           # local to global map for each element's nodes (size is (Np, nel))
    glbvertex::Array{Int,2}         # global indices of each elements' vertices (size is (Nvertex, nel))
    # faces (For CG, G=1. For DG, G=2)
    face2glb::Array{Int,3}          # local to global map for faces (size is (Nfp, G, Nfaces))
    element2face::Array{Int,2}      # face indices for each element (size is (Nfaces, nel))
    face2element::Array{Int,2}      # elements on both sides of a face, 0=boundary (size is (2, Nfaces))
    facenormals::Array{Float64,2}   # normal vector for each face
    faceRefelInd::Array{Int,2}      # Index for face within the refel for each side
    
    facebid::Array{Int,1}           # BID of each face (0=interior face)
end

etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5]; # number of vertices for each type
etypetodim= [1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3]; # dimension of each type
etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for element types
etypetoftype=[0,1, 1, 2, 3, 3, 3, 0, 1, 1, 2, 3, 3, 3, 0, 0, 0, 0, 0]; # type of faces for this element type

# Build a grid from a mesh
function grid_from_mesh(mesh)
    dim = config.dimension;
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]];;
    #totalfaces = nfaces*mesh.nel;
    totalfaces = size(mesh.normals,2);
    nel = mesh.nel;
    if dim == 1
        facenvtx = 1
    else
        facenvtx = etypetonv[etypetoftype[mesh.etypes[1]]]; # Assumes one element type
    end
    nvtx = etypetonv[mesh.etypes[1]]; # Assumes one element type
    
    refel = build_refel(dim, ord, nfaces, config.elemental_nodes);
    
    if config.solver_type == DG # || config.solver_type == FV ????
        Gness = 2;
    else
        Gness = 1;
    end
    
    tol = 1e-12; # tolerance for is_same_node
    
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
    f2glb = zeros(Int, refel.Nfp[1], Gness, totalfaces);  # face node local to global
    element2face = zeros(Int, nfaces, nel);  # element to face map
    face2element = zeros(Int, 2, size(mesh.face2element,2));  # face to element map
    facenormals = zeros(dim, totalfaces); # normal vectors for every face
    faceRefelInd = zeros(Int, 2, totalfaces); # Index in refel for this face for elements on both sides
    facebid = zeros(Int, totalfaces); # BID of each face
    
    tmpallnodes = zeros(dim, mesh.nel*refel.Np);
    for ei=1:nel
        # Find this element's nodes
        n_vert = etypetonv[mesh.etypes[ei]];
        e_vert = mesh.nodes[1:dim, mesh.elements[1:n_vert, ei]];
        
        if dim == 1
            e_x = line_refel_to_x(refel.r[:,1], e_vert);
        elseif dim == 2
            if nfaces == 3 # triangles
                (e_x, e_y) = triangle_refel_to_xy_(refel.r[:,1], refel.r[:,2], e_vert);
            else # quads
                (e_x, e_y) = quad_refel_to_xy(refel.r[:,1], refel.r[:,2], e_vert);
            end
        elseif dim == 3
            if nvtx == 8 # hexes
                (e_x, e_y, e_z) = hex_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
            else # tets
                (e_x, e_y, e_z) = tetrahedron_refel_to_xyz(refel.r[:,1], refel.r[:,2], refel.r[:,3], e_vert);
            end
        end
        
        # Add them to the tmp global nodes
        tmpallnodes[1, ((ei-1)*Np+1):(ei*Np)] = e_x;
        if dim > 1
            tmpallnodes[2, ((ei-1)*Np+1):(ei*Np)] = e_y;
            if dim > 2
                tmpallnodes[3, ((ei-1)*Np+1):(ei*Np)] = e_z;
            end
        end
        
        # temporary mapping
        loc2glb[:,ei] = ((ei-1)*Np+1):(ei*Np);
    end
    
    if Gness == 1 # CG
        # Go back and remove duplicate nodes. Adjust loc2glb.
        remove_count = 0;
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
        
    else
        allnodes = tmpallnodes; # DG grid is already made
    end
    
    # if config.solver_type == FV
    #     centers = zeros(dim, nel);# For finite volume, store cell centers
    # else
    #     centers = zeros(0);
    # end
    
    # vertices, faces and boundary
    for ei=1:nel
        n_vert = etypetonv[mesh.etypes[ei]];
        mfids = mesh.element2face[:,ei];
        normals = mesh.normals[:,mfids];
        el_center = zeros(dim);
        
        # vertices and center
        for ni=1:Np
            el_center .+= allnodes[:,loc2glb[ni,ei]];
            
            for vi=1:n_vert
                if is_same_node(mesh.nodes[:, mesh.elements[vi,ei]], allnodes[:,loc2glb[ni,ei]], tol)
                    glbvertex[vi, ei] = loc2glb[ni,ei];
                end
            end
        end
        el_center ./= Np;
        # # For finite volumes, store the cell center
        # if config.solver_type == FV
        #     centers[:, ei] = el_center;
        # end
        
        # f2glb has duplicates. Compare to mesh faces and keep same ordering as mesh.
        # Copy normals and bdry info.
        # Set element2face map.
        meshfaces = mesh.element2face[:,ei];    # index from mesh
        if dim == 1
            test_same_face = is_same_node;
        elseif dim == 2
            test_same_face = is_same_line;
        elseif dim == 3
            test_same_face = is_same_plane;
        end
        
        for gfi=1:nfaces
            tmpf2glb = loc2glb[refel.face2local[gfi], ei];
            
            for mfi=1:nfaces
                thisfaceind = meshfaces[mfi];
                if test_same_face(mesh.nodes[:,mesh.face2vertex[:,thisfaceind]], allnodes[:, tmpf2glb], tol)
                    # This mesh face corresponds to this tmpf2glb face
                    # Put the tmpf2glb map into f2glb at the mesh index(thisfaceind).
                    # Move the f2glb[:,1,ind] to f2glb[:,2,ind] first if DG (Gness==2)
                    # Set element2face according to gfi(not mfi)
                    # Move face2element[1] to face2element[2] and put this one in [1]
                    if (Gness == 2) f2glb[:, 2, thisfaceind] = f2glb[:, 1, thisfaceind]; end
                    f2glb[:, 1, thisfaceind] = tmpf2glb;
                    element2face[gfi, ei] = thisfaceind;
                    face2element[2, thisfaceind] = face2element[1, thisfaceind];
                    face2element[1, thisfaceind] = ei;
                    
                    # Find the normal for every face. The normal points from e1 to e2 or outward for boundary.
                    # Note that the normal stored in mesh_data could be pointing either way.
                    thisnormal = normals[:, mfi];
                    f_center = zeros(dim);
                    for ni=1:length(tmpf2glb)
                        f_center .+= allnodes[:, tmpf2glb[ni]];
                    end
                    f_center ./= length(tmpf2glb)
                    d1 = norm(f_center .+ thisnormal .- el_center);
                    d2 = norm(f_center .- thisnormal .- el_center);
                    if d1 < d2 # normal is pointing in the wrong direction
                        thisnormal = -thisnormal;
                    end
                    facenormals[:, thisfaceind] = thisnormal;
                    
                    # Copy boundary info: bdry, bdryface, bdrynorm
                    mbid = mesh.bdryID[thisfaceind];
                    gbid = indexin([mbid], bids)[1];
                    nfacenodes = length(tmpf2glb);
                    if !(gbid === nothing) # This is a boundary face
                        append!(bdry[gbid], tmpf2glb);
                        push!(bdryfc[gbid], thisfaceind);
                        facebid[thisfaceind] = gbid;
                        thisnormal = normals[:, mfi];
                        normchunk = zeros(config.dimension, nfacenodes);
                        for ni=1:nfacenodes
                            normchunk[:,ni] = thisnormal;
                        end
                        bdrynorm[gbid] = hcat(bdrynorm[gbid], normchunk);
                    end
                end
            end
        end
        
    end # element loop
    
    # There are duplicates in the bdry info. Remove them
    newbdry = similar(bdry);
    newbdrynorm = similar(bdrynorm);
    for i=1:length(bdry)
        newbdry[i] = zeros(length(bdry[i]));
        newbdrynorm[i] = zeros(size(bdrynorm[i]));
    end
    
    for bidi=1:length(bids)
        nextbdryind = 1;
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
                newbdry[bidi][nextbdryind] = bdry[bidi][bi];
                newbdrynorm[bidi][:,nextbdryind] = bdrynorm[bidi][:,bi];
                nextbdryind += 1;
            end
        end
        newbdry[bidi] = newbdry[bidi][1:nextbdryind-1];
        newbdrynorm[bidi] = newbdrynorm[bidi][:, 1:nextbdryind-1];
    end
    bdry = newbdry;
    bdrynorm = newbdrynorm;
    
    # Refel index for each face
    for fi=1:totalfaces
        faceRefelInd[1,fi] = which_refel_face(f2glb[:,1,fi], allnodes, refel, loc2glb[:,face2element[1,fi]]);
        if face2element[2,fi] > 0
            faceRefelInd[2,fi] = which_refel_face(f2glb[:,Gness,fi], allnodes, refel, loc2glb[:,face2element[2,fi]]);
        end
    end
    
    return (refel, Grid(allnodes, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, element2face, face2element, facenormals, faceRefelInd, facebid));
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

function line_refel_to_x(r, v)
    x = v[1] .+ 0.5 .* (v[2]-v[1]) .* (1 .+ r);
    
    return x;
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
    d = v[:,1];
    A = [v[:,2].-d v[:,3].-d v[:,4].-d];
    
    np = length(r);
    mv = zeros(3,np);
    for i=1:np
        tmp = [(r[i]+1)/2, (s[i]+1)/2, (t[i]+1)/2];
        mv[:,i] = A*tmp + d;
    end
    
    x = mv[1,:]
    y = mv[2,:]
    z = mv[3,:]
    
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

function which_refel_face(f2glb, nodes, ref, elem)
    # Check f2glb against the face2local in refel
    fnodes = nodes[:,f2glb];
    for fi=1:ref.Nfaces
        refnodes = nodes[:, elem[ref.face2local[fi]]];
        
        if is_same_face(fnodes, refnodes, ref.dim)
            return fi;
        end
    end
    
    printerr("Couldn't match face when building grid (see which_refel_face() in grid.jl)");
end

function is_same_face(f1, f2, dim)
    tol = 1e-12;
    if dim == 1 # one point
        return is_same_node(f1, f2, tol);
    elseif dim == 2 # same line(two same points)
        return is_same_line(f1, f2, tol);
    elseif dim == 3 # same plane(three same points)
        return is_same_plane(f1, f2, tol);
    end
end

# Adds a boundary ID to some region. Find boundary points satifying on_bdry and moves them to a new set for this bid.
function add_boundary_ID_to_grid(bid, on_bdry, grid)
    # Find if this bid exists. If so, just add points to it, removing from others.
    ind = indexin([bid], grid.bids)[1];
    nbids = length(grid.bids);
    if ind === nothing
        # This is a new bid, add to bids, bdry, bdryface, bdrynorm
        ind = nbids + 1;
        nbids += 1;
        push!(grid.bids, bid);
        push!(grid.bdry, zeros(Int, 0));
        push!(grid.bdryface, zeros(Int, 0));
        push!(grid.bdrynorm, zeros(config.dimension, 0));
    end
    
    # Search all other bids for nodes and faces on this segment. Remove them there and add them here.
    # First find indices and count them. Then move.
    move_nodes = Array{Array{Int,1},1}(undef,nbids);
    node_count = zeros(Int, nbids);
    move_faces = Array{Array{Int,1},1}(undef,nbids);
    face_count = zeros(Int, nbids);
    for i=1:nbids
        bi = grid.bids[i];
        move_nodes[i] = [];
        move_faces[i] = [];
        if bi != bid
            # First the nodes
            for j=1:length(grid.bdry[i])
                nj = grid.bdry[i][j];
                if config.dimension == 1
                    if on_bdry(grid.allnodes[1, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif config.dimension == 2
                    if on_bdry(grid.allnodes[1, nj], grid.allnodes[2, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                elseif config.dimension == 3
                    if on_bdry(grid.allnodes[1, nj], grid.allnodes[2, nj], grid.allnodes[3, nj])
                        push!(move_nodes[i], nj);
                        node_count[i] += 1;
                    end
                end
            end
            # Then the faces
            for j=1:length(grid.bdryface[i])
                fj = grid.bdryface[i][j];
                nfp = size(grid.face2glb,1)
                isbdryface = true
                for ni=1:nfp
                    fx = grid.allnodes[:,grid.face2glb[ni,1,fj]];
                    if config.dimension == 1
                        if !on_bdry(fx[1])
                            isbdryface = false;
                        end
                    elseif config.dimension == 2
                        if !on_bdry(fx[1], fx[2])
                            isbdryface = false;
                        end
                    elseif config.dimension == 3
                        if !on_bdry(fx[1],fx[2],fx[3])
                            isbdryface = false;
                        end
                    end
                    if !isbdryface
                        break;
                    end
                end
                
                if isbdryface
                    push!(move_faces[i], fj);
                    face_count[i] += 1;
                end
            end
        end
    end # find indices
    
    # Move things from other bids to this one
    for i=1:nbids
        if i != ind
            # Add to this bid
            append!(grid.bdry[ind], move_nodes[i]);
            append!(grid.bdryface[ind], move_faces[i]);
            grid.bdrynorm[ind] = hcat(grid.bdrynorm[ind], grid.bdrynorm[i][:,indexin(move_nodes[i], grid.bdry[i])]);
            
            # Make sure all of the norms correspond to the face on this bdry
            startnodeind = length(grid.bdry[ind]) - length(move_nodes[i]);
            for ni=1:length(move_nodes[i])
                for fi=1:length(move_faces[i])
                    # Does this node lie on this face?
                    facenodeindex = indexin(move_nodes[i][ni], grid.face2glb[:,1,move_faces[i][fi]]);
                    if length(facenodeindex) > 0
                        grid.bdrynorm[ind][:,startnodeind + ni] = grid.facenormals[:,move_faces[i][fi]];
                        break;
                    end
                end
            end
            
            # Remove things from other bids
            # Remove bdrynorm and bdryfacenorm
            numremove = length(move_nodes[i]);
            if numremove > 0
                newbdrynorm = zeros(config.dimension, size(grid.bdrynorm[i],2) - numremove);
                nextind = 1;
                for j=1:length(grid.bdry[i])
                    keepit = true;
                    for k=1:numremove
                        if grid.bdry[i][j] == move_nodes[i][k]
                            keepit = false;
                            break;
                        end
                    end
                    if keepit
                        newbdrynorm[:,nextind] = grid.bdrynorm[i][:,j];
                        nextind += 1;
                    end
                end
                grid.bdrynorm[i] = newbdrynorm;
            end
            
            # Remove nodes
            deleteat!(grid.bdry[i], indexin(move_nodes[i], grid.bdry[i]));
            
            # Remove bdryface
            deleteat!(grid.bdryface[i], indexin(move_faces[i], grid.bdryface[i]));
            
            # Change the bid of moved faces in facebid
            for fi=1:length(move_faces[i])
                grid.facebid[move_faces[i][fi]] = bid;
            end
        end
    end
    
    log_entry("Added boundary ID: "*string(bid)*" including "*string(node_count)*" nodes, "*string(face_count)*" faces.");
    
end