#=
# Contains info about all nodes on the domain
=#
struct Grid
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
    face2glb::Array{Int,2}          # local to global map for faces
    faceVertex2glb::Array{Int,2}    # global indices of face vertices
end

etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces

# Build a grid from a mesh
function grid_from_mesh(mesh)
    if config.dimension == 1
        return grid_from_mesh_1d(mesh);
    elseif config.dimension == 2
        return grid_from_mesh_2d(mesh);
    elseif config.dimension == 3
        return grid_from_mesh_3d(mesh);
    end
end

function grid_from_mesh_1d(mesh)
    ord = config.basis_order_min;
    nfaces = etypetonf[mesh.etypes[1]];
    Nf = size(mesh.face2vertex, 2);
    nx = mesh.nx;
    nel = mesh.nel;
    
    refel = build_refel(1, ord, nfaces, config.elemental_nodes);
    refelfc = build_refel(0, ord, 0, config.elemental_nodes);
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
    
    
    return (refel, refelfc, Grid(x, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb));
end

function grid_from_mesh_2d(mesh)
    
end

function grid_from_mesh_3d(mesh)
    
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