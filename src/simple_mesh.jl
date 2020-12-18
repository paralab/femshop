# A set of simple mesh makers for testing
export simple_line_mesh, simple_quad_mesh, simple_hex_mesh

#=
# Builds a 1D interval mesh.
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_line_mesh(nx, bn, interval)
    # mesh only
    Nv = nx;                    # number of vertex nodes
    xv = zeros(1,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = nx-1;                 # number of elements
    el = zeros(Int, 2, nel);    # element vertex maps
    etypes = ones(Int, nel);    # element types (gmsh number)
    numvert = 2*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    f2n = zeros(Int, 1,Nv);     # face2node mapping
    f2e = zeros(Int, 2,Nv);     # face2element mapping
    e2f = zeros(Int, 2,nel);    # element2face mapping
    normals = ones(1,Nv);      # normals of faces
    bdryID = zeros(Int, Nv);    # BID of faces (0=interior)
    
    scale = interval[2]-interval[1];
    h = scale/(nx-1);               # uniformly divided
    
    # vertex nodes are initially ordered lexicographically (that's a 17 letter word!)
    for j=1:nx
        xv[1,j] = interval[1] + (j-1)*h;
    end
    
    # Elements are ordered the same way
    for ei=1:nel
        el[1,ei] = ei;
        el[2,ei] = ei+1;
        e2f[1,ei] = ei;
        e2f[2,ei] = ei+1;
    end
    
    # face2node, face2element, normals
    for fi=1:Nv
        f2n[1,fi] = fi;
        f2e[1,fi] = fi-1;
        f2e[2,fi] = fi;
        # normals[fi] = 1; # already set to 1
        # bdryID[fi] = 0; # already set to 0
    end
    f2e[2,Nv] = 0;
    bdryID[1] = 1;
    f2e[:,1] = [1,0];
    normals[1] = -1;
    if bn == 2
        bdryID[Nv] = 2; # for two BIDs
    else
        bdryID[Nv] = 1; # for one BID
    end
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    ###########################
    # The following is commented because the 1D grid is built automatically in grid.jl
    ###########################
    
    # # grid
    # ord = config.basis_order_min;
    # refel = build_refel(1, ord, 2, LOBATTO);
    # N = (nx-1)*ord + 1;         # number of total nodes
    # Np = refel.Np;              # number of nodes per element
    # x = zeros(N,1);             # coordinates of all nodes
    # bdry = [];                  # index(in x) of boundary nodes for each BID
    # bdryel = [];                # index of elements touching each BID
    # bdrynorm = [];              # normal at boundary nodes
    # if bn == 2
    #     bids = [1,2]; # for two BIDs
    # else
    #     bids = [1]; # for one BID
    # end
    # loc2glb = zeros(Int, nel,Np)# local to global index map for each element's nodes
    # glbvertex = zeros(Int, nel,2);# local to global for vertices

    # # Elements are ordered the same way
    # for ei=1:(nx-1)
    #     x1 = interval[1] + (ei-1)*h; # left vertex
    #     glbvertex[ei,1] = (ei-1)*(Np-1) + 1;
    #     glbvertex[ei,2] = ei*(Np-1) + 1;
        
    #     for ni=1:Np-1
    #         gi = (ei-1)*(Np-1) + ni; # global index of this node
    #         x[gi,1] = x1 .+ h*0.5 .* (refel.r[ni] + 1); # coordinates of this node
    #         loc2glb[ei,ni] = gi; # local to global map
    #     end
    #     loc2glb[ei,Np] = ei*(Np-1) + 1;

    #     el[ei,1] = ei;
    #     el[ei,2] = ei+1;
    # end
    # # The last node was left out
    # x[N] = interval[2];
    
    # # boundary nodes
    # if bn == 2
    #     bdry = [[1],[N]];
    #     bdryel = [[1],[nel]];
    #     bdry1 = Array{Float64,2}(undef,1,1);
    #     bdry1[1] = -1;
    #     bdry2 = Array{Float64,2}(undef,1,1);
    #     bdry2[1] = 1;
    #     bdrynorm = [bdry1,bdry2];
    # else
    #     bdry = [[1,N]]; # for one BID
    #     bdryel = [[1,nel]];
    #     bdry1 = Array{Float64,2}(undef,2,1);
    #     bdry1[1] = -1;
    #     bdry1[2] = 1;
    #     bdrynorm = [bdry1];
    # end
    
    
    # grid = Grid(x, bdry, bdryel, bdrynorm, bids, loc2glb, glbvertex);
    
    (refel, refelfc, grid) = grid_from_mesh_1d(mesh);
    return (mesh, refel, refelfc, grid);
end

#=
# Builds a 2D quad mesh
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_quad_mesh(nxy, bn, interval)
    if length(nxy) == 2
        nx = nxy[1];
        ny = nxy[2];
    else
        nx = nxy;
        ny = nx;
    end
    if length(interval) == 2
        interval = [interval[1], interval[2], interval[1], interval[2]];
    end
    
    # mesh
    Nv = nx*ny;                 # number of vertex nodes
    xv = zeros(2,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = (nx-1)*(ny-1);        # number of elements
    el = zeros(Int, 4, nel);    # element vertex maps
    etypes = 3*ones(Int, nel);  # element types (gmsh number)
    numvert = 4*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    Nf = nx*(ny-1) + ny*(nx-1); # number of faces
    f2n = zeros(Int, 2,Nf);     # face2node mapping
    f2e = zeros(Int, 2,Nf);     # face2element mapping
    e2f = zeros(Int, 4,nel);    # element2face mapping
    normals = ones(2,Nf);       # normals of faces
    bdryID = zeros(Int, Nf);    # BID of faces (0=interior)
    
    if bn == 4
        bids = [1,2,3,4]; # x=0, x=1, y=0, y=1
        allbids = [1,2,3,4];
    elseif bn == 3
        bids = [1,2,3]; # x=0 , x=1, y=0,1
        allbids = [1,2,3,3];
    elseif bn == 2
        bids = [1,2]; # x=0,1 , y=0,1
        allbids = [1,1,2,2];
    else
        bids = [1]; # everywhere
        allbids = [1,1,1,1];
    end
    
    scalex = interval[2]-interval[1];
    hx = scalex/(nx-1);         # uniformly divided
    scaley = interval[4]-interval[3];
    hy = scaley/(ny-1);         # uniformly divided
    
    # vertex nodes are ordered lexicographically
    for j=1:ny
        for i=1:nx
            k = i + (j-1)*nx;
            xv[1, k] = interval[1] + (i-1)*hx;
            xv[2, k] = interval[3] + (j-1)*hy;
        end
    end
    
    # Elements are ordered the same way
    for j=1:(ny-1)
        for i=1:(nx-1)
            ei = i + (j-1)*(nx-1); # element index
            
            el[1,ei] = i + (j-1)*nx;
            el[2,ei] = i + (j-1)*nx + 1;
            el[3,ei] = i + (j)*nx;
            el[4,ei] = i + (j)*nx + 1;
            
            # face2node, face2element, element2face, normals
            f1 = 2*(ei-1)+1 + (j-1); # left face index
            f2 = f1+1; # bottom
            f3 = f1+2; # right
            f4 = j==(ny-1) ? Nf-(nx-1)+i : j*2*(nx-1) + j + i*2; # top
            
            f2n[1,f1] = el[1,ei];
            f2n[2,f1] = el[3,ei];
            f2n[1,f2] = el[1,ei];
            f2n[2,f2] = el[2,ei];
            f2n[1,f3] = el[2,ei];
            f2n[2,f3] = el[4,ei];
            f2n[1,f4] = el[3,ei];
            f2n[2,f4] = el[4,ei];
            
            f2e[2,f1] = ei;
            f2e[2,f2] = ei;
            f2e[1,f3] = ei;
            f2e[1,f4] = ei;
            
            e2f[1,ei] = f1;
            e2f[2,ei] = f2;
            e2f[3,ei] = f3;
            e2f[4,ei] = f4;
            
            normals[:,f1] = [1,0];
            normals[:,f2] = [0,1];
            normals[:,f3] = [1,0];
            normals[:,f4] = [0,1];
        end
    end
    
    # boundaries
    for j=1:(ny-1)
        eleft = (j-1)*(nx-1) + 1;
        eright = j*(nx-1);
        
        fleft = e2f[1,eleft];
        fright = e2f[3,eright];
        
        bdryID[fleft] = 1; # always 1
        bdryID[fright] = allbids[2];
        
        # need to change normals, f2e for left side
        normals[:,fleft] = [-1,0];
        f2e[:,fleft] = [eleft, 0];
    end
    for i=1:(nx-1)
        ebottom = i;
        etop = nel - (nx-1) + i;
        
        fbottom = e2f[2,ebottom];
        ftop = e2f[4,etop];
        
        bdryID[fbottom] = allbids[3];
        bdryID[ftop] = allbids[4];
        
        # need to change normals, f2e for bottom side
        normals[:,fbottom] = [0,-1];
        f2e[:,fbottom] = [ebottom, 0];
    end
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    # grid
    ord = config.basis_order_min;
    refel = build_refel(2, ord, 4, config.elemental_nodes);  # refel for elements
    tmprefel = build_refel(1, ord, 2, config.elemental_nodes); # 1D refel for getting face nodes
    leftnodes =    map_face_nodes_2d(tmprefel.g, [-1,-1], [-1,1]); # maps 1D gauss nodes to 2D face
    rightnodes =   map_face_nodes_2d(tmprefel.g, [1,-1], [1,1]); # maps 1D gauss nodes to 2D face
    bottomnodes =  map_face_nodes_2d(tmprefel.g, [-1,-1], [1,-1]); # maps 1D gauss nodes to 2D face
    topnodes =     map_face_nodes_2d(tmprefel.g, [-1,1], [1,1]); # maps 1D gauss nodes to 2D face
    frefelLeft =   custom_quadrature_refel(refel, leftnodes, tmprefel.wg); # refel for left face
    frefelRight =  custom_quadrature_refel(refel, rightnodes, tmprefel.wg); # refel for right face
    frefelBottom = custom_quadrature_refel(refel, bottomnodes, tmprefel.wg); # refel for bottom face
    frefelTop =    custom_quadrature_refel(refel, topnodes, tmprefel.wg); # refel for top face
    refelfc = [frefelLeft, frefelRight, frefelBottom, frefelTop];
    
    N = ((nx-1)*ord + 1)*((ny-1)*ord + 1); # number of total nodes
    Np = refel.Np;              # number of nodes per element
    x = zeros(2,N);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryfc = [];                # index of faces touching each BID
    bdrynorm = [];              # normal at boundary nodes
    loc2glb = zeros(Int, Np, nel)# local to global index map for each element's nodes
    glbvertex = zeros(Int, 4, nel);# local to global for vertices
    f2glb = zeros(Int, size(leftnodes,2), Nf);  # face node local to global
    fvtx2glb = zeros(Int, 2, Nf);# face vertex local to global
    
    # Elements are ordered along x then y
    n1d = ord+1;
    rowsize = (nx-1)*(n1d-1) + 1;
    colsize = (ny-1)*(n1d-1) + 1;
    for j=1:(ny-1)
        for i=1:(nx-1)
            ei = i + (j-1)*(nx-1); # element index
            glbvertex[1,ei] = ((j-1)*(n1d-1))*rowsize + (i-1)*(n1d-1) + 1;
            glbvertex[2,ei] = ((j-1)*(n1d-1))*rowsize + i*(n1d-1) + 1;
            glbvertex[3,ei] = (j*(n1d-1))*rowsize + (i-1)*(n1d-1) + 1;
            glbvertex[4,ei] = (j*(n1d-1))*rowsize + i*(n1d-1) + 1;
            
            x1 = [interval[1] + (i-1)*hx ; interval[3] + (j-1)*hy]; # southwest corner
            for jj=1:n1d-1
                row = (j-1)*(n1d-1) + jj;
                for ii=1:n1d-1
                    gi = (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                    li = (jj-1)*n1d + ii; # local index of this node
                    x[1,gi] = x1[1] .+ hx*0.5 .* (refel.r[li,1] + 1); # coordinates of this node
                    x[2,gi] = x1[2] .+ hy*0.5 .* (refel.r[li,2] + 1);
                    loc2glb[li, ei] = gi; # local to global map
                end
                loc2glb[jj*n1d, ei] = (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
            end
            row = j*(n1d-1)+1;
            for ii=1:n1d
                loc2glb[n1d*(n1d-1)+ii, ei] = ((row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
            end
        end
        # The last node on the row was not set
        x1 = [interval[1] + (nx-2)*hx ; interval[3] + (j-1)*hy]; # southwest corner of last element
        for jj=1:n1d-1
            li = jj*n1d; # local index of this node
            x[1, ((j-1)*(n1d-1)+jj)*rowsize] = x1[1] + hx*0.5 * (refel.r[li,1] + 1); # coordinates of this node
            x[2, ((j-1)*(n1d-1)+jj)*rowsize] = x1[2] + hy*0.5 * (refel.r[li,2] + 1);
        end
    end
    # The top row was not set
    rowstart = ((n1d-1)*(ny-1))*rowsize;
    for i=1:(nx-1)
        x1 = [interval[1] + (i-1)*hx ; interval[3] + (ny-2)*hy]; # southwest corner of top element
        for ii=1:n1d-1
            gi = rowstart + (i-1)*(n1d-1) + ii;
            li = (n1d-1)*n1d + ii; # local index of this node
            x[1, gi] = x1[1] + hx*0.5 * (refel.r[li,1] + 1); # coordinates of this node
            x[2, gi] = x1[2] + hy*0.5 * (refel.r[li,2] + 1);
        end
    end
    # corner
    x[:,N] = [interval[2]; interval[4]];
    
    # face global maps
    # It's easier to loop over elements than faces, though less efficient
    for ei=1:nel
        f1 = e2f[1,ei];
        f2 = e2f[2,ei];
        f3 = e2f[3,ei];
        f4 = e2f[4,ei];
        
        fvtx2glb[1,f1] = glbvertex[1,ei];
        fvtx2glb[2,f1] = glbvertex[3,ei];
        fvtx2glb[1,f2] = glbvertex[1,ei];
        fvtx2glb[2,f2] = glbvertex[2,ei];
        fvtx2glb[1,f3] = glbvertex[2,ei];
        fvtx2glb[2,f3] = glbvertex[4,ei];
        fvtx2glb[1,f4] = glbvertex[3,ei];
        fvtx2glb[2,f4] = glbvertex[4,ei];
        
        f2glb[:, f1] = loc2glb[1:n1d:(Np-n1d+1), ei];
        f2glb[:, f2] = loc2glb[1:n1d, ei];
        f2glb[:, f3] = loc2glb[n1d:n1d:Np, ei];
        f2glb[:, f4] = loc2glb[(Np-n1d+1):Np, ei];
    end
    
    # boundaries
    if bn == 4
        # boundary nodes and norms
        bdry1 = zeros(Int, colsize); # x=0
        bdry2 = zeros(Int, colsize); # x=1
        bdry3 = zeros(Int, rowsize - 2); # y=0
        bdry4 = zeros(Int, rowsize - 2); # y=1
        bdrynorm1 = zeros(2, length(bdry1));
        bdrynorm2 = zeros(2, length(bdry2));
        bdrynorm3 = zeros(2, length(bdry3));
        bdrynorm4 = zeros(2, length(bdry4));
        for i=2:rowsize-1
            bdry3[i-1] = i; # bottom
            bdrynorm3[:, i-1] = [0,-1];
            bdry4[i-1] = N-rowsize + i; # top
            bdrynorm4[:, i-1] = [0,1];
        end
        for i=1:colsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[:, i] = [-1,0];
            bdry2[i] = i*rowsize; # right
            bdrynorm2[:, i] = [1,0];
        end
        
        bdry = [bdry1, bdry2, bdry3, bdry4];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3, bdrynorm4];
        
        # boundary faces
        bdryfc1 = zeros(Int, (ny-1));
        bdryfc2 = zeros(Int, (ny-1));
        bdryfc3 = zeros(Int, (nx-1));
        bdryfc4 = zeros(Int, (nx-1));
        for i=1:(nx-1)
            bdryfc3[i] = i*2; # bottom
            bdryfc4[i] = Nf-(nx-1)+i; # top
        end
        for i=1:(ny-1)
            bdryfc1[i] = 2*(nx-1)*(i-1) + i; # left
            bdryfc2[i] = 2*(nx-1)*i + i; # right
        end
        bdryfc = [bdryfc1, bdryfc2, bdryfc3, bdryfc4];
        
    elseif bn == 3
        # boundary nodes and norms
        bdry1 = zeros(Int, colsize); # x=0
        bdry2 = zeros(Int, colsize); # x=1
        bdry3 = zeros(Int, rowsize*2 - 4); # y=0,1
        bdrynorm1 = zeros(2, length(bdry1));
        bdrynorm2 = zeros(2, length(bdry2));
        bdrynorm3 = zeros(2, length(bdry3));
        for i=2:rowsize-1
            bdry3[i-1] = i; # bottom
            bdrynorm3[:, i-1] = [0,-1];
            bdry3[i-3+rowsize] = N-rowsize + i; # top
            bdrynorm3[:, i-3+rowsize] = [0,1];
        end
        for i=1:colsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[:, i] = [-1,0];
            bdry2[i] = i*rowsize; # right
            bdrynorm2[:, i] = [1,0];
        end
        
        bdry = [bdry1, bdry2, bdry3];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3];
        
        # boundary faces
        bdryfc1 = zeros(Int, (ny-1));
        bdryfc2 = zeros(Int, (ny-1));
        bdryfc3 = zeros(Int, (nx-1)*2);
        for i=1:(nx-1)
            bdryfc3[i] = i*2; # bottom
            bdryfc3[i+nx-1] = Nf-(nx-1)+i; # top
        end
        for i=1:(ny-1)
            bdryfc1[i] = 2*(nx-1)*(i-1) + i; # left
            bdryfc2[i] = 2*(nx-1)*i + i; # right
        end
        bdryfc = [bdryfc1, bdryfc2, bdryfc3];

    elseif bn == 2
        # boundary nodes and norms
        bdry1 = zeros(Int, colsize*2); # x=0,1
        bdry2 = zeros(Int, rowsize*2 - 4); # y=0,1
        bdrynorm1 = zeros(2, length(bdry1));
        bdrynorm2 = zeros(2, length(bdry2));
        for i=2:rowsize-1
            bdry2[i-1] = i; # bottom
            bdrynorm2[:, i-1] = [0,-1];
            bdry2[i-3+rowsize] = N-rowsize + i; # top
            bdrynorm2[:, i-3+rowsize] = [0,1];
        end
        for i=1:colsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[:, i] = [-1,0];
            bdry1[i + colsize] = i*rowsize; # right
            bdrynorm1[:, i + colsize] = [1,0];
        end
        
        bdry = [bdry1, bdry2];
        bdrynorm = [bdrynorm1, bdrynorm2];
        
        # boundary faces
        bdryfc1 = zeros(Int, (ny-1)*2);
        bdryfc2 = zeros(Int, (nx-1)*2);
        for i=1:(nx-1)
            bdryfc2[i] = i*2; # bottom
            bdryfc2[i+nx-1] = Nf-(nx-1)+i; # top
        end
        for i=1:(ny-1)
            bdryfc1[i] = 2*(nx-1)*(i-1) + i; # left
            bdryfc1[i + ny-1] = 2*(nx-1)*i + i; # right
        end
        bdryfc = [bdryfc1, bdryfc2];
        
    else # bn=1
        # boundary nodes and norms
        bdry = zeros(Int, rowsize*2 + colsize*2 - 4);
        bdrynorm1 = zeros(2, length(bdry));
        for i=1:rowsize
            bdry[i] = i; # bottom
            bdrynorm1[:, i] = [0,-1];
            bdry[i+rowsize] = N-rowsize + i; # top
            bdrynorm1[:, i+rowsize] = [0,1];
        end
        for i=2:colsize-1
            bdry[i-1 + rowsize*2] = (i-1)*rowsize + 1; # left
            bdrynorm1[:, i-1 + rowsize*2] = [-1,0];
            bdry[i-3 + rowsize*2 + colsize] = i*rowsize; # right
            bdrynorm1[:, i-3 + rowsize*2 + colsize] = [1,0];
        end
        
        bdry = [bdry];
        bdrynorm = [bdrynorm1];
        
        # boundary faces
        bdryfc = zeros(Int, (nx-1)*2 + (ny-1)*2);
        for i=1:(nx-1)
            bdryfc[i] = i*2; # bottom
            bdryfc[i+nx-1] = Nf-(nx-1)+i; # top
        end
        for i=1:(ny-1)
            bdryfc[(nx-1)*2 + i] = 2*(nx-1)*(i-1) + i; # left
            bdryfc[(nx-1)*2 + i + ny-1] = 2*(nx-1)*i + i; # right
        end
        bdryfc = [bdryfc];
    end
    
    grid = Grid(x, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb);
    
    return (mesh, refel, refelfc, grid);
end

#=
# Builds a 3D hex mesh
# nx = number of vertices
# bn = number of boundary regions
# interval = limits of the square domain
=#
function simple_hex_mesh(nxyz, bn, interval)
    if length(nxyz) == 3
        nx = nxyz[1];
        ny = nxyz[2];
        nz = nxyz[3]
    else
        nx = nxyz;
        ny = nx;
        nz = nx;
    end
    if length(interval) == 2
        interval = [interval[1], interval[2], interval[1], interval[2], interval[1], interval[2]];
    end
    
    # mesh
    Nv = nx*ny*nz;              # number of vertex nodes
    xv = zeros(3,Nv);           # coordinates of vertices
    ind = Array(1:Nv);          # indices are in order
    nel = (nx-1)*(ny-1)*(nz-1); # number of elements
    el = zeros(Int, 8, nel);    # element vertex maps
    etypes = 5*ones(Int, nel);  # element types (gmsh number)
    numvert = 8*ones(Int, nel); # number of vertices per element
    invind = invert_index(ind); # inverse index ind[i] = j -> invind[j] = i
    Nf = nz*(nx-1)*(ny-1) + ny*(nx-1)*(nz-1) + nx*(nz-1)*(ny-1); # number of faces
    f2n = zeros(Int, 4,Nf);     # face2node mapping
    f2e = zeros(Int, 2,Nf);     # face2element mapping
    e2f = zeros(Int, 6,nel);    # element2face mapping
    normals = ones(3,Nf);       # normals of faces
    bdryID = zeros(Int, Nf);    # BID of faces (0=interior)
    
    if bn == 6
        bids = [1,2,3,4,5,6]; # all separate
        allbids = [1,2,3,4,5,6];
    elseif bn == 5
        bids = [1,2,3,4,5]; # combine z
        allbids = [1,2,3,4,5,5];
    elseif bn == 4
        bids = [1,2,3,4]; # combine y and z
        allbids = [1,2,3,3,4,4];
    elseif bn == 3
        bids = [1,2,3]; # combine x,y,z
        allbids = [1,1,2,2,3,3];
    elseif bn == 2
        bids = [1,2]; # x=0, other
        allbids = [1,2,2,2,2,2];
    else
        bids = [1]; # everywhere
        allbids = [1,1,1,1,1,1];
    end
    
    scalex = interval[2]-interval[1];
    hx = scalex/(nx-1);         # uniformly divided
    scaley = interval[4]-interval[3];
    hy = scaley/(ny-1);         # uniformly divided
    scalez = interval[6]-interval[5];
    hz = scalez/(nz-1);         # uniformly divided
    
    # vertex nodes are ordered lexicographically
    for k=1:nz
        for j=1:ny
            for i=1:nx
                ni = i + (j-1)*nx;
                xv[1, ni] = interval[1] + (i-1)*hx;
                xv[2, ni] = interval[3] + (j-1)*hy;
                xv[3, ni] = interval[5] + (k-1)*hz;
            end
        end
    end
    
    # Elements are ordered the same way
    nextface = 1;
    for k=1:(nz-1)
        for j=1:(ny-1)
            for i=1:(nx-1)
                ei = i + (j-1)*(nx-1) + (k-1)*(ny-1)*(nx-1); # element index
                
                el[1,ei] = i + (j-1)*nx + (k-1)*nx*nx;
                el[2,ei] = i + (j-1)*nx + (k-1)*nx*nx + 1;
                el[3,ei] = i + (j)*nx + (k-1)*nx*nx;
                el[4,ei] = i + (j)*nx + (k-1)*nx*nx + 1;
                el[5,ei] = i + (j-1)*nx + (k)*nx*nx;
                el[6,ei] = i + (j-1)*nx + (k)*nx*nx + 1;
                el[7,ei] = i + (j)*nx + (k)*nx*nx;
                el[8,ei] = i + (j)*nx + (k)*nx*nx + 1;
                
                # face2node, face2element, element2face, normals
                if i==1
                    f1 = nextface;
                    nextface = nextface+1;
                else
                    f1 = e2f[2,ei-1];# left face index
                end
                f2 = nextface; # right
                nextface = nextface+1;
                
                if j==1
                    f3 = nextface;
                    nextface = nextface+1;
                else
                    f3 = e2f[4,ei-(nx-1)]; # bottom
                end
                f4 = nextface; # top
                nextface = nextface+1;
                
                if k==1
                    f5 = nextface;
                    nextface = nextface+1;
                else
                    f5 = e2f[6,ei-((nx-1)*(ny-1))]; # front
                end
                f6 = nextface; # back
                nextface = nextface+1;
                
                f2n[:,f1] = el[[1,3,5,7],ei];
                f2n[:,f2] = el[[2,4,6,8],ei];
                f2n[:,f3] = el[[1,2,5,6],ei];
                f2n[:,f4] = el[[3,4,7,8],ei];
                f2n[:,f5] = el[[1,2,3,4],ei];
                f2n[:,f6] = el[[5,6,7,8],ei];
                
                f2e[2,f1] = ei;
                f2e[1,f2] = ei;
                f2e[2,f3] = ei;
                f2e[1,f4] = ei;
                f2e[2,f5] = ei;
                f2e[1,f6] = ei;
                
                e2f[1,ei] = f1;
                e2f[2,ei] = f2;
                e2f[3,ei] = f3;
                e2f[4,ei] = f4;
                e2f[5,ei] = f5;
                e2f[6,ei] = f6;
                
                normals[:,f1] = [1,0,0];
                normals[:,f2] = [1,0,0];
                normals[:,f3] = [0,1,0];
                normals[:,f4] = [0,1,0];
                normals[:,f5] = [0,0,1];
                normals[:,f6] = [0,0,1];
            end
        end
    end
    
    # boundaries
    for k=1:(nz-1)
        for j=1:(ny-1)
            eleft = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
            eright = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
            
            fleft = e2f[1,eleft];
            fright = e2f[2,eright];
            
            bdryID[fleft] = 1; # always 1
            bdryID[fright] = allbids[2];
            
            # need to change normals, f2e for left side
            normals[:,fleft] = [-1,0,0];
            f2e[:,fleft] = [eleft, 0];
        end
    end
    
    for k=1:(nz-1)
        for i=1:(nx-1)
            ebottom = (k-1)*(nx-1)*(ny-1) + i;
            etop = (k-1)*(nx-1)*(ny-1) + (nx-1)*(ny-2) + i;
            
            fbottom = e2f[3,ebottom];
            ftop = e2f[4,etop];
            
            bdryID[fbottom] = allbids[3];
            bdryID[ftop] = allbids[4];
            
            # need to change normals, f2e for bottom side
            normals[:,fbottom] = [0,-1, 0];
            f2e[:,fbottom] = [ebottom, 0];
        end
    end
    
    for j=1:(ny-1)
        for i=1:(nx-1)
            efront = (j-1)*(nx-1) + i;
            eback = (nz-2)*(nx-1)*(ny-1) + efront;
            
            ffront = e2f[5,efront];
            fback = e2f[6,eback];
            
            bdryID[ffront] = allbids[5];
            bdryID[fback] = allbids[6];
            
            # need to change normals, f2e for bottom side
            normals[:,ffront] = [0,0,-1];
            f2e[:,ffront] = [efront, 0];
        end
    end
    
    
    mesh = MeshData(Nv, xv, ind, nel, el, etypes, numvert, invind, f2n, f2e, e2f, normals, bdryID); # MeshData struct
    
    
    #################################################################################################
    
    ord = config.basis_order_min;
    refel = build_refel(3, ord, 6, config.elemental_nodes);
    tmprefel = build_refel(2, ord, 4, config.elemental_nodes); # 2D refel for getting face nodes
    # leftnodes =    map_face_nodes_3d(tmprefel.g, [-1,-1,-1], [-1,1,-1], [-1,-1,1], [-1,1,1]); # maps 2D gauss nodes to 3D face
    # rightnodes =   map_face_nodes_3d(tmprefel.g, [1,-1,-1], [1,1,-1], [1,-1,1], [1,1,1]); 
    # bottomnodes =  map_face_nodes_3d(tmprefel.g, [-1,-1,-1], [1,-1,-1], [-1,-1,1], [1,-1,1]); 
    # topnodes =     map_face_nodes_3d(tmprefel.g, [-1,1,-1], [1,1,-1], [-1,1,1], [1,1,1]); 
    # frontnodes =   map_face_nodes_3d(tmprefel.g, [-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1]); 
    # backnodes =    map_face_nodes_3d(tmprefel.g, [-1,-1,1], [1,-1,1], [-1,1,1], [1,1,1]); 
    g2d= tmprefel.g;
    leftnodes =    -ones(3, length(tmprefel.wg));
    rightnodes =   ones(3, length(tmprefel.wg));
    bottomnodes =  -ones(3, length(tmprefel.wg));
    topnodes =     ones(3, length(tmprefel.wg));
    frontnodes =   -ones(3, length(tmprefel.wg));
    backnodes =    ones(3, length(tmprefel.wg));
    
    leftnodes[[2,3],:] = g2d';
    rightnodes[[2,3],:] = g2d';
    bottomnodes[[1,3],:] = g2d';
    topnodes[[1,3],:] = g2d';
    frontnodes[[1,2],:] = g2d';
    backnodes[[1,2],:] = g2d';
    
    frefelLeft =   custom_quadrature_refel(refel, leftnodes, tmprefel.wg); # refel for left face
    frefelRight =  custom_quadrature_refel(refel, rightnodes, tmprefel.wg); # refel for right face
    frefelBottom = custom_quadrature_refel(refel, bottomnodes, tmprefel.wg); # refel for bottom face
    frefelTop =    custom_quadrature_refel(refel, topnodes, tmprefel.wg); # refel for top face
    frefelFront =  custom_quadrature_refel(refel, frontnodes, tmprefel.wg); # refel for front face
    frefelBack =   custom_quadrature_refel(refel, backnodes, tmprefel.wg); # refel for back face
    refelfc = [frefelLeft, frefelRight, frefelBottom, frefelTop, frefelFront, frefelBack];
    
    N = ((nx-1)*ord + 1)*((ny-1)*ord + 1)*((nz-1)*ord + 1);# number of total nodes
    Np = refel.Np;              # number of nodes per element
    x = zeros(3,N);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryfc = [];                # index of elements touching each BID
    bdrynorm = [];              # normal at boundary nodes
    loc2glb = zeros(Int, Np, nel)# local to global index map for each element's nodes
    glbvertex = zeros(Int, 8, nel);# local to global for vertices
    f2glb = zeros(Int, (ord+1)*(ord+1), Nf);  # face node local to global
    fvtx2glb = zeros(Int, 4, Nf);# face vertex local to global
    
    # Start with a 2d quad mesh
    (mesh2d, refel2d, refelfc2d, grid2d) = simple_quad_mesh([nx,ny], 1, interval[1:4]);
    (mesh1d, refel1d, refelfc1d, grid1d) = simple_line_mesh(nz, 1, interval[5:6]);

    # Use 2d quad mesh to build hex mesh nodes
    n1d = ord+1;
    rowsize = (nx-1)*(n1d-1) + 1;
    slicesize = ((nx-1)*(n1d-1) + 1)*((ny-1)*(n1d-1) + 1);
    zvals = grid1d.allnodes[:];
    for k=1:length(zvals);
        range = ((k-1)*slicesize+1):(k*slicesize);
        x[1:2, range] = grid2d.allnodes;
        x[3, range] = zvals[k].*ones(1,slicesize);
    end

    # Elements are ordered lexicographically
    for k=1:(nz-1)
        for j=1:(ny-1) # element indices
            for i=1:(nx-1)
                ei = i + (j-1)*(nx-1) + (k-1)*(nx-1)*(ny-1); # element index
                glbvertex[1,ei] = (k-1)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[2,ei] = (k-1)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[3,ei] = (k-1)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[4,ei] = (k-1)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[5,ei] = (k)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[6,ei] = (k)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[7,ei] = (k)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[8,ei] = (k)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                
                for kk=1:n1d-1
                    slice = (k-1)*(n1d-1) + kk;
                    for jj=1:n1d-1 # elemental node indices
                        row = (j-1)*(n1d-1) + jj;
                        for ii=1:n1d-1
                            gi = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                            li = (kk-1)*n1d*n1d + (jj-1)*n1d + ii; # local index of this node
                            loc2glb[li,ei] = gi; # local to global map
                        end
                        loc2glb[(kk-1)*n1d*n1d + jj*n1d, ei] = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
                    end

                    row = j*(n1d-1)+1;
                    for ii=1:n1d
                        loc2glb[(kk-1)*n1d*n1d + n1d*(n1d-1)+ii, ei] = ((slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
                    end
                end

                slice = k*(n1d-1)+1;  # last in local z
                for jj=1:n1d-1
                    row = (j-1)*(n1d-1) + jj;
                    for ii=1:n1d-1
                        gi = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                        li = (n1d-1)*n1d*n1d + (jj-1)*n1d + ii; # local index of this node
                        loc2glb[li, ei] = gi; # local to global map
                    end
                    loc2glb[(n1d-1)*n1d*n1d + jj*n1d, ei] = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
                end
                row = j*(n1d-1)+1;
                for ii=1:n1d
                    loc2glb[(n1d-1)*n1d*n1d + n1d*(n1d-1)+ii, ei] = ((slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
                end
            end
        end
    end
    
    # face global maps
    # It's easier to loop over elements than faces, though less efficient
    for ei=1:nel
        f1 = e2f[1,ei];
        f2 = e2f[2,ei];
        f3 = e2f[3,ei];
        f4 = e2f[4,ei];
        f5 = e2f[5,ei];
        f6 = e2f[6,ei];
        
        fvtx2glb[1,f1] = glbvertex[1,ei];
        fvtx2glb[2,f1] = glbvertex[3,ei];
        fvtx2glb[3,f1] = glbvertex[5,ei];
        fvtx2glb[4,f1] = glbvertex[7,ei];
        
        fvtx2glb[1,f2] = glbvertex[2,ei];
        fvtx2glb[2,f2] = glbvertex[4,ei];
        fvtx2glb[3,f2] = glbvertex[6,ei];
        fvtx2glb[4,f2] = glbvertex[8,ei];
        
        fvtx2glb[1,f3] = glbvertex[1,ei];
        fvtx2glb[2,f3] = glbvertex[2,ei];
        fvtx2glb[3,f3] = glbvertex[5,ei];
        fvtx2glb[4,f3] = glbvertex[6,ei];
        
        fvtx2glb[1,f4] = glbvertex[3,ei];
        fvtx2glb[2,f4] = glbvertex[4,ei];
        fvtx2glb[3,f4] = glbvertex[7,ei];
        fvtx2glb[4,f4] = glbvertex[8,ei];
        
        fvtx2glb[1,f5] = glbvertex[1,ei];
        fvtx2glb[2,f5] = glbvertex[2,ei];
        fvtx2glb[3,f5] = glbvertex[3,ei];
        fvtx2glb[4,f5] = glbvertex[4,ei];
        
        fvtx2glb[1,f6] = glbvertex[5,ei];
        fvtx2glb[2,f6] = glbvertex[6,ei];
        fvtx2glb[3,f6] = glbvertex[7,ei];
        fvtx2glb[4,f6] = glbvertex[8,ei];
        
        f2glb[:, f1] = loc2glb[1:n1d:Np, ei];
        f2glb[:, f2] = loc2glb[n1d:n1d:Np, ei];
        for k=1:n1d
            f2glb[((k-1)*n1d+1):(k*n1d), f3] = loc2glb[((k-1)*n1d*n1d).+(1:n1d), ei];
            f2glb[((k-1)*n1d+1):(k*n1d), f4] = loc2glb[(k*n1d*n1d-n1d).+(1:n1d), ei];
        end
        
        f2glb[:, f5] = loc2glb[1:(n1d*n1d), ei];
        f2glb[:, f6] = loc2glb[(Np-n1d*n1d+1):Np, ei];
    end

    # boundaries
    xlen = (nx-1)*(n1d-1) + 1;
    ylen = (ny-1)*(n1d-1) + 1;
    zlen = (nz-1)*(n1d-1) + 1;
    xbdrysize = zlen*ylen;
    ybdrysize = (xlen-2)*zlen;
    zbdrysize = (xlen-2)*(ylen-2);
    zslice = xlen*ylen;
    if bn == 6 # all separate
        bdry1 = zeros(Int, xbdrysize); # x=0
        bdry2 = zeros(Int, xbdrysize); # x=1
        bdry3 = zeros(Int, ybdrysize); # y=0
        bdry4 = zeros(Int, ybdrysize); # y=1
        bdry5 = zeros(Int, zbdrysize); # z=0
        bdry6 = zeros(Int, zbdrysize); # z=1
        
        bdrynorm1 = zeros(3, xbdrysize); # x=0
        bdrynorm2 = zeros(3, xbdrysize); # x=1
        bdrynorm3 = zeros(3, ybdrysize); # y=0
        bdrynorm4 = zeros(3, ybdrysize); # y=1
        bdrynorm5 = zeros(3, zbdrysize); # z=0
        bdrynorm6 = zeros(3, zbdrysize); # z=1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry2[ind1] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm2[:,ind1] = [1,0,0];
                ind1 = ind1+1;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry3[ind2] = (k-1)*zslice + i; # y=0
                bdry4[ind2] = k*zslice - xlen + i; # y = 1
                bdrynorm3[:,ind2] = [0,-1,0];
                bdrynorm4[:,ind2] = [0,1,0];
                ind2 = ind2+1;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry5[ind3] = (j-1)*xlen + i; # z=0
                bdry6[ind3] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm5[:,ind3] = [0,0,-1];
                bdrynorm6[:,ind3] = [0,0,1];
                ind3 = ind3+1;
            end
        end
        
        bdry = [bdry1, bdry2, bdry3, bdry4, bdry5, bdry6];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3, bdrynorm4, bdrynorm5, bdrynorm6];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1));
        bdryfc2 = zeros(Int, (ny-1)*(nz-1));
        bdryfc3 = zeros(Int, (nx-1)*(nz-1));
        bdryfc4 = zeros(Int, (nx-1)*(nz-1));
        bdryfc5 = zeros(Int, (nx-1)*(ny-1));
        bdryfc6 = zeros(Int, (nx-1)*(ny-1));
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc2[ind1] = e2f[2,e2];
                ind1 = ind1+1;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc3[ind2] = e2f[3,e1];
                bdryfc4[ind2] = e2f[4,e2];
                ind2 = ind2+1;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc5[ind3] = e2f[5,e1];
                bdryfc6[ind3] = e2f[6,e2];
                ind3 = ind3+1;
            end
        end
        
        bdryfc = [bdryfc1, bdryfc2, bdryfc3, bdryfc4, bdryfc5, bdryfc6];
        
    elseif bn == 5 # combine zs
        bdry1 = zeros(Int, xbdrysize); # x=0
        bdry2 = zeros(Int, xbdrysize); # x=1
        bdry3 = zeros(Int, ybdrysize); # y=0
        bdry4 = zeros(Int, ybdrysize); # y=1
        bdry5 = zeros(Int, zbdrysize*2); # z=0,1
        
        bdrynorm1 = zeros(3, xbdrysize); # x=0
        bdrynorm2 = zeros(3, xbdrysize); # x=1
        bdrynorm3 = zeros(3, ybdrysize); # y=0
        bdrynorm4 = zeros(3, ybdrysize); # y=1
        bdrynorm5 = zeros(3, zbdrysize*2); # z=0,1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry2[ind1] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm2[:,ind1] = [1,0,0];
                ind1 = ind1+1;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry3[ind2] = (k-1)*zslice + i; # y=0
                bdry4[ind2] = k*zslice - xlen + i; # y = 1
                bdrynorm3[:,ind2] = [0,-1,0];
                bdrynorm4[:,ind2] = [0,1,0];
                ind2 = ind2+1;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry5[ind3] = (j-1)*xlen + i; # z=0
                bdry5[ind3+1] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm5[:,ind3] = [0,0,-1];
                bdrynorm5[:,ind3+1] = [0,0,1];
                ind3 = ind3+2;
            end
        end
        
        bdry = [bdry1, bdry2, bdry3, bdry4, bdry5];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3, bdrynorm4, bdrynorm5];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1));
        bdryfc2 = zeros(Int, (ny-1)*(nz-1));
        bdryfc3 = zeros(Int, (nx-1)*(nz-1));
        bdryfc4 = zeros(Int, (nx-1)*(nz-1));
        bdryfc5 = zeros(Int, (nx-1)*(ny-1)*2);
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc2[ind1] = e2f[2,e2];
                ind1 = ind1+1;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc3[ind2] = e2f[3,e1];
                bdryfc4[ind2] = e2f[4,e2];
                ind2 = ind2+1;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc5[ind3] = e2f[5,e1];
                bdryfc5[ind3+1] = e2f[6,e2];
                ind3 = ind3+2;
            end
        end
        
        bdryfc = [bdryfc1, bdryfc2, bdryfc3, bdryfc4, bdryfc5];
        
    elseif bn == 4 # combine ys and zs
        bdry1 = zeros(Int, xbdrysize); # x=0
        bdry2 = zeros(Int, xbdrysize); # x=1
        bdry3 = zeros(Int, ybdrysize*2); # y=0,1
        bdry4 = zeros(Int, zbdrysize*2); # z=0,1
        
        bdrynorm1 = zeros(3, xbdrysize); # x=0
        bdrynorm2 = zeros(3, xbdrysize); # x=1
        bdrynorm3 = zeros(3, ybdrysize*2); # y=0,1
        bdrynorm4 = zeros(3, zbdrysize*2); # z=0,1
        
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry2[ind1] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm2[:,ind1] = [1,0,0];
                ind1 = ind1+1;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry3[ind2] = (k-1)*zslice + i; # y=0
                bdry3[ind2+1] = k*zslice - xlen + i; # y = 1
                bdrynorm3[:,ind2] = [0,-1,0];
                bdrynorm3[:,ind2+1] = [0,1,0];
                ind2 = ind2+2;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry4[ind3] = (j-1)*xlen + i; # z=0
                bdry4[ind3+1] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm4[:,ind3] = [0,0,-1];
                bdrynorm4[:,ind3+1] = [0,0,1];
                ind3 = ind3+2;
            end
        end
        
        bdry = [bdry1, bdry2, bdry3, bdry4];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3, bdrynorm4];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1));
        bdryfc2 = zeros(Int, (ny-1)*(nz-1));
        bdryfc3 = zeros(Int, (nx-1)*(nz-1)*2);
        bdryfc4 = zeros(Int, (nx-1)*(ny-1)*2);
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc2[ind1] = e2f[2,e2];
                ind1 = ind1+1;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc3[ind2] = e2f[3,e1];
                bdryfc3[ind2+1] = e2f[4,e2];
                ind2 = ind2+2;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc4[ind3] = e2f[5,e1];
                bdryfc4[ind3+1] = e2f[6,e2];
                ind3 = ind3+2;
            end
        end
        
        bdryfc = [bdryfc1, bdryfc2, bdryfc3, bdryfc4];
        
    elseif bn == 3 # xs, ys, zs
        bdry1 = zeros(Int, xbdrysize*2); # x=0,1
        bdry2 = zeros(Int, ybdrysize*2); # y=0,1
        bdry3 = zeros(Int, zbdrysize*2); # z=0,1
        
        bdrynorm1 = zeros(3, xbdrysize*2); # x=0,1
        bdrynorm2 = zeros(3, ybdrysize*2); # y=0,1
        bdrynorm3 = zeros(3, zbdrysize*2); # z=0,1
        
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry1[ind1+1] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm1[:,ind1+1] = [1,0,0];
                ind1 = ind1+2;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry2[ind2] = (k-1)*zslice + i; # y=0
                bdry2[ind2+1] = k*zslice - xlen + i; # y = 1
                bdrynorm2[:,ind2] = [0,-1,0];
                bdrynorm2[:,ind2+1] = [0,1,0];
                ind2 = ind2+2;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry3[ind3] = (j-1)*xlen + i; # z=0
                bdry3[ind3+1] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm3[:,ind3] = [0,0,-1];
                bdrynorm3[:,ind3+1] = [0,0,1];
                ind3 = ind3+2;
            end
        end
        
        bdry = [bdry1, bdry2, bdry3];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1)*2);
        bdryfc2 = zeros(Int, (nx-1)*(nz-1)*2);
        bdryfc3 = zeros(Int, (nx-1)*(ny-1)*2);
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc1[ind1+1] = e2f[2,e2];
                ind1 = ind1+2;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc2[ind2] = e2f[3,e1];
                bdryfc2[ind2+1] = e2f[4,e2];
                ind2 = ind2+2;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc3[ind3] = e2f[5,e1];
                bdryfc3[ind3+1] = e2f[6,e2];
                ind3 = ind3+2;
            end
        end
        
        bdryfc = [bdryfc1, bdryfc2, bdryfc3];
        
    elseif bn == 2 # x=0, everything else
        bdry1 = zeros(Int, xbdrysize); # x=0
        bdry2 = zeros(Int, xbdrysize + ybdrysize*2 + zbdrysize*2); # other
        
        bdrynorm1 = zeros(3, xbdrysize); # x=0
        bdrynorm2 = zeros(3, xbdrysize + ybdrysize*2 + zbdrysize*2); # other
        
        ind1 = 1;
        ind2 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry2[ind2] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm2[:,ind2] = [1,0,0];
                ind1 = ind1+1;
                ind2 = ind2+1;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry2[ind2] = (k-1)*zslice + i; # y=0
                bdry2[ind2+1] = k*zslice - xlen + i; # y = 1
                bdrynorm2[:,ind2] = [0,-1,0];
                bdrynorm2[:,ind2+1] = [0,1,0];
                ind2 = ind2+2;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry2[ind2] = (j-1)*xlen + i; # z=0
                bdry2[ind2+1] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm2[:,ind2] = [0,0,-1];
                bdrynorm2[:,ind2+1] = [0,0,1];
                ind2 = ind2+2;
            end
        end
        
        bdry = [bdry1, bdry2];
        bdrynorm = [bdrynorm1, bdrynorm2];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1));
        bdryfc2 = zeros(Int, (ny-1)*(nz-1) + (nx-1)*(nz-1)*2 + (nx-1)*(ny-1)*2);
        ind1 = 1;
        ind2 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc2[ind2] = e2f[2,e2];
                ind1 = ind1+1;
                ind2 = ind2+1;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc2[ind2] = e2f[3,e1];
                bdryfc2[ind2+1] = e2f[4,e2];
                ind2 = ind2+2;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc2[ind2] = e2f[5,e1];
                bdryfc2[ind2+1] = e2f[6,e2];
                ind2 = ind2+2;
            end
        end
        
        bdryfc = [bdryfc1, bdryfc2];
        
    else #bn=1
        bdry1 = zeros(Int, xbdrysize*2 + ybdrysize*2 + zbdrysize*2); # all
        bdrynorm1 = zeros(3, xbdrysize*2 + ybdrysize*2 + zbdrysize*2); # all

        ind1 = 1;
        for k=1:zlen
            for j=1:ylen
                bdry1[ind1] = (k-1)*zslice + (j-1)*xlen + 1; # x=0
                bdry1[ind1+1] = (k-1)*zslice + j*xlen; # x = 1
                bdrynorm1[:,ind1] = [-1,0,0];
                bdrynorm1[:,ind1+1] = [1,0,0];
                ind1 = ind1+2;
            end
        end
        for k=1:zlen
            for i=2:xlen-1
                bdry1[ind1] = (k-1)*zslice + i; # y=0
                bdry1[ind1+1] = k*zslice - xlen + i; # y = 1
                bdrynorm1[:,ind1] = [0,-1,0];
                bdrynorm1[:,ind1+1] = [0,1,0];
                ind1 = ind1+2;
            end
        end
        for j=2:ylen-1
            for i=2:xlen-1
                bdry1[ind1] = (j-1)*xlen + i; # z=0
                bdry1[ind1+1] = N - zslice + (j-1)*xlen + i; # z = 1
                bdrynorm1[:,ind1] = [0,0,-1];
                bdrynorm1[:,ind1+1] = [0,0,1];
                ind1 = ind1+2;
            end
        end
        
        bdry = [bdry1];
        bdrynorm = [bdrynorm1];
        
        bdryfc1 = zeros(Int, (ny-1)*(nz-1)*2 + (nx-1)*(nz-1)*2 + (nx-1)*(ny-1)*2);
        ind1 = 1;
        for k=1:nz-1
            for j=1:ny-1
                e1 = (k-1)*(nx-1)*(ny-1) + (j-1)*(nx-1) + 1;
                e2 = (k-1)*(nx-1)*(ny-1) + j*(nx-1);
                bdryfc1[ind1] = e2f[1,e1];
                bdryfc1[ind1+1] = e2f[2,e2];
                ind1 = ind1+2;
            end
        end
        for k=1:nz-1
            for i=1:nx-1
                e1 = (k-1)*(nx-1)*(ny-1) + i;
                e2 = (k-1)*(nx-1)*(ny-1) + (ny-2)*(nx-1) + i;
                bdryfc1[ind1] = e2f[3,e1];
                bdryfc1[ind1+1] = e2f[4,e2];
                ind1 = ind1+2;
            end
        end
        for j=1:ny-1
            for i=1:nx-1
                e1 = (j-1)*(nx-1) + i;
                e2 = (nz-2)*(nx-1)*(ny-1) + (j-1)*(nx-1) + i;
                bdryfc1[ind1] = e2f[5,e1];
                bdryfc1[ind1+1] = e2f[6,e2];
                ind1 = ind1+2;
            end
        end
        
        bdryfc = [bdryfc1];
    end
    
    grid = Grid(x, bdry, bdryfc, bdrynorm, bids, loc2glb, glbvertex, f2glb, fvtx2glb);
    
    return (mesh, refel, refelfc, grid);
end
