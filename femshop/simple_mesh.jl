# A set of simple mesh makers for testing
export simple_line_mesh, simple_quad_mesh, simple_hex_mesh

#=
# Conveniently builds a uniform unit interval mesh.
# Simple just to get things working.
=#
function simple_line_mesh(nx, bn, interval)
    ord = config.basis_order_min;
    refel = build_refel(1, ord, 2, LOBATTO);
    Nv = nx;                    # number of vertex nodes
    N = (nx-1)*ord + 1;         # number of total nodes
    nel = nx-1;                 # number of elements
    Np = refel.Np;              # number of nodes per element
    el = zeros(Int, nel, 2);    # element vertex maps
    xv = zeros(Nv,1);           # coordinates of vertices
    x = zeros(N,1);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryel = [];                # index of elements touching each BID
    bdrynorm = [];              # normal at boundary nodes
    if bn == 2
        bids = [1,2]; # for two BIDs
    else
        bids = [1]; # for one BID
    end
    loc2glb = zeros(Int, nel,Np)# local to global index map for each element's nodes
    glbvertex = zeros(Int, nel,2);# local to global for vertices
    scale = interval[2]-interval[1];
    h = scale/(nx-1);               # uniformly divided

    # vertex nodes are ordered
    for j=1:nx
        xv[j,1] = interval[1] + (j-1)*h;
    end

    # Elements are ordered the same way
    for ei=1:(nx-1)
        x1 = interval[1] + (ei-1)*h; # left vertex
        glbvertex[ei,1] = (ei-1)*(Np-1) + 1;
        glbvertex[ei,2] = ei*(Np-1) + 1;
        
        for ni=1:Np-1
            gi = (ei-1)*(Np-1) + ni; # global index of this node
            x[gi,1] = x1 .+ h*0.5 .* (refel.r[ni] + 1); # coordinates of this node
            loc2glb[ei,ni] = gi; # local to global map
        end
        loc2glb[ei,Np] = ei*(Np-1) + 1;

        el[ei,1] = ei;
        el[ei,2] = ei+1;
    end
    # The last node was left out
    x[N] = interval[2];
    
    # indices are in order
    ind = Array(1:Nv);

    # boundary nodes
    if bn == 2
        bdry = [[1],[N]];
        bdryel = [[1],[nel]];
        bdry1 = Array{Float64,2}(undef,1,1);
        bdry1[1] = -1;
        bdry2 = Array{Float64,2}(undef,1,1);
        bdry2[1] = 1;
        bdrynorm = [bdry1,bdry2];
    else
        bdry = [[1,N]]; # for one BID
        bdryel = [[1,nel]];
        bdry1 = Array{Float64,2}(undef,2,1);
        bdry1[1] = -1;
        bdry1[2] = 1;
        bdrynorm = [bdry1];
    end

    mesh = MeshData(Nv, xv, ind, nel, el, ones(nel), 2*ones(nel)); # MeshData struct
    grid = Grid(x, bdry, bdryel, bdrynorm, bids, loc2glb, glbvertex);

    return (mesh, refel, grid);
end

#=
# Conveniently builds a square, uniform, unit quad mesh.
# Simple just to get things working.
=#
function simple_quad_mesh(nx, bn, interval)
    ord = config.basis_order_min;
    refel = build_refel(2, ord, 4, LOBATTO);
    Nv = nx*nx;                 # number of vertex nodes
    N = (nx-1)*ord + 1;
    N = N*N;                    # number of total nodes
    nel = (nx-1)*(nx-1);        # number of elements
    Np = refel.Np;              # number of nodes per element
    el = zeros(Int, nel, 4);    # element vertex maps
    xv = zeros(Nv,2);           # coordinates of vertices
    x = zeros(N,2);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryel = [];                # index of elements touching each BID
    bdrynorm = [];              # normal at boundary nodes
    if bn == 4
        bids = [1,2,3,4]; # x=0, x=1, y=0, y=1
    elseif bn == 3
        bids = [1,2,3]; # x=0 , x=1, y=0,1
    elseif bn == 2
        bids = [1,2]; # x=0,1 , y=0,1
    else
        bids = [1]; # everywhere
    end

    loc2glb = zeros(Int, nel,Np)# local to global index map for each element's nodes
    glbvertex = zeros(Int, nel,4);# local to global for vertices
    
    scale = interval[2]-interval[1];
    h = scale/(nx-1); # uniformly divided

    # vertex nodes are ordered along x then y
    for j=1:nx
        for i=1:nx
            k = i + (j-1)*nx;
            xv[k,1] = interval[1] + (i-1)*h;
            xv[k,2] = interval[1] + (j-1)*h;
        end
    end

    # Elements are ordered along x then y
    n1d = ord+1;
    rowsize = (nx-1)*(n1d-1) + 1;
    for j=1:(nx-1)
        for i=1:(nx-1)
            ei = i + (j-1)*(nx-1); # element index
            glbvertex[ei,1] = ((j-1)*(n1d-1))*rowsize + (i-1)*(n1d-1) + 1;
            glbvertex[ei,2] = ((j-1)*(n1d-1))*rowsize + i*(n1d-1) + 1;
            glbvertex[ei,3] = (j*(n1d-1))*rowsize + (i-1)*(n1d-1) + 1;
            glbvertex[ei,4] = (j*(n1d-1))*rowsize + i*(n1d-1) + 1;
            
            x1 = [interval[1] + (i-1)*h ; interval[1] + (j-1)*h]; # southwest corner
            for jj=1:n1d-1
                row = (j-1)*(n1d-1) + jj;
                for ii=1:n1d-1
                    gi = (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                    li = (jj-1)*n1d + ii; # local index of this node
                    x[gi,:] = x1 .+ h*0.5 .* (refel.r[li,:] + [1; 1]); # coordinates of this node
                    loc2glb[ei,li] = gi; # local to global map
                end
                loc2glb[ei,jj*n1d] = (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
            end
            row = j*(n1d-1)+1;
            for ii=1:n1d
                loc2glb[ei, n1d*(n1d-1)+ii] = ((row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
            end

            el[ei,1] = i + (j-1)*nx;
            el[ei,2] = i + (j-1)*nx + 1;
            el[ei,3] = i + (j)*nx;
            el[ei,4] = i + (j)*nx + 1;
        end
        # The last node on the row was not set
        x1 = [interval[1] + (nx-2)*h ; interval[1] + (j-1)*h]; # southwest corner of last element
        for jj=1:n1d-1
            li = jj*n1d; # local index of this node
            x[((j-1)*(n1d-1)+jj)*rowsize,:] = x1 .+ h*0.5 .* (refel.r[li,:] + [1; 1]); # coordinates of this node
        end
    end
    # The top row was not set
    rowstart = ((n1d-1)*(nx-1))*rowsize;
    for i=1:(nx-1)
        x1 = [interval[1] + (i-1)*h ; interval[1] + (nx-2)*h]; # southwest corner of top element
        for ii=1:n1d-1
            gi = rowstart + (i-1)*(n1d-1) + ii;
            li = (n1d-1)*n1d + ii; # local index of this node
            x[gi, :] = x1 .+ h*0.5 .* (refel.r[li,:] + [1; 1]); # coordinates of this node
        end
    end
    # corner
    x[N,:] = [interval[2]; interval[2]];

    # indices are in order
    ind = Array(1:Nv);

    # boundaries
    if bn == 4
        bdry1 = zeros(Int, rowsize); # x=0
        bdry2 = zeros(Int, rowsize); # x=1
        bdry3 = zeros(Int, rowsize - 2); # y=0
        bdry4 = zeros(Int, rowsize - 2); # y=1
        bdrynorm1 = Array{Float64,2}(undef,length(bdry1),2);
        bdrynorm2 = Array{Float64,2}(undef,length(bdry2),2);
        bdrynorm3 = Array{Float64,2}(undef,length(bdry3),2);
        bdrynorm4 = Array{Float64,2}(undef,length(bdry4),2);
        for i=1:rowsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[i,:] = [-1,0];
            bdry2[i] = i*rowsize; # right
            bdrynorm2[i,:] = [1,0];
            if i > 1 && i < rowsize
                bdry3[i-1] = i; # bottom
                bdrynorm3[i-1,:] = [0,-1];
                bdry4[i-1] = N-rowsize + i; # top
                bdrynorm4[i-1,:] = [0,1];
            end
        end
        bdry = [bdry1, bdry2, bdry3, bdry4];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3, bdrynorm4];
        
        bdryel1 = zeros(Int, nx-1);
        bdryel2 = zeros(Int, nx-1);
        bdryel3 = zeros(Int, nx-1);
        bdryel4 = zeros(Int, nx-1);
        for i=1:(nx-1)
            bdryel1[i] = (i-1)*(nx-1) + 1;
            bdryel2[i] = i*(nx-1);
            bdryel3[i] = i;
            bdryel4[i] = (nx-2)*(nx-1) + i;
        end
        bdryel = [bdryel1, bdryel2, bdryel3, bdryel4];

    elseif bn == 3
        bdry1 = zeros(Int, rowsize); # x=0
        bdry2 = zeros(Int, rowsize); # x=1
        bdry3 = zeros(Int, rowsize*2 - 4); # y=0,1
        bdrynorm1 = Array{Float64,2}(undef,length(bdry1),2);
        bdrynorm2 = Array{Float64,2}(undef,length(bdry2),2);
        bdrynorm3 = Array{Float64,2}(undef,length(bdry3),2);
        for i=1:rowsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[i,:] = [-1,0];
            bdry2[i] = i*rowsize; # right
            bdrynorm2[i,:] = [1,0];
            if i > 1 && i < rowsize
                bdry3[i-1] = i; # bottom
                bdrynorm3[i-1,:] = [0,-1];
                bdry3[i+rowsize-3] = N-rowsize + i; # top
                bdrynorm3[i+rowsize-3,:] = [0,1];
            end
        end
        bdry = [bdry1, bdry2, bdry3];
        bdrynorm = [bdrynorm1, bdrynorm2, bdrynorm3];
        
        bdryel1 = zeros(Int, nx-1);
        bdryel2 = zeros(Int, nx-1);
        bdryel3 = zeros(Int, (nx-1)*2);
        for i=1:(nx-1)
            bdryel1[i] = (i-1)*(nx-1) + 1;
            bdryel2[i] = i*(nx-1);
            bdryel3[i] = i;
            bdryel3[i+nx-1] = (nx-2)*(nx-1) + i;
        end
        bdryel = [bdryel1, bdryel2, bdryel3];

    elseif bn == 2
        bdry1 = zeros(Int, rowsize*2); # x=0,1
        bdry2 = zeros(Int, rowsize*2 - 4); # y=0,1
        bdrynorm1 = Array{Float64,2}(undef,length(bdry1),2);
        bdrynorm2 = Array{Float64,2}(undef,length(bdry2),2);
        for i=1:rowsize
            bdry1[i] = (i-1)*rowsize + 1; # left
            bdrynorm1[i,:] = [-1,0];
            bdry1[i + rowsize] = i*rowsize; # right
            bdrynorm1[i + rowsize,:] = [1,0];
            if i > 1 && i < rowsize
                bdry2[i-1] = i; # bottom
                bdrynorm2[i,:] = [0,-1];
                bdry2[i+rowsize-3] = N-rowsize + i; # top
                bdrynorm2[i+rowsize-3,:] = [0,1];
            end
        end
        bdry = [bdry1, bdry2];
        bdrynorm = [bdrynorm1, bdrynorm2];
        
        bdryel1 = zeros(Int, (nx-1)*2);
        bdryel2 = zeros(Int, (nx-1)*2);
        for i=1:(nx-1)
            bdryel1[i] = (i-1)*(nx-1) + 1;
            bdryel1[i+nx-1] = i*(nx-1);
            bdryel2[i] = i;
            bdryel2[i+nx-1] = (nx-2)*(nx-1) + i;
        end
        bdryel = [bdryel1, bdryel2];

    else
        bdry = zeros(Int, rowsize*4 - 4);
        bdrynorm1 = Array{Float64,2}(undef,length(bdry),2);
        for i=1:rowsize
            bdry[i] = i; # bottom
            bdrynorm1[i,:] = [0,-1];
            bdry[i+rowsize] = N-rowsize + i; # top
            bdrynorm1[i,:] = [0,1];
            if i > 1 && i < rowsize
                bdry[i-1 + rowsize*2] = (i-1)*rowsize + 1; # left
                bdrynorm1[i-1 + rowsize*2,:] = [-1,0];
                bdry[i-3 + rowsize*3] = i*rowsize; # right
                bdrynorm1[i-3 + rowsize*3,:] = [1,0];
            end
        end
        bdry = [bdry];
        bdrynorm = [bdrynorm1];
        
        bdryel = zeros(Int, (nx-1)*4 - 4);
        for i=1:(nx-2)
            bdryel[i] = i;
            bdryel[i+nx-2] = (nx-1)*(nx-2) + i + 1;
            bdryel[i+2*(nx-2)] = i*(nx-1) + 1;
            bdryel[i+3*(nx-2)] = i*(nx-1);
        end
        bdryel = [bdryel];
    end


    mesh = MeshData(Nv, xv, ind, nel, el, 3*ones(nel), 4*ones(nel)); # MeshData struct
    grid = Grid(x, bdry, bdryel, bdrynorm, bids, loc2glb, glbvertex);

    return (mesh, refel, grid);
end

#=
# Conveniently builds a uniform, unit cube, hex mesh.
# Simple just to get things working.
=#
function simple_hex_mesh(nx, bn, interval)
    ord = config.basis_order_min;
    refel = build_refel(3, ord, 6, LOBATTO);
    Nv = nx*nx*nx;                 # number of vertex nodes
    N = (nx-1)*ord + 1;
    N = N*N*N;                    # number of total nodes
    nel = (nx-1)*(nx-1)*(nx-1);        # number of elements
    Np = refel.Np;              # number of nodes per element
    el = zeros(Int, nel, 8);    # element vertex maps
    xv = zeros(Nv,3);           # coordinates of vertices
    x = zeros(N,3);             # coordinates of all nodes
    bdry = [];                  # index(in x) of boundary nodes for each BID
    bdryel = [];                # index of elements touching each BID
    bdrynorm = [];              # normal at boundary nodes
    bids = [];
    for i=1:bn
        push!(bids,i);
    end
    loc2glb = zeros(Int, nel, Np)# local to global index map for each element's nodes
    glbvertex = zeros(Int, nel,8);# local to global for vertices
    
    scale = interval[2]-interval[1];
    h = 1/(nx-1); # uniformly divided

    # Start with a 2d quad mesh
    (mesh2d, refel2d, grid2d) = simple_quad_mesh(nx, 1, interval);

    # vertex nodes are ordered along x then y then z
    for k=1:nx
        for j=1:nx
            for i=1:nx
                ind = i + (j-1)*nx + (k-1)*nx*nx;
                xv[ind,1] = interval[1] + (i-1)*h;
                xv[ind,2] = interval[1] + (j-1)*h;
                xv[ind,3] = interval[1] + (k-1)*h;
            end
        end
    end

    # Use 2d quad mesh to build hex mesh nodes
    n1d = ord+1;
    rowsize = (nx-1)*(n1d-1) + 1;
    slicesize = rowsize*rowsize;
    zvals = grid2d.allnodes[1:rowsize,1];
    for k=1:rowsize;
        range = ((k-1)*slicesize+1):(k*slicesize);
        x[range, :] = [grid2d.allnodes  zvals[k].*ones(slicesize)];
    end

    # Elements are ordered along x then y then z
    for k=1:(nx-1)
        for j=1:(nx-1) # element indices
            for i=1:(nx-1)
                ei = i + (j-1)*(nx-1) + (k-1)*(nx-1)*(nx-1); # element index
                glbvertex[ei,1] = (k-1)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[ei,2] = (k-1)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[ei,3] = (k-1)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[ei,4] = (k-1)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[ei,5] = (k)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[ei,6] = (k)*(n1d-1)*slicesize + (j-1)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                glbvertex[ei,7] = (k)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i-1)*(n1d-1) + 1;
                glbvertex[ei,8] = (k)*(n1d-1)*slicesize + (j)*(n1d-1)*rowsize + (i)*(n1d-1) + 1;
                
                for kk=1:n1d-1
                    slice = (k-1)*(n1d-1) + kk;
                    for jj=1:n1d-1 # elemental node indices
                        row = (j-1)*(n1d-1) + jj;
                        for ii=1:n1d-1
                            gi = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                            li = (kk-1)*n1d*n1d + (jj-1)*n1d + ii; # local index of this node
                            loc2glb[ei,li] = gi; # local to global map
                        end
                        loc2glb[ei,(kk-1)*n1d*n1d + jj*n1d] = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
                    end

                    row = j*(n1d-1)+1;
                    for ii=1:n1d
                        loc2glb[ei, (kk-1)*n1d*n1d + n1d*(n1d-1)+ii] = ((slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
                    end
                end

                slice = k*(n1d-1)+1;  # last in local z
                for jj=1:n1d-1
                    row = (j-1)*(n1d-1) + jj;
                    for ii=1:n1d-1
                        gi = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + ii; # global index of this node
                        li = (n1d-1)*n1d*n1d + (jj-1)*n1d + ii; # local index of this node
                        loc2glb[ei,li] = gi; # local to global map
                    end
                    loc2glb[ei,(n1d-1)*n1d*n1d + jj*n1d] = (slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1) + n1d; # last in local x
                end
                row = j*(n1d-1)+1;
                for ii=1:n1d
                    loc2glb[ei, (n1d-1)*n1d*n1d + n1d*(n1d-1)+ii] = ((slice-1)*slicesize + (row-1)*rowsize + (i-1)*(n1d-1)) + ii; # last in local y
                end

                el[ei,1] = i + (j-1)*nx + (k-1)*nx*nx;
                el[ei,2] = i + (j-1)*nx + (k-1)*nx*nx + 1;
                el[ei,3] = i + (j)*nx + (k-1)*nx*nx;
                el[ei,4] = i + (j)*nx + (k-1)*nx*nx + 1;
                el[ei,5] = i + (j-1)*nx + (k)*nx*nx;
                el[ei,6] = i + (j-1)*nx + (k)*nx*nx + 1;
                el[ei,7] = i + (j)*nx + (k)*nx*nx;
                el[ei,8] = i + (j)*nx + (k)*nx*nx + 1;
            end
        end
    end

    # indices are in order
    indexorder = Array(1:Nv);

    # boundaries
    if bn == 6 # all separate
        bdry1 = zeros(Int, slicesize); # x=0
        bdry2 = zeros(Int, slicesize); # x=1
        bdry3 = zeros(Int, slicesize - 2*rowsize); # y=0
        bdry4 = zeros(Int, slicesize - 2*rowsize); # y=1
        bdry5 = zeros(Int, slicesize - 4*rowsize); # z=0
        bdry6 = zeros(Int, slicesize - 4*rowsize); # z=1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for j=1:rowsize
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry1[ind1] = offset + (i-1)*rowsize + 1; # x=0
                bdry2[ind1] = offset + i*rowsize; # x = 1
                ind1 = ind1+1;
                if i > 1 && i < rowsize
                    bdry3[ind2] = offset + i;# y=0
                    bdry4[ind2] = offset + (rowsize-1)*rowsize + i;# y=1
                    ind2 = ind2+1;
                end
                if j > 1 && j < rowsize && i > 1 && i < rowsize
                    bdry5[ind3] = (j-1)*rowsize + i;# z=0
                    bdry6[ind3] = (rowsize-1)*slicesize + (j-1)*rowsize + i;# z=1
                    ind3 = ind3+1;
                end

            end
        end
        bdry = [bdry1, bdry2, bdry3, bdry4, bdry5, bdry6];
        
        bdryel1 = zeros(Int, (nx-1)*(nx-1));
        bdryel2 = zeros(Int, (nx-1)*(nx-1));
        bdryel3 = zeros(Int, (nx-1)*(nx-1));
        bdryel4 = zeros(Int, (nx-1)*(nx-1));
        bdryel5 = zeros(Int, (nx-1)*(nx-1));
        bdryel6 = zeros(Int, (nx-1)*(nx-1));
        ind1 = 1;
        ind3 = 1;
        ind5 = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel1[ind1] = xoffset; # x=0
                bdryel2[ind1] = xoffset + nx-2; # x=1
                ind1 = ind1 + 1;
                
                yoffset = (j-1)*NN + i;
                bdryel3[ind3] = yoffset; # y=0
                bdryel4[ind3] = yoffset + (nx-1)*(nx-2); # y=1
                ind3 = ind3 + 1;
                
                zoffset = (j-1)*(nx-1) + i;
                bdryel5[ind5] = zoffset; # z=0
                bdryel6[ind5] = NNN - NN + zoffset; # z=1
                ind5 = ind5 + 1;
            end
        end
        
        bdryel = [bdryel1, bdryel2, bdryel3, bdryel4, bdryel5, bdryel6];
        
    elseif bn == 5 # combine zs
        bdry1 = zeros(Int, slicesize); # x=0
        bdry2 = zeros(Int, slicesize); # x=1
        bdry3 = zeros(Int, slicesize - 2*rowsize); # y=0
        bdry4 = zeros(Int, slicesize - 2*rowsize); # y=1
        bdry5 = zeros(Int, 2*slicesize - 8*rowsize); # z=0,1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for j=1:rowsize
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry1[ind1] = offset + (i-1)*rowsize + 1; # x=0
                bdry2[ind1] = offset + i*rowsize; # x = 1
                ind1 = ind1+1;
                if i > 1 && i < rowsize
                    bdry3[ind2] = offset + i;# y=0
                    bdry4[ind2] = offset + (rowsize-1)*rowsize + i;# y=1
                    ind2 = ind2+1;
                end
                if j > 1 && j < rowsize && i > 1 && i < rowsize
                    bdry5[ind3] = (j-1)*rowsize + i;# z=0
                    bdry5[ind3+1] = (rowsize-1)*slicesize + (j-1)*rowsize + i;# z=1
                    ind3 = ind3+2;
                end

            end
        end
        bdry = [bdry1, bdry2, bdry3, bdry4, bdry5];
        
        bdryel1 = zeros(Int, (nx-1)*(nx-1));
        bdryel2 = zeros(Int, (nx-1)*(nx-1));
        bdryel3 = zeros(Int, (nx-1)*(nx-1));
        bdryel4 = zeros(Int, (nx-1)*(nx-1));
        bdryel5 = zeros(Int, (nx-1)*(nx-1)*2);
        ind1 = 1;
        ind3 = 1;
        ind5 = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel1[ind1] = xoffset; # x=0
                bdryel2[ind1] = xoffset + nx-2; # x=1
                ind1 = ind1 + 1;
                
                yoffset = (j-1)*NN + i;
                bdryel3[ind3] = yoffset; # y=0
                bdryel4[ind3] = yoffset + (nx-1)*(nx-2); # y=1
                ind3 = ind3 + 1;
                
                zoffset = (j-1)*(nx-1) + i;
                bdryel5[ind5] = zoffset; # z=0
                bdryel5[ind5+1] = NNN - NN + zoffset; # z=1
                ind5 = ind5 + 2;
            end
        end
        
        bdryel = [bdryel1, bdryel2, bdryel3, bdryel4, bdryel5];
        
    elseif bn == 4 # combine ys and zs
        bdry1 = zeros(Int, slicesize); # x=0
        bdry2 = zeros(Int, slicesize); # x=1
        bdry3 = zeros(Int, 2*slicesize - 4*rowsize); # y=0,1
        bdry4 = zeros(Int, 2*slicesize - 8*rowsize); # z=0,1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for j=1:rowsize
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry1[ind1] = offset + (i-1)*rowsize + 1; # x=0
                bdry2[ind1] = offset + i*rowsize; # x = 1
                ind1 = ind1+1;
                if i > 1 && i < rowsize
                    bdry3[ind2] = offset + i;# y=0
                    bdry3[ind2+1] = offset + (rowsize-1)*rowsize + i;# y=1
                    ind2 = ind2+2;
                end
                if j > 1 && j < rowsize && i > 1 && i < rowsize
                    bdry4[ind3] = (j-1)*rowsize + i;# z=0
                    bdry4[ind3+1] = (rowsize-1)*slicesize + (j-1)*rowsize + i;# z=1
                    ind3 = ind3+2;
                end

            end
        end
        bdry = [bdry1, bdry2, bdry3, bdry4];
        
        bdryel1 = zeros(Int, (nx-1)*(nx-1));
        bdryel2 = zeros(Int, (nx-1)*(nx-1));
        bdryel3 = zeros(Int, (nx-1)*(nx-1)*2);
        bdryel4 = zeros(Int, (nx-1)*(nx-1)*2);
        ind1 = 1;
        ind3 = 1;
        ind4 = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel1[ind1] = xoffset; # x=0
                bdryel2[ind1] = xoffset + nx-2; # x=1
                ind1 = ind1 + 1;
                
                yoffset = (j-1)*NN + i;
                bdryel3[ind3] = yoffset; # y=0
                bdryel3[ind3+1] = yoffset + (nx-1)*(nx-2); # y=1
                ind3 = ind3 + 2;
                
                zoffset = (j-1)*(nx-1) + i;
                bdryel4[ind4] = zoffset; # z=0
                bdryel4[ind4+1] = NNN - NN + zoffset; # z=1
                ind4 = ind4 + 2;
            end
        end
        
        bdryel = [bdryel1, bdryel2, bdryel3, bdryel4];
        
    elseif bn == 3 # xs, ys, zs
        bdry1 = zeros(Int, 2*slicesize); # x=0,1
        bdry2 = zeros(Int, 2*slicesize - 4*rowsize); # y=0,1
        bdry3 = zeros(Int, 2*slicesize - 8*rowsize); # z=0,1

        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        for j=1:rowsize
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry1[ind1] = offset + (i-1)*rowsize + 1; # x=0
                bdry1[ind1+1] = offset + i*rowsize; # x = 1
                ind1 = ind1+2;
                if i > 1 && i < rowsize
                    bdry2[ind2] = offset + i;# y=0
                    bdry2[ind2+1] = offset + (rowsize-1)*rowsize + i;# y=1
                    ind2 = ind2+2;
                end
                if j > 1 && j < rowsize && i > 1 && i < rowsize
                    bdry3[ind3] = (j-1)*rowsize + i;# z=0
                    bdry3[ind3+1] = (rowsize-1)*slicesize + (j-1)*rowsize + i;# z=1
                    ind3 = ind3+2;
                end

            end
        end
        bdry = [bdry1, bdry2, bdry3];
        
        bdryel1 = zeros(Int, (nx-1)*(nx-1)*2);
        bdryel2 = zeros(Int, (nx-1)*(nx-1)*2);
        bdryel3 = zeros(Int, (nx-1)*(nx-1)*2);
        ind1 = 1;
        ind2 = 1;
        ind3 = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel1[ind1] = xoffset; # x=0
                bdryel1[ind1+1] = xoffset + nx-2; # x=1
                ind1 = ind1 + 2;
                
                yoffset = (j-1)*NN + i;
                bdryel2[ind2] = yoffset; # y=0
                bdryel2[ind2+1] = yoffset + (nx-1)*(nx-2); # y=1
                ind2 = ind2 + 2;
                
                zoffset = (j-1)*(nx-1) + i;
                bdryel3[ind3] = zoffset; # z=0
                bdryel3[ind3+1] = NNN - NN + zoffset; # z=1
                ind3 = ind3 + 2;
            end
        end
        
        bdryel = [bdryel1, bdryel2, bdryel3];
        
    elseif bn == 2 # x=0, everything else
        bdry1 = zeros(Int, slicesize); # x=0
        bdry2 = zeros(Int, 5*slicesize - 12*rowsize + 8); # the rest
        bdrynorm1 = Array{Float64,2}(undef,length(bdry1),3);
        bdrynorm2 = Array{Float64,2}(undef,length(bdry2),3);

        ind1 = 1;
        ind2 = 1;
        for j=1:rowsize
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry1[ind1] = offset + (i-1)*rowsize + 1; # x=0
                bdrynorm1[ind1,:] = [-1,0,0];
                ind1 = ind1+1;
                bdry2[ind2] = offset + i*rowsize; # x = 1
                bdrynorm2[ind2,:] = [1,0,0];
                ind2 = ind2+1;
                if i > 1 && i < rowsize
                    bdry2[ind2] = offset + i;# y=0
                    bdry2[ind2+1] = offset + (rowsize-1)*rowsize + i;# y=1
                    bdrynorm2[ind2,:] = [0,-1,0];
                    bdrynorm2[ind2+1,:] = [0,1,0];
                    ind2 = ind2+2;
                    
                    if j > 1 && j < rowsize 
                        bdry2[ind2] = (j-1)*rowsize + i;# z=0
                        bdry2[ind2+1] = (rowsize-1)*slicesize + (j-1)*rowsize + i;# z=1
                        bdrynorm2[ind2,:] = [0,0,-1];
                        bdrynorm2[ind2+1,:] = [0,0,1];
                        ind2 = ind2+2;
                    end
                end

            end
        end
        bdry = [bdry1, bdry2];
        bdrynorm = [bdrynorm1, bdrynorm2];
        
        bdryel1 = zeros(Int, (nx-1)*(nx-1));
        bdryel2 = zeros(Int, (nx-1)*(nx-1) + (nx-1)*(nx-2)*2 + (nx-3)*(nx-2)*2);
        ind1 = 1;
        ind2 = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel1[ind1] = xoffset; # x=0
                bdryel2[ind2] = xoffset + nx-2; # x=1
                ind2 = ind2 + 1;
                ind1 = ind1 + 1;
                
                if i < (nx-1)
                    yoffset = (j-1)*NN + i;
                    bdryel2[ind2] = yoffset; # y=0
                    bdryel2[ind2+1] = yoffset + (nx-1)*(nx-2); # y=1
                    ind2 = ind2 + 2;
                    
                    if j > 1 && j < (nx-1)
                        zoffset = (j-1)*(nx-1) + i;
                        bdryel2[ind2] = zoffset; # z=0
                        bdryel2[ind2+1] = NNN - NN + zoffset; # z=1
                        ind2 = ind2 + 2;
                    end
                end
            end
        end
        
        bdryel = [bdryel1, bdryel2];
    else
        N1d = (nx-1)*ord + 1;
        bdry = zeros(Int, N - (N1d-2)*(N1d-2)*(N1d-2));
        bdrynorm1 = Array{Float64,2}(undef,length(bdry),3);
        ind = 1;
        for i=1:slicesize
            bdry[i] = i; # z=0
            bdry[i+slicesize] = N-i+1; # z=1
            bdrynorm1[i,:] = [0,0,-1];
            bdrynorm1[i+slicesize,:] = [0,0,1];
        end
        ind = 2*slicesize + 1;
        for j=2:rowsize-1
            offset = (j-1)*slicesize;
            for i=1:rowsize
                bdry[ind] = offset + i; # y=0
                bdry[ind+1] = offset + slicesize-rowsize + i; # y=1
                bdrynorm1[ind,:] = [0,-1,0];
                bdrynorm1[ind+1,:] = [0,1,0];
                ind += 2;
                if i > 1 && i < rowsize
                    bdry[ind] = offset + (i-1)*rowsize + 1; # x=0
                    bdry[ind + 1] = offset + i*rowsize; # x=1
                    bdrynorm1[ind,:] = [-1,0,0];
                    bdrynorm1[ind+1,:] = [1,0,0];
                    ind += 2;
                end
            end
        end
        bdry = [bdry];
        bdrynorm = [bdrynorm1];
        
        bdryel = zeros(Int, (nx-1)*(nx-1)*2 + (nx-1)*(nx-3)*2 + (nx-3)*(nx-3)*2);
        ind = 1;
        NN = (nx-1)^2;
        NNN = (nx-1)^3;
        for j=1:nx-1
            for i=1:nx-1
                xoffset = (j-1)*NN + (i-1)*(nx-1) + 1;
                bdryel[ind] = xoffset; # x=0
                bdryel[ind+1] = xoffset + nx-2; # x=1
                ind = ind + 2;
                
                if i > 1 && i < (nx-1)
                    yoffset = (j-1)*NN + i;
                    bdryel[ind] = yoffset; # y=0
                    bdryel[ind+1] = yoffset + (nx-1)*(nx-2); # y=1
                    ind = ind + 2;
                    
                    if j > 1 && j < (nx-1)
                        zoffset = (j-1)*(nx-1) + i;
                        bdryel[ind] = zoffset; # z=0
                        bdryel[ind+1] = NNN - NN + zoffset; # z=1
                        ind = ind + 2;
                    end
                end
            end
        end
        
        bdryel = [bdryel];
    end

    mesh = MeshData(Nv, xv, indexorder, nel, el, 5*ones(nel), 8*ones(nel)); # MeshData struct
    grid = Grid(x, bdry, bdryel, bdrynorm, bids, loc2glb, glbvertex);

    return (mesh, refel, grid);
end
