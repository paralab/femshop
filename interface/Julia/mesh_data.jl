#=
# A struct containing mesh information
=#
export MeshData, build_faces, find_boundaries, find_normals

struct MeshData
    nx::Int;                    # Number of nodes
    nodes::Array{Float64,2};    # node locations
    indices::Array{Int,1};      # node indices may not be in order
    
    nel::Int;                   # Number of elements
    elements::Array{Int,2};     # Element vertex to node mapping
    etypes::Array{Int,1};       # Element types as defined by GMSH
    nv::Array{Int,1};           # Elements can have different numbers of vertices
    
    face2node::Array{Int,3}     # Nodes defining each face of each element
    normals::Array{Float64,3}   # Normal vectors for each face of each element
    bdry::Array{Int,2};         # Boundary ID for each face (0=interior face)
end

# numbers of nodes and faces for first order elements as defined by GMSH
# line, triangle, quad, tet, hex, prism, 5-pyramid
etypetonv = [2, 3, 4, 4, 8, 6, 5]; # number of vertices
etypetonf = [2, 3, 4, 4, 6, 5, 5]; # number of faces
etypetonfn= [1, 3, 4, 3, 4, 4, 4]; # number of nodes for each face (except prism and 5-pyramids!)


function build_faces(nel, elements, etypes)
    f2n = zeros(Int, nel, 6, 4); #(elementind, faceind, nodeind)
    for ei=1:nel
        if etypes[ei] == 1 # line
            f2n[ei,1,1] = elements[ei,1];
            f2n[ei,2,1] = elements[ei,2];
        elseif etypes[ei] == 2 # triangle
            f2n[ei,1,1:2] = elements[ei,1:2];
            f2n[ei,2,1:2] = elements[ei,2:3];
            f2n[ei,3,1:2] = elements[ei,[3 1]];
        elseif etypes[ei] == 3 # quad
            f2n[ei,1,1:2] = elements[ei,1:2];
            f2n[ei,2,1:2] = elements[ei,2:3];
            f2n[ei,3,1:2] = elements[ei,3:4];
            f2n[ei,4,1:2] = elements[ei,[4 1]];
        elseif etypes[ei] == 4 # tet
            f2n[ei,1,1:3] = elements[ei,[0 2 1]];
            f2n[ei,2,1:3] = elements[ei,[1 2 3]];
            f2n[ei,3,1:3] = elements[ei,[0 1 3]];
            f2n[ei,4,1:3] = elements[ei,[0 3 2]];
        elseif etypes[ei] == 5 # hex
            f2n[ei,1,1:4] = elements[ei,[0 4 7 3]];
            f2n[ei,2,1:4] = elements[ei,[1 2 6 5]];
            f2n[ei,3,1:4] = elements[ei,[0 1 5 4]];
            f2n[ei,4,1:4] = elements[ei,[2 3 7 6]];
            f2n[ei,5,1:4] = elements[ei,[0 3 2 1]];
            f2n[ei,6,1:4] = elements[ei,[4 5 6 7]];
        elseif etypes[ei] == 6 # prism
            f2n[ei,1,1:3] = elements[ei,[0 2 1]];
            f2n[ei,2,1:3] = elements[ei,[3 4 5]];
            f2n[ei,3,1:4] = elements[ei,[0 1 4 3]];
            f2n[ei,4,1:4] = elements[ei,[0 3 5 2]];
            f2n[ei,5,1:4] = elements[ei,[1 2 5 4]];
        elseif etypes[ei] == 7 # 5-pyramid
            f2n[ei,1,1:4] = elements[ei,[0 3 2 1]];
            f2n[ei,2,1:3] = elements[ei,[0 1 4]];
            f2n[ei,3,1:3] = elements[ei,[2 3 4]];
            f2n[ei,4,1:3] = elements[ei,[1 2 4]];
            f2n[ei,5,1:3] = elements[ei,[3 0 4]];
        end
    end
    return f2n;
end

function find_normals(nel, elements, etypes, faces, nodes)
    normal = zeros(nel, 6, 3); #(elementind, faceind, vector component)
    for ei=1:nel
        if etypes[ei] == 1 # line
            normal[ei,1,1] = -1;
            normal[ei,2,1] = 1;
        elseif etypes[ei] == 2 # triangle
            normal[ei,1,1:2] = normal2(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:]);
            normal[ei,2,1:2] = normal2(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:]);
            normal[ei,3,1:2] = normal2(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:]);
        elseif etypes[ei] == 3 # quad
            normal[ei,1,1:2] = normal2(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:]);
            normal[ei,2,1:2] = normal2(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:]);
            normal[ei,3,1:2] = normal2(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:]);
            normal[ei,4,1:2] = normal2(nodes[faces[ei,4,1],:], nodes[faces[ei,4,2],:]);
        elseif etypes[ei] == 4 # tet
            normal[ei,1,1:3] = normal3(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:], nodes[faces[ei,1,3],:]);
            normal[ei,2,1:3] = normal3(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:], nodes[faces[ei,2,3],:]);
            normal[ei,3,1:3] = normal3(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:], nodes[faces[ei,3,3],:]);
            normal[ei,4,1:3] = normal3(nodes[faces[ei,4,1],:], nodes[faces[ei,4,2],:], nodes[faces[ei,4,3],:]);
        elseif etypes[ei] == 5 # hex
            normal[ei,1,1:3] = normal3(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:], nodes[faces[ei,1,3],:]);
            normal[ei,2,1:3] = normal3(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:], nodes[faces[ei,2,3],:]);
            normal[ei,3,1:3] = normal3(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:], nodes[faces[ei,3,3],:]);
            normal[ei,4,1:3] = normal3(nodes[faces[ei,4,1],:], nodes[faces[ei,4,2],:], nodes[faces[ei,4,3],:]);
            normal[ei,5,1:3] = normal3(nodes[faces[ei,5,1],:], nodes[faces[ei,5,2],:], nodes[faces[ei,5,3],:]);
            normal[ei,6,1:3] = normal3(nodes[faces[ei,6,1],:], nodes[faces[ei,6,2],:], nodes[faces[ei,6,3],:]);
        elseif etypes[ei] == 6 # prism
            normal[ei,1,1:3] = normal3(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:], nodes[faces[ei,1,3],:]);
            normal[ei,2,1:3] = normal3(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:], nodes[faces[ei,2,3],:]);
            normal[ei,3,1:3] = normal3(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:], nodes[faces[ei,3,3],:]);
            normal[ei,4,1:3] = normal3(nodes[faces[ei,4,1],:], nodes[faces[ei,4,2],:], nodes[faces[ei,4,3],:]);
            normal[ei,5,1:3] = normal3(nodes[faces[ei,5,1],:], nodes[faces[ei,5,2],:], nodes[faces[ei,5,3],:]);
        elseif etypes[ei] == 7 # 5-pyramid
            normal[ei,1,1:3] = normal3(nodes[faces[ei,1,1],:], nodes[faces[ei,1,2],:], nodes[faces[ei,1,3],:]);
            normal[ei,2,1:3] = normal3(nodes[faces[ei,2,1],:], nodes[faces[ei,2,2],:], nodes[faces[ei,2,3],:]);
            normal[ei,3,1:3] = normal3(nodes[faces[ei,3,1],:], nodes[faces[ei,3,2],:], nodes[faces[ei,3,3],:]);
            normal[ei,4,1:3] = normal3(nodes[faces[ei,4,1],:], nodes[faces[ei,4,2],:], nodes[faces[ei,4,3],:]);
            normal[ei,5,1:3] = normal3(nodes[faces[ei,5,1],:], nodes[faces[ei,5,2],:], nodes[faces[ei,5,3],:]);
        end
    end
    return normal;
end

function normal2(a, b)
    nx = b[2]-a[2];
    ny = a[1]-b[1];
    d = sqrt(nx*nx + ny*ny);
    return [nx/d; ny/d];
end

function normal3(a, b, c)
    v = b.-a;
    w = c.-a;
    nx = v[2]*w[3]-v[3]*w[2];
    ny = v[3]*w[1]-v[1]*w[3];
    nz = v[1]*w[2]-v[3]*w[1];
    d = sqrt(nx*nx + ny*ny + nz*nz);
    return [nx/d; ny/d; nz/d];
end

function find_boundaries(nel, faces)
    bdry = zeros(Int, nel, 6);
    foundmatch = false;
    for ei=1:nel
        for fi=1:6
            if faces[ei,fi,1] > 0
                foundmatch = false;
                for ej=1:nel
                    if ei != ej
                        for fj=1:6
                            if faces[ej,fj,1] > 0
                                # This huge set of loops compares each face of each element with every other.
                                # If there are no two faes with the same set of nodes, it is a boundary.
                                if shared_face(faces[ei,fi,:], faces[ej,fj,:])
                                    foundmatch = true;
                                end
                            else
                                break;
                            end
                        end
                    end
                    if foundmatch
                        break;
                    end
                end
                if !foundmatch
                    bdry[ei, fi] = 612;
                end
            else
                break;
            end
        end
    end
    return bdry;
end

function shared_face(f1, f2)
    h = 0;
    s = 0;
    for i=1:4
        if f1[i] > 0
            h += 1;
            for j=1:4
                if f1[i] == f2[j]
                    s += 1;
                end
            end
        end
    end
    
    return s == h;
end