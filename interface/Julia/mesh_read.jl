#=
#
# Reads a .msh file and builds a MeshData struct
=#

export read_mesh, MeshData

struct MeshData
    nx::Int;                    # Number of nodes
    nodes::Array{Float64,2};    # node locations
    indices::Array{Int,1};      # node indices may not be in order
    
    nel::Int;                   # Number of elements
    elements::Array{Int,2};     # Element vertex to node mapping
    etypes::Array{Int,1};       # Element types as defined by GMSH
    nv::Array{Int,1};           # Elements can have different numbers of vertices
end

# numbers of nodes for first and second order elements as defined by GMSH
# TODO add higher order types
etypetonv = [2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13];

# Reads from the file stream
# Returns a MeshData struct
function read_mesh(file)
    # Determine the MSH format version
    msh_version = 2;
    while !eof(file)
        line = readline(file);
        if occursin("\$MeshFormat", line)
            # Check if the MSH version is old(2) or new(4)
            vals = split(readline(file), " ", keepempty=false);
            if parse(Float64, vals[1]) >= 4
                msh_version = 4;
            end
            break;
        end
    end
    
    # use the appropriate reader
    seekstart(file);
    if msh_version == 2
        return read_msh_v2(file);
    elseif msh_version == 4
        return read_msh_v4(file);
    end
end

function read_msh_v2(file)
    # Find the beginning of the Nodes and Elements sections.
    # Then read in the numbers.
    nodes_done = false;
    elements_done = false;
    msh_version = 2;
    nx = 0;
    nel = 0;
    nodes = [];
    indices = [];
    elements = [];
    etypes = [];
    nv = [];
    while((!nodes_done || !elements_done) && !eof(file))
        line = readline(file);
        if occursin("\$Nodes", line)
            # The Nodes section
            line = readline(file);
            nx = parse(Int, split(line, " ", keepempty=false)[1]);
            if nx > 0
                # parse each node's coordinates
                nodes = zeros(nx, 3);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file);
                while !occursin("\$EndNodes", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    indices[i] = parse(Int, vals[1]);
                    nodes[i,1] = parse(Float64, vals[2]);
                    nodes[i,2] = parse(Float64, vals[3]);
                    nodes[i,3] = parse(Float64, vals[4]);
                    i += 1;
                    line = readline(file);
                end
                
                nodes_done = true;
            end
        elseif occursin("\$Elements", line)
            # The Elements section
            line = readline(file);
            nel = parse(Int, split(line, " ", keepempty=false)[1]);
            if nel > 0
                # parse each element's numbers
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, nel, 8);
                i = 1;
                line = readline(file);
                while !occursin("\$EndElements", line) && !eof(file)
                    vals = split(line, " ", keepempty=false);
                    etypes[i] = parse(Int, vals[2]);
                    nv[i] = etypetonv[etypes[i]];
                    offset = parse(Int, vals[3]) + 3;
                    for j=1:nv[i]
                        elements[i,j] = parse(Int, vals[offset + j]);
                    end
                    i += 1;
                    line = readline(file);
                end
                
                elements_done = true;
            end
        end
    end
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end

function read_msh_v4(file)
    # Not ready yet
    return MeshData(0, [], [], 0, [], [], []);
end