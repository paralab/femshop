#=
#
# Reads a .msh file and builds a MeshData struct
=#
export read_mesh

# numbers of CORNER nodes for first and second order elements as defined by GMSH
# TODO add higher order types
etypetonv = [2, 3, 4, 4, 8, 6, 5, 2, 3, 4, 4, 8, 6, 5, 1, 4, 8, 6, 5];

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
            #=
            $Nodes
            numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
            entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
                nodeTag(size_t)
                ...
                x(double) y(double) z(double)
                < u(double; if parametric and entityDim >= 1) >
                < v(double; if parametric and entityDim >= 2) >
                < w(double; if parametric and entityDim == 3) >
                ...
            ...
            $EndNodes
            =#
            line = readline(file);
            nx = parse(Int, split(line, " ", keepempty=false)[2]);
            if nx > 0
                # parse node info
                nodes = zeros(nx, 3);
                indices = zeros(Int, nx);
                i = 1;
                line = readline(file);
                while !occursin("\$EndNodes", line) && !eof(file)
                    # process entity blocks 
                    # The first line has entity info. 
                    entnx = parse(Int, split(line, " ", keepempty=false)[4]);
                    # Node indices
                    for ni=1:entnx
                        line = readline(file);
                        vals = split(line, " ", keepempty=false);
                        indices[i-1 + ni] = parse(Int, vals[1]);
                    end
                    # Node coordinates
                    for ni=1:entnx
                        line = readline(file);
                        vals = split(line, " ", keepempty=false);
                        indices[i-1 + ni] = parse(Int, vals[1]);
                        nodes[i-1 + ni,1] = parse(Float64, vals[1]);
                        nodes[i-1 + ni,2] = parse(Float64, vals[2]);
                        nodes[i-1 + ni,3] = parse(Float64, vals[3]);
                    end
                    i += entnx;
                    line = readline(file);
                end
                
                nodes_done = true;
            end
        elseif occursin("\$Elements", line)
            # The Elements section
            #=
            $Elements
            1 2 1 2          1 entity bloc, 2 elements total, min/max element tags: 1 and 2
            2 1 3 2          2D entity (surface) 1, element type 3 (4-node quad), 2 elements
            1 1 2 3 4          quad tag #1, nodes 1 2 3 4
            2 2 5 6 3          quad tag #2, nodes 2 5 6 3
            $EndElements
            =#
            line = readline(file);
            nel = parse(Int, split(line, " ", keepempty=false)[2]);
            if nel > 0
                # parse elements
                nv = zeros(Int, nel);
                etypes = zeros(Int, nel);
                elements = zeros(Int, nel, 8);
                i = 1;
                line = readline(file);
                while !occursin("\$EndElements", line) && !eof(file)
                    # Process entity blocks
                    # The first line has entity info. 
                    vals = split(line, " ", keepempty=false);
                    entnel = parse(Int, vals[4]);
                    enttype = parse(Int, vals[3]);
                    entnv = etypetonv[etypes[i]];
                    # read element info
                    for ei=1:entnel
                        line = readline(file);
                        etypes[i-1 + ei] = enttype;
                        nv[i-1 + ei] = entnv;
                        for j=1:entnv
                            elements[i-1 + ei,j] = parse(Int, vals[1 + j]);
                        end
                    end
                    
                    i += entnel;
                    line = readline(file);
                end
                
                elements_done = true;
            end
        end
    end
    
    return MeshData(nx, nodes, indices, nel, elements, etypes, nv);
end
