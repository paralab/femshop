#=
# Builds a local to global mapping.
# Each element has a mapping from its nodes to a global index.
# mesh is a MeshData struct, refel is a Refel struct
# For now assume all elements share the same refel.
=#
function build_global_map(mesh, refel)
    nel = mesh.nel;
    Np = refel.Np;
    glbl = zeros(nel, Np);
    offset = 0;
    etypetonf = [2, 3, 4, 4, 6, 5, 5, 2, 3, 4, 4, 6, 5, 5, 1, 4, 6, 5, 5]; # number of faces for gmsh element types
    nfaces = etypetonf[mesh.etypes[1]]; # assume all similar
    
    # The first element gets the first Np spots
    glbl[1,1:Np] = Array(1:Np);
    
    # Successive elements check for neighbors with lower element index.
    # If nodes are shared with a lower neighbor, use those values.
    for ei=2:nel
        # check neighbors
        for ni=1:nfaces
            if mesh.neighbors[ei,ni] < ei
                # wait a minute.... what am I doing? This is not the point.
            end
        end
    end
    
    return glbl
end