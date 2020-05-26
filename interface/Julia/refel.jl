#=
# A reference element.
# Notation follows that used in Nodal Discontinuous Galerkin Methods
# by Hesthaven and Warburton.
#  https://link.springer.com/book/10.1007/978-0-387-72067-8
#
# Can be used for 1, 2, 3, 4 dimensions
=#
include("refel_nodes.jl");
include("jacobi_polynomial.jl");

mutable struct Refel
    dim::Int                # Dimension
    N::Int                  # Order of polynomials
    Np::Int                 # Number of nodes
    Nfaces::Int             # Number of faces
    
    r::Array{Float64}       # Node coordinates in r (dim 1)
    s::Array{Float64}       # Node coordinates in s (dim 2)
    t::Array{Float64}       # Node coordinates in t (dim 3)
    u::Array{Float64}       # Node coordinates in u (dim 4)
    
    V::Array{Float64}       # Vandermonde matrix
    invV::Array{Float64}    # Inverse V
    mass::Array{Float64}    # Mass matrix
    Dr::Array{Float64}      # Differentiation matrix wrt r
    Ds::Array{Float64}      # Differentiation matrix wrt s (empty for 1D)
    Dt::Array{Float64}      # Differentiation matrix wrt t (empty for 1,2D)
    Du::Array{Float64}      # Differentiation matrix wrt u (empty for 1~3D)
    
    # for DG
    lift::Array{Float64}    # Surface integral matrix
    
    # Constructor needs at least this information
    Refel(dim, order, nnodes, nfaces) = new(
        dim,
        order,
        nnodes,
        nfaces,
        [],[],[],[],
        [],[],[],[],[],[],[],[]
    )
end

function build_refel(dimension, order, nfaces, nodetype)
    # Check for errors
    if (dimension == 1 && nfaces != 2) || nfaces < dimension+1 || order < 1
        println("Error: buildRefel(dimension, order, nfaces, nodetype), check for valid parameters.");
        return nothing;
    end
    # Number of points determined by order and dimension
    if (dimension == 1)     Np = order+1;
    elseif (dimension == 2) Np = (Int)((order+1)*(order+2)/2);
    elseif (dimension == 3) Np = (Int)((order+1)*(order+2)*(order+3)/6);
    elseif (dimension == 4) Np = (Int)((order+1)*(order+2)*(order+3)*(order+4)/24);
    end
    
    refel = Refel(dimension, order, Np, nfaces);
    
    # Get nodes on the reference element
    refel_nodes!(refel, nodetype);
    
    # Build Vandermonde matrix and inverse
    if refel.dim == 1
        refel.V = zeros(refel.Np, refel.order+1);
        for i=1:refel.order+1
            if config.trial_function == LEGENDRE
                refel.V[:,i] = jacobi_polynomial(refel.r, 0, 0, i-1);
            end
        end
        refel.invV = inv(V);
    else
        # not ready
    end
    
    # Build Mass matrix
    refel.mass = inv(V*V');
    
    # Build differentiation matrices
    if dimension == 1
        refel.Dr = zeros(refel.Np, refel.order+1);
        refel.Dr[:,1] = zeros(refel.Np);
        for i=2:refel.order+1
            if config.trial_function == LEGENDRE
                n = refel.order-1;
                refel.Dr[:,i] = sqrt(n*(n+1)) .* jacobi_polynomial(refel.r, 1, 1, i-2);
            end
        end
    else
        # not ready
    end
    
    # Build surface integral matrix
    if dimension == 1
        edgemat = zeros(refel.Np, 2);
        edgemat[1,1] = 1;
        edgemat[refel.Np,2] = 1;
        refel.lift = refel.V * (refel.V' * edgemat);
    else
        # not ready
    end
    
    return refel;
end

