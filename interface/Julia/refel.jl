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
    
    r1d::Array{Float64}     # Node coordinates in 1D
    r::Array{Float64}       # dim-dim Node coordinates
    
    wr1d::Array{Float64}    # r1d gll Quadrature weights
    wr::Array{Float64}      # r gll Quadrature weights
    
    g1d::Array{Float64}     # 1D Gauss points
    wg1d::Array{Float64}    # 1D Gauss weights
    
    g::Array{Float64}       # dim-dim Gauss points
    wg::Array{Float64}      # dim-dim Gauss weights
    
    V::Array{Float64}       # Vandermonde matrix at r
    gradV::Array{Float64}   # grad of basis at r
    invV::Array{Float64}    # Inverse V
    
    Vg::Array{Float64}      # Vandermonde matrix at Gauss
    gradVg::Array{Float64}  # grad of basis at g
    invVg::Array{Float64}   # Inverse Vg
    
    Dr::Array{Float64}      # Differentiation matrix for r
    Dg::Array{Float64}      # Differentiation matrix for g
    
    # Useful matrices, use them like so
    # for integral(basis_j*basis_i) -> invV' * integral(modal_j*modal_i) * invV -> invV'*Vg'*diag(w*detJ)*Vg*invV -> Q'*diag(w*detJ)*Q
    # for integral(gradbasis_j dot gradbasis_i) -> invV'*gradVg'*Jacobian'*diag(w*detJ)*Jacobian*gradVg*invV 
    #                                           -> [QrQsQt]*Jacobian' * diag(w*detJ) * Jacobian*[QrQsQt]'
    Q1d::Array{Float64}     # 1D quadrature matrix: like Vg*invV
    Q::Array{Float64}       # dim-dim quadrature matrix
    Qr::Array{Float64}      # quad of derivative matrix: like gradVg*invV
    Qs::Array{Float64}      # 
    Qt::Array{Float64}      # 
    
    # for DG
    lift::Array{Float64}    # Surface integral matrix
    
    # Constructor needs at least this information
    Refel(dim, order, nnodes, nfaces) = new(
        dim,
        order,
        nnodes,
        nfaces,
        [],[],
        [],[],
        [],[],
        [],[],
        [],[],[],
        [],[],[],
        [],[],
        [],[],[],[],[],
        []
    )
end

function build_refel(dimension, order, nfaces, nodetype)
    # Check for errors
    if (dimension == 1 && nfaces != 2) || nfaces < dimension+1 || order < 1
        println("Error: buildRefel(dimension, order, nfaces, nodetype), check for valid parameters.");
        return nothing;
    end
    # Number of points determined by order element type
    if (dimension == 1)     Np = order+1; # line segment
    elseif (dimension == 2 && nfaces == 3) Np = (Int)((order+1)*(order+2)/2); # triangle
    elseif (dimension == 2 && nfaces == 4) Np = (Int)((order+1)*(order+1)); # quad
    elseif (dimension == 3 && nfaces == 4) Np = (Int)((order+1)*(order+2)*(order+3)/6); # tet
    elseif (dimension == 3 && nfaces == 6) Np = (Int)((order+1)*(order+1)*(order+1)); # hex
    elseif (dimension == 4) Np = (Int)((order+1)*(order+2)*(order+3)*(order+4)/24); # ??
    end
    
    refel = Refel(dimension, order, Np, nfaces);
    
    # Get nodes on the reference element
    refel_nodes!(refel, nodetype);
    
    # Vandermonde matrix and grad,inv
    refel.V = zeros(order+1, order+1);
    refel.gradV = zeros(order+1, order+1);
    for i=1:refel.N+1
        refel.V[:,i] = jacobi_polynomial(refel.r1d, 0, 0, i-1);
    end
    for i=1:refel.N
        refel.gradV[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(refel.r1d, 1, 1, i-1);
    end
    refel.invV = inv(refel.V);
    
    # Gauss versions
    refel.Vg = zeros(order+1, order+1);
    refel.gradVg = zeros(order+1, order+1);
    for i=1:refel.N+1
        refel.Vg[:,i] = jacobi_polynomial(refel.g1d, 0, 0, i-1);
    end
    for i=1:refel.N
        refel.gradVg[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(refel.g1d, 1, 1, i-1);
    end
    refel.invVg = inv(refel.Vg);
    
    # Differentiation matrices
    refel.Dr = refel.gradV*refel.invV;
    refel.Dg = refel.gradVg*refel.invV;
    
    refel.Q1d = refel.Vg*refel.invV;
    
    if dimension == 1
        refel.Q = refel.Q1d;
        refel.Qr = refel.Dg;
    elseif dimension == 2
        refel.Q = kron(refel.Q1d,refel.Q1d);
        refel.Qr = kron(refel.Q1d,refel.Dg);
        refel.Qs = kron(refel.Dg,refel.Q1d);
    elseif dimension == 3
        refel.Q = kron(kron(refel.Q1d, refel.Q1d), refel.Q1d);
        refel.Qr = kron(kron(refel.Q1d, refel.Q1d), refel.Dg);
        refel.Qs = kron(kron(refel.Q1d, refel.Dg), refel.Q1d);
        refel.Qt = kron(kron(refel.Dg, refel.Q1d), refel.Q1d);
    end
    
    # # Build surface integral matrix
    # if dimension == 1
    #     edgemat = zeros(refel.Np, 2);
    #     edgemat[1,1] = 1;
    #     edgemat[refel.Np,2] = 1;
    #     refel.lift = refel.V * (refel.V' * edgemat);
    # else
    #     # not ready
    # end
    
    return refel;
end

