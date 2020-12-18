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
    
    Ddr::Array{Float64}      # Similar to Qr, but for the elemental nodes, not quadrature nodes
    Dds::Array{Float64}      # 
    Ddt::Array{Float64}      # 
    
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
        [],[],[]
    )
end

import Base.copy
function copy(ref::Refel)
    newref = Refel(ref.dim, ref.N, ref.Np, ref.Nfaces);
    newref.r1d = copy(ref.r1d);
    newref.r = copy(ref.r);
    newref.wr1d = copy(ref.wr1d);
    newref.wr = copy(ref.wr);
    newref.g1d = copy(ref.g1d);
    newref.wg1d = copy(ref.wg1d);
    newref.g = copy(ref.g);
    newref.wg = copy(ref.wg);
    newref.V = copy(ref.V);
    newref.gradV = copy(ref.gradV);
    newref.invV = copy(ref.invV);
    newref.Vg = copy(ref.Vg);
    newref.gradVg = copy(ref.gradVg);
    newref.invVg = copy(ref.invVg);
    newref.Dr = copy(ref.Dr);
    newref.Dg = copy(ref.Dg);
    newref.Q1d = copy(ref.Q1d);
    newref.Q = copy(ref.Q);
    newref.Qr = copy(ref.Qr);
    newref.Qs = copy(ref.Qs);
    newref.Qt = copy(ref.Qt);
    newref.Ddr = copy(ref.Ddr);
    newref.Dds = copy(ref.Dds);
    newref.Ddt = copy(ref.Ddt);
    
    return newref;
end

function build_refel(dimension, order, nfaces, nodetype)
    # Check for errors
    if (dimension == 1 && nfaces != 2) || (dimension > 1 && nfaces < dimension+1)
        println("Error: buildRefel(dimension, order, nfaces, nodetype), check for valid parameters.");
        return nothing;
    end
    
    # Number of points determined by order element type
    if (dimension == 0) Np = 1; #point
    elseif (dimension == 1)     Np = order+1; # line segment
    elseif (dimension == 2 && nfaces == 3) Np = (Int)((order+1)*(order+2)/2); # triangle
    elseif (dimension == 2 && nfaces == 4) Np = (Int)((order+1)*(order+1)); # quad
    elseif (dimension == 3 && nfaces == 4) Np = (Int)((order+1)*(order+2)*(order+3)/6); # tet
    elseif (dimension == 3 && nfaces == 6) Np = (Int)((order+1)*(order+1)*(order+1)); # hex
    elseif (dimension == 4) Np = (Int)((order+1)*(order+2)*(order+3)*(order+4)/24); # ??
    end
    
    refel = Refel(dimension, order, Np, nfaces);
    
    # Get nodes on the reference element
    refel_nodes!(refel, nodetype);
    
    #  0D refels
    if dimension == 0
        refel.V = ones(1,1);
        refel.gradV = zeros(1,1);
        refel.invV = ones(1,1);
        refel.Vg = ones(1,1);
        refel.gradVg = zeros(1,1);
        refel.invVg = ones(1,1);
        refel.Dr = zeros(1,1);
        refel.Dg = zeros(1,1);
        refel.Q1d = zeros(1,1);
        refel.Q = ones(1,1);
        refel.Qr = zeros(1,1);
        refel.Qs = zeros(1,1);
        refel.Qt = zeros(1,1);
        refel.Ddr = zeros(1,1);
        refel.Dds = zeros(1,1);
        refel.Ddt = zeros(1,1);
        
        return refel;
    end
    
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
        refel.Ddr = refel.Dr;
    elseif dimension == 2
        ident = Matrix(1.0*I,order+1,order+1);
        refel.Q = kron(refel.Q1d,refel.Q1d);
        refel.Qr = kron(refel.Q1d,refel.Dg);
        refel.Qs = kron(refel.Dg,refel.Q1d);
        refel.Ddr = kron(ident,refel.Dr);
        refel.Dds = kron(refel.Dr,ident);
    elseif dimension == 3
        ident = Matrix(1.0*I,order+1,order+1);
        refel.Q = kron(kron(refel.Q1d, refel.Q1d), refel.Q1d);
        refel.Qr = kron(kron(refel.Q1d, refel.Q1d), refel.Dg);
        refel.Qs = kron(kron(refel.Q1d, refel.Dg), refel.Q1d);
        refel.Qt = kron(kron(refel.Dg, refel.Q1d), refel.Q1d);
        refel.Ddr = kron(kron(ident, ident), refel.Dg);
        refel.Dds = kron(kron(ident, refel.Dg), ident);
        refel.Ddt = kron(kron(refel.Dg, ident), ident);
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

function custom_quadrature_refel(oldrefel, nodes, weights)
    # The supplied nodes will replace the gauss quadrature nodes and weights r.g and r.wg
    # NOTE: Elemental node and weight arrays will not be changed, so quadrature must be set to GAUSS.
    # NOTE: Only quadrature matrices are changed. Nothing else.
    refel = copy(oldrefel);
    refel.g = copy(nodes);
    refel.wg = copy(weights);
    
    # Vandermonde matrix and grad,inv for r will not change
    
    # Gauss versions will change
    Nn = size(nodes,2);
    refel.Vg = zeros(Nn, refel.N+1);
    refel.gradVg = zeros(Nn, refel.N+1);
    if refel.dim == 1
        for i=1:refel.N+1
            refel.Vg[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
        end
        for i=1:refel.N
            refel.gradVg[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
        end
        
        # Differentiation matrices
        refel.Dg = refel.gradVg*refel.invV;
        refel.Q1d = refel.Vg*refel.invV;
        
        # Quadrature matrices
        refel.Q = refel.Q1d;
        refel.Qr = refel.Dg;
        refel.Ddr = refel.Dr;
        
    elseif refel.dim == 2
        # Now it's a little trickier, so forget the efficiency and do it carefully.
        Nn = size(nodes,2);
        tmp1 = zeros(Nn, refel.N+1);
        tmp2 = zeros(Nn, refel.N+1);
        Dtmp1 = zeros(Nn, refel.N+1);
        Dtmp2 = zeros(Nn, refel.N+1);
        
        for i=1:refel.N+1
            tmp1[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
            tmp2[:,i] = jacobi_polynomial(nodes[2,:], 0, 0, i-1);
        end
        for i=1:refel.N
            Dtmp1[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
            Dtmp2[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[2,:], 1, 1, i-1);
        end
        
        # Use kron(Vg*invV, Vg,invV) = kron(Vg,Vg)*kron(invV,invV)
        VgXVg = zeros(Nn, (refel.N+1)*(refel.N+1));
        VgXDg = zeros(Nn, (refel.N+1)*(refel.N+1));
        DgXVg = zeros(Nn, (refel.N+1)*(refel.N+1));
        for i=1:Nn
            for j=1:(refel.N+1)
                for k=1:(refel.N+1)
                    ind = (k-1)*(refel.N+1) + j;
                    VgXVg[i,ind] = tmp1[i,j]*tmp2[i,k];
                    VgXDg[i,ind] = Dtmp1[i,j]*tmp2[i,k];
                    DgXVg[i,ind] = tmp1[i,j]*Dtmp2[i,k];
                end
            end
        end
        viXvi = kron(refel.invV, refel.invV);
        
        refel.Q = VgXVg * viXvi;
        refel.Qr = VgXDg * viXvi;
        refel.Qs = DgXVg * viXvi;
        
    elseif refel.dim == 3
        Nn = size(nodes,2);
        tmp1 = zeros(Nn, refel.N+1);
        tmp2 = zeros(Nn, refel.N+1);
        tmp3 = zeros(Nn, refel.N+1);
        Dtmp1 = zeros(Nn, refel.N+1);
        Dtmp2 = zeros(Nn, refel.N+1);
        Dtmp3 = zeros(Nn, refel.N+1);
        
        for i=1:refel.N+1
            tmp1[:,i] = jacobi_polynomial(nodes[1,:], 0, 0, i-1);
            tmp2[:,i] = jacobi_polynomial(nodes[2,:], 0, 0, i-1);
            tmp3[:,i] = jacobi_polynomial(nodes[3,:], 0, 0, i-1);
        end
        for i=1:refel.N
            Dtmp1[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[1,:], 1, 1, i-1);
            Dtmp2[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[2,:], 1, 1, i-1);
            Dtmp3[:,i+1] = sqrt(i*(i+1)) .* jacobi_polynomial(nodes[3,:], 1, 1, i-1);
        end
        
        viXviXvi = kron(kron(refel.invV, refel.invV), refel.invV);
        
        VgXVgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        VgXVgXDg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        VgXDgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        DgXVgXVg = zeros(Nn, (refel.N+1)*(refel.N+1)*(refel.N+1));
        for ni=1:Nn
            for i=1:(refel.N+1)
                for j=1:(refel.N+1)
                    for k=1:(refel.N+1)
                        ind = ((k-1)*(refel.N+1) + (j-1))*(refel.N+1) + i;
                        VgXVgXVg[ni,ind] = tmp1[ni,i]*tmp2[ni,j]*tmp2[ni,k];
                        VgXVgXDg[ni,ind] = Dtmp1[ni,i]*tmp2[ni,j]*tmp2[ni,k];
                        VgXDgXVg[ni,ind] = tmp1[ni,i]*Dtmp2[ni,j]*tmp2[ni,k];
                        DgXVgXVg[ni,ind] = tmp1[ni,i]*tmp2[ni,j]*Dtmp2[ni,k];
                    end
                end
            end
        end
        
        refel.Q = VgXVgXVg * viXviXvi;
        refel.Qr = VgXVgXDg * viXviXvi;
        refel.Qs = VgXDgXVg * viXviXvi;
        refel.Qt = DgXVgXVg * viXviXvi;
    end
    
    return refel;
    
end

# Map 1D nodes to 2D face
# Assumes interval is scaled the same 
function map_face_nodes_2d(g1d, v1, v2)
    fnodes = zeros(2,length(g1d));
    dx = v2[1] - v1[1];
    dy = v2[2] - v1[2];
    dist = sqrt(dx*dx+dy*dy);
    xmult = dx/dist;
    ymult = dy/dist;
    
    fnodes[1,:] = v1[1] .+ (g1d .- v1[1]) .*xmult;
    fnodes[2,:] = v1[2] .+ (g1d .- v1[2]) .*ymult;
    
    return fnodes;
end

# Map 2D nodes to 3D face
# Assumes g2d[1] aligns with v1 and interval is scaled the same 
function map_face_nodes_3d(g2d, v1, v2, v3, v4=[])
    fnodes = zeros(3,length(g2d));
    if length(v4) > 0 # quads
        dx = 0
        dy = 0
        dz = 0
        dist = sqrt(dx*dx+dy*dy+dz*dz);
        xmult = dx/dist;
        ymult = dy/dist;
        zmult = dz/dist;
        
        fnodes[1,:] = v1[1] .+ (g2d .- v1[1]) .*xmult;
        fnodes[2,:] = v1[2] .+ (g2d .- v1[2]) .*ymult;
        fnodes[3,:] = v1[3] .+ (g2d .- v1[3]) .*zmult;
    else #triangles
        
    end
    
    
    return fnodes;
end