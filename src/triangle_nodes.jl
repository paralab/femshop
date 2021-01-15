#=
Find nodes for triangular elements.
Use the book.
=#
include("triangle_quadrature_table.jl");

# Build global nodes for a triangular element with refel and vertices v
function triangle_element_nodes(refel, v)
    return  triangle_refel_to_xy(refel.r[:,1], refel.r[:,2], v);
end

# Set up refel nodal array
function triangle_refel_nodes!(refel)
    (eqx, eqy) = triangle_equilateral_nodes(refel.N);
    (r, s) = triangle_equilateral_to_rs(eqx,eqy);
    
    refel.r = zeros(refel.Np,2);
    refel.wr = zeros(refel.Np);
    #refel.g = zeros(refel.Np,2);
    #refel.wg = zeros(refel.Np);
    
    refel.r[:,1] = r;
    refel.r[:,2] = s;
    
    # quadrature nodes/weights from a table
    xyw = triangle_quadrature_nodes_weights(refel.N);
    refel.g = xyw[:,1:2];
    refel.wg = xyw[:,3];
    
end

# Purpose  : Compute (x,y) nodes in equilateral triangle for polynomial of order N  
function triangle_equilateral_nodes(N)
    # optimal alpha values for N up to 16 (from book)
    alpopt = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];
            
    # Set optimized parameter, alpha, depending on order N
    if (N<16)
        alpha = alpopt[N];
    else
        alpha = 5/3;
    end

    # total number of nodes
    Np = (N+1)*(N+2)/2;
    Np = Int64(Np)
    # Create equidistributed nodes on equilateral triangle
    L1 = zeros(Np); 
    L2 = zeros(Np); 
    L3 = zeros(Np);
    sk = 1;
    for n=1:(N+1)
        for m=1:(N+2-n)
            L1[sk] = (n-1)/N; 
            L3[sk] = (m-1)/N;
            sk = sk+1;
        end
    end
    L2 = 1.0 .- L1 .- L3;
    x = -L2+L3; 
    y = (-L2-L3+2*L1) * 0.5773502691896258; # 1/sqrt(3)=0.5773502691896258

    # Compute blending function at each node for each edge
    blend1 = 4*L2.*L3; 
    blend2 = 4*L1.*L3; 
    blend3 = 4*L1.*L2;

    # Amount of warp for each node, for each edge
    warpf1 = triangle_warpfactor(N,L3-L2); 
    warpf2 = triangle_warpfactor(N,L1-L3); 
    warpf3 = triangle_warpfactor(N,L2-L1);

    # Combine blend & warp
    warp1 = blend1.*warpf1.*(1 .+ (alpha*L1).^2);
    warp2 = blend2.*warpf2.*(1 .+ (alpha*L2).^2);
    warp3 = blend3.*warpf3.*(1 .+ (alpha*L3).^2);

    # Accumulate deformations associated with each edge
    # x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
    # y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
    x = x + warp1 - 0.5*warp2 - 0.5*warp3;
    y = y + 0.8660254037844387*warp2 - 0.8660254037844387*warp3;
    
    return (x, y);
end

# Compute scaled warp function at order N based on rout interpolation nodes
function triangle_warpfactor(N, rout)
    # Compute LGL and equidistant node distribution
    (LGLr, w) = jacobi_gauss_quad(0,0,N);
    h = 2/N;
    req = -1:h:1;
    
    # Compute Vandermonde based on req
    Veq = zeros(N+1, N+1);
    for i=1:N+1
        Veq[:,i] = jacobi_polynomial(req, 0, 0, i-1);
    end
    
    # Evaluate Lagrange polynomial at rout
    Nr = length(rout); 
    Pmat = zeros(N+1,Nr);
    for i=1:N+1
      Pmat[i,:] = jacobi_polynomial(rout, 0, 0, i-1);
    end
    Lmat = Veq'\Pmat;
    
    # Compute warp factor
    warp = Lmat'*(LGLr - req);
    
    # Scale factor
    sf = zeros(length(rout))
    for i=1:length(rout)
        sf[i] = (abs(rout[i]) < 0.9999999999) ? (1-rout[i]*rout[i]) : 1
    end
    ## warp = warp./sf + warp.*(zerof-1);   zerof?
    
    return warp;
end

# From (x,y) in equilateral triangle to (r,s) coordinates in standard triangle
function triangle_equilateral_to_rs(x,y)
    # 1/sqrt(3)=0.5773502691896258
    L1 = (0.5773502691896258*y .+ 1.0)/3.0;
    L2 = (-3.0*x - 0.5773502691896258*y .+ 2.0)/6.0;
    L3 = ( 3.0*x - 0.5773502691896258*y .+ 2.0)/6.0;
    
    r = -L2 + L3 - L1; 
    s = -L2 - L3 + L1;
    
    return (r, s);
end

# From (r,s) coordinates in reference triangle to (x,y) in triangle with vertices v
# v is a 2x3 array [x1 x2 x3; y1 y2 y3]
function triangle_refel_to_xy(r, s, v)
    x = 0.5*(-(r+s)*v[1,1]+(1+r)*v[1,2]+(1+s)*v[1,3]);
    y = 0.5*(-(r+s)*v[2,1]+(1+r)*v[2,2]+(1+s)*v[2,3]);
    
    return (x, y);
end