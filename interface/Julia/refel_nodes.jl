#=
# Compute the node locations for the reference element
=#
include("jacobi_gauss_quad.jl");

function refel_nodes!(refel, nodetype)
    if refel.dim == 1
        # 1D has segments
        if nodetype == UNIFORM
            refel.r1d = Array(-1:(2/(refel.Np-1)):1);
            refel.wr1d = ones(length(refel.r1d)) ./ length(refel.r1d);
        elseif nodetype == GAUSS
            (r,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.r1d = r;
            refel.wr1d = w;
        elseif nodetype == LOBATTO
            if refel.N == 1
                refel.r1d = [-1; 1];
                refel.wr1d = [1; 1];
            else
                (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                refel.r1d = [-1; r ; 1];
                
                # compute the weights
                w = jacobi_polynomial(refel.r1d, 0, 0, refel.N)
                adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1))
                w = w.*w
                w = adgammaN./w
                
                refel.wr1d = w;
            end
        end
        # Then find Gauss points
        (g,w) = jacobi_gauss_quad(0,0,refel.N);
        refel.g1d = g;
        refel.wg1d = w;
        
        # 1D is easy
        refel.g = refel.g1d;
        refel.wg = refel.wg1d;
        refel.r = refel.r1d;
        refel.wr = refel.wr1d;
        
    elseif refel.dim == 2
        # 2D has triangles and quads
        
        # for now assume rectangular quads
        if nodetype == UNIFORM
            refel.r1d = Array(-1:(2/(refel.Np-1)):1);
            refel.wr1d = ones(length(refel.r1d)) ./ length(refel.r1d);
        elseif nodetype == GAUSS
            (r,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.r1d = r;
            refel.wr1d = w;
        elseif nodetype == LOBATTO
            if refel.N == 1
                refel.r1d = [-1; 1];
                refel.wr1d = [1; 1];
            else
                (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                refel.r1d = [-1; r ; 1];
                
                # compute the weights
                w = jacobi_polynomial(refel.r1d, 0, 0, refel.N);
                adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1));
                w = w.*w;
                w = adgammaN./w;
                
                refel.wr1d = w;
            end
        end
        # Then find Gauss points
        (g,w) = jacobi_gauss_quad(0,0,refel.N);
        refel.g1d = g;
        refel.wg1d = w;
        
        # r and w
        refel.r = zeros(refel.Np,2);
        refel.wr = zeros(refel.Np);
        refel.g = zeros(refel.Np,2);
        refel.wg = zeros(refel.Np);
        n1d = length(refel.r1d);
        for j=1:n1d
            for i=1:n1d
                k = (j-1)*n1d + i;
                refel.r[k,:] = [refel.r1d[i]; refel.r1d[j]];
                refel.wr[k] = refel.wr1d[i] * refel.wr1d[j];
                refel.g[k,:] = [refel.g1d[i]; refel.g1d[j]];
                refel.wg[k] = refel.wg1d[i] * refel.wg1d[j];
            end
        end
        
    elseif refel.dim == 3
        # 3D has tets, hexs and prisms
        
    else
        # Not ready
    end
end

# Returns the global node locations for elemental nodes
# Has size Np,dim
function get_node_coords(vx, elnodes)
    x = [];
    if refel.dim == 1
        hx = vx[2] - vx[1];
        x = vx[1] .+ (elnodes .+ 1) .* (hx*0.5);
    elseif refel.dim == 2
        # for now assume rectangular quads
        x = zeros(size(elnodes));
        hx = max(vx[2,1]-vx[1,1], vx[3,1]-vx[1,1]);
        hy = max(vx[2,2]-vx[1,2], vx[3,2]-vx[1,2]);
        x[:,1] = vx[1,1] .+ (elnodes[:,1] .+ 1) .*(hx*0.5);
        x[:,2] = vx[1,2] .+ (elnodes[:,2] .+ 1) .*(hy*0.5);
    elseif refel.dim == 3
        
    end
    
    return x;
end

# Returns the volume of an element
function get_volume(refel, vx)
    vol = 0;
    if refel.dim == 1
        vol = abs(vx[2]-vx[1]);
    elseif refel.dim == 2
        # for now assume rectangular quads
        hx = max(vx[2,1]-vx[1,1], vx[3,1]-vx[1,1]);
        hy = max(vx[2,2]-vx[1,2], vx[3,2]-vx[1,2]);
        vol = abs(hx*hy);
    elseif refel.dim == 3
        
    end
    
    return vol;
end