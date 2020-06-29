#=
# Compute the node locations for the reference element
=#
include("jacobi_gauss_quad.jl");

function refel_nodes!(refel, nodetype)
    if refel.dim == 1
        # 1D has segments
        if nodetype == UNIFORM
            refel.r = Array(-1:(2/(refel.Np-1)):1);
        elseif nodetype == GAUSS
            (r,w) = jacobi_gauss_quad(0,0,refel.N);
            refel.r = r;
            refel.wr = w;
        elseif nodetype == LOBATTO
            if refel.N == 1
                refel.r = [-1; 1];
                refel.wr = [1; 1];
            else
                (r,w) = jacobi_gauss_quad(1,1,refel.N-2);
                refel.r = [-1; r ; 1];
                
                # compute the weights
                w = jacobi_polynomial(r, 0, 0, refel.N)
                adgammaN = (2*refel.N + 1) / (refel.N * (refel.N + 1))
                w = w.*w
                w = adgammaN./w
                
                refel.wr = w;
            end
        end
    elseif refel.dim == 2
        # 2D has triangles and quads
        
    elseif refel.dim == 3
        # 3D has tets, hexs and prisms
        
    else
        # Not ready
    end
end

# Returns the global node locations for elemental nodes
# Has size Np*dim
function get_node_coords(vx, refel)
    x = [];
    if refel.dim == 1
        hx = vx[2] - vx[1];
        x = vx[1] .+ (refel.r .+ 1) .* (hx*0.5);
    elseif refel.dim == 2
        
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
        
    elseif refel.dim == 3
        
    end
    
    return vol;
end