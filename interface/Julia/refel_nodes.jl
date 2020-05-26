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
            refel.r = jacobi_gauss_quad(0,0,refel.order)[1];
        elseif nodetype == LOBATTO
            if order == 1
                refel.r = [-1; 1];
            else
                refel.r = [-1; jacobi_gauss_quad(1,1,refel.order-2)[1] ; 1];
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
