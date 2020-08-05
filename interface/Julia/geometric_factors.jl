# Geometric factors
export geometric_factors, Jacobian

include("tensor_ops.jl");

#=
# Stores the elemental jacobian
# Used by the "geometric_factors" function
=#
struct Jacobian
    rx; ry; rz; sx; sy; sz; tx; ty; tz;
end

function geometric_factors(refel, pts)
    # pts = element node global coords
    # J = detJ
    # D = Jacobian
    if refel.dim == 1
        xr  = refel.Dg*pts;
        J = xr[:];
        rx = 1 ./ J;
        D = Jacobian(rx,[],[],[],[],[],[],[],[]);
        
    elseif refel.dim == 2
        (xr, xs) = tensor_grad2(refel.Dg, pts[:,1]);
        (yr, ys) = tensor_grad2(refel.Dg, pts[:,2]);
        J = -xs.*yr + xr.*ys;
        
        rx =  ys./J;
        sx = -yr./J;
        ry = -xs./J;
        sy =  xr./J;
        D = Jacobian(rx,ry,[],sx,sy,[],[],[],[]);
        
    else
        (xr, xs, xt) = tensor_grad3(refel.Dg, pts[:,1]);
        (yr, ys, yt) = tensor_grad3(refel.Dg, pts[:,2]);
        (zr, zs, zt) = tensor_grad3(refel.Dg, pts[:,3]);
        J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
        
        rx =  (ys.*zt - zs.*yt)./J;
        ry = -(xs.*zt - zs.*xt)./J;
        rz =  (xs.*yt - ys.*xt)./J;
        
        sx = -(yr.*zt - zr.*yt)./J;
        sy =  (xr.*zt - zr.*xt)./J;
        sz = -(xr.*yt - yr.*xt)./J;
        
        tx =  (yr.*zs - zr.*ys)./J;
        ty = -(xr.*zs - zr.*xs)./J;
        tz =  (xr.*ys - yr.*xs)./J;
        D = Jacobian(rx,ry,rz,sx,sy,sz,tx,ty,tz);
    end
    
    return (J,D);
end