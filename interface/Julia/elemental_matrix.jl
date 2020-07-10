#=
# Elemental matrices
# This code is based on HOMG
=#
include("tensor_ops.jl");

function elemental_mass(refel, detJ)
    return refel.Q' * diagm(refel.wg .* detJ) * refel.Q;
end

function elemental_stiffness(refel, detJ, J, gpts, coef=0)
    #             | Qx Qy Qz || rx ry rz |     | rx sx tx || Qx |
    #    Ke =                 | sx sy sz | J W | ry sy ty || Qy |
    #                         | tx ty tz |     | rz sz tz || Qz |
    
    #gpts = gauss points for this element, used for eval of coefficients
    
    if coef == 0
        mu = 1;
    else
        #TODO coef
        mu = 1;
    end
    
    if refel.dim == 1
        factor = (J.rx.*J.rx) .* detJ .* refel.wg .* mu ; # d2u/dx^2
        
        Ke =  refel.Qx' * diagm(factor) * refel.Qx;
        
    elseif refel.dim == 2
        nn = length(detJ);
        factor = zeros(nn, 3);
        #             1  3
        # 2D factor   3  2
        factor[:,1] = (J.rx.*J.rx + J.ry.*J.ry ) .* detJ .* refel.wg .* mu ; # d2u/dx^2
        factor[:,2] = (J.sx.*J.sx + J.sy.*J.sy ) .* detJ .* refel.wg .* mu ; # d2u/dy^2
        factor[:,3] = (J.rx.*J.sx + J.ry.*J.sy ) .* detJ .* refel.wg .* mu ; # d2u/dxdy
        
        Ke =  refel.Qx' * diagm(factor[:,1]) * refel.Qx +
              refel.Qy' * diagm(factor[:,2]) * refel.Qy +
              refel.Qx' * diagm(factor[:,3]) * refel.Qy +
              refel.Qy' * diagm(factor[:,3]) * refel.Qx ;
            
    else
        nn = length(detJ);
        factor = zeros(nn, 6);
        #             1  4  5
        # 3D factor   4  2  6
        #             5  6  3
        # first compute dj.w.J.J'
        factor[:,1] = (J.rx.*J.rx + J.ry.*J.ry + J.rz.*J.rz ) .* detJ .* refel.wg .* mu ; # d2u/dx^2
        factor[:,2] = (J.sx.*J.sx + J.sy.*J.sy + J.sz.*J.sz ) .* detJ .* refel.wg .* mu ; # d2u/dy^2
        factor[:,3] = (J.tx.*J.tx + J.ty.*J.ty + J.tz.*J.tz ) .* detJ .* refel.wg .* mu ; # d2u/dz^2
        
        factor[:,4] = (J.rx.*J.sx + J.ry.*J.sy + J.rz.*J.sz ) .* detJ .* refel.wg .* mu ; # d2u/dxdy
        factor[:,5] = (J.rx.*J.tx + J.ry.*J.ty + J.rz.*J.tz ) .* detJ .* refel.wg .* mu ; # d2u/dxdz
        factor[:,6] = (J.sx.*J.tx + J.sy.*J.ty + J.sz.*J.tz ) .* detJ .* refel.wg .* mu ; # d2u/dydz
        
        Ke =  refel.Qx' * diagm(factor[:,1]) * refel.Qx +
              refel.Qy' * diagm(factor[:,2]) * refel.Qy +
              refel.Qz' * diagm(factor[:,3]) * refel.Qz +
              refel.Qx' * diagm(factor[:,4]) * refel.Qy +
              refel.Qy' * diagm(factor[:,4]) * refel.Qx +
              refel.Qx' * diagm(factor[:,5]) * refel.Qz +
              refel.Qz' * diagm(factor[:,5]) * refel.Qx +
              refel.Qz' * diagm(factor[:,6]) * refel.Qy +
              refel.Qy' * diagm(factor[:,6]) * refel.Qz ;
    end
    
    return Ke;
end

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
