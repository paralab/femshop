# Geometric factors
export geometric_factors, Jacobian, build_deriv_matrix

include("tensor_ops.jl");

#=
# Stores the elemental jacobian
# Used by the "geometric_factors" function
=#
struct Jacobian
    rx; ry; rz; sx; sy; sz; tx; ty; tz;
end
import Base.copy
function copy(j::Jacobian)
    return Jacobian(copy(j.rx),copy(j.ry),copy(j.rz),copy(j.sx),copy(j.sy),copy(j.sz),copy(j.tx),copy(j.ty),copy(j.tz));
end

function geometric_factors(refel, pts)
    # pts = element node global coords
    # J = detJ
    # D = Jacobian
    if refel.dim == 0
        J = [1];
        D = Jacobian([1],[],[],[],[],[],[],[],[]);
        
    elseif refel.dim == 1
        xr  = refel.Dg*pts[:];
        J = xr[:];
        rx = 1 ./ J;
        D = Jacobian(rx,[],[],[],[],[],[],[],[]);
        
    elseif refel.dim == 2
        (xr, xs) = tensor_grad2(refel.Dg, pts[1,:][:]);
        (yr, ys) = tensor_grad2(refel.Dg, pts[2,:][:]);
        J = -xs.*yr + xr.*ys;
        
        rx =  ys./J;
        sx = -yr./J;
        ry = -xs./J;
        sy =  xr./J;
        D = Jacobian(rx,ry,[],sx,sy,[],[],[],[]);
        
    else
        (xr, xs, xt) = tensor_grad3(refel.Dg, pts[1,:][:]);
        (yr, ys, yt) = tensor_grad3(refel.Dg, pts[2,:][:]);
        (zr, zs, zt) = tensor_grad3(refel.Dg, pts[3,:][:]);
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

function build_deriv_matrix(refel, J)
    if refel.dim == 1
        RQ1 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j];
            end
        end
        return (RQ1,RD1);
        
    elseif refel.dim == 2
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j];
            end
        end
        return (RQ1, RQ2, RD1, RD2);
        
    elseif refel.dim == 3
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RQ3 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        RD3 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j] + J.tx[i]*refel.Qt[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j] + J.ty[i]*refel.Qt[i,j];
                RQ3[i,j] = J.rz[i]*refel.Qr[i,j] + J.sz[i]*refel.Qs[i,j] + J.tz[i]*refel.Qt[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j] + J.tx[i]*refel.Ddt[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j] + J.ty[i]*refel.Ddt[i,j];
                RD3[i,j] = J.rz[i]*refel.Ddr[i,j] + J.sz[i]*refel.Dds[i,j] + J.tz[i]*refel.Ddt[i,j];
            end
        end
        return (RQ1, RQ2, RQ3, RD1, RD2, RD3);
    end
end

function build_face_deriv_matrix(refel, J)
    if refel.dim == 1
        RQ1 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j];
            end
        end
        return (RQ1,);
        
    elseif refel.dim == 2
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j];
            end
        end
        return (RQ1, RQ2);
        
    elseif refel.dim == 3
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RQ3 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j] + J.tx[i]*refel.Qt[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j] + J.ty[i]*refel.Qt[i,j];
                RQ3[i,j] = J.rz[i]*refel.Qr[i,j] + J.sz[i]*refel.Qs[i,j] + J.tz[i]*refel.Qt[i,j];
            end
        end
        return (RQ1, RQ2, RQ3);
    end
end

###################################################################################################################################
function geometric_factors_cachesim(refel, pts)
    # pts = element node global coords
    # J = detJ
    # D = Jacobian
    if refel.dim == 1
        xr  = refel.Dg*pts[:];
        J = xr[:];
        rx = 1 ./ J;
        D = Jacobian(rx,[],[],[],[],[],[],[],[]);
        cachesim_load_range(10);
        cachesim_store_range(13);
        cachesim_store_range(16);
        
    elseif refel.dim == 2
        (xr, xs) = tensor_grad2(refel.Dg, pts[1,:][:]);
        (yr, ys) = tensor_grad2(refel.Dg, pts[2,:][:]);
        J = -xs.*yr + xr.*ys;
        
        rx =  ys./J;
        sx = -yr./J;
        ry = -xs./J;
        sy =  xr./J;
        D = Jacobian(rx,ry,[],sx,sy,[],[],[],[]);
        
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_store_range(13);
        cachesim_store_range(14);
        cachesim_store_range(16);
        
    else
        (xr, xs, xt) = tensor_grad3(refel.Dg, pts[1,:][:]);
        (yr, ys, yt) = tensor_grad3(refel.Dg, pts[2,:][:]);
        (zr, zs, zt) = tensor_grad3(refel.Dg, pts[3,:][:]);
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
        
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_load_range(12);
        cachesim_store_range(13);
        cachesim_store_range(14);
        cachesim_store_range(15);
        cachesim_store_range(16);
    end
    
    return (J,D);
end

function build_deriv_matrix_cachesim(refel, J)
    if refel.dim == 1
        RQ1 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j];
            end
        end
        return (RQ1,RD1);
        
        cachesim_load_range(13);
        cachesim_load_range(7);
        cachesim_load_range(10);
        
    elseif refel.dim == 2
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j];
            end
        end
        return (RQ1, RQ2, RD1, RD2);
        
        cachesim_load_range(7);
        cachesim_load_range(8);
        cachesim_load_range(9);
        cachesim_load_range(10);
        cachesim_load_range(11);
        cachesim_load_range(12);
        cachesim_load_range(13);
        cachesim_load_range(14);
        cachesim_load_range(15);
        
    elseif refel.dim == 3
        RQ1 = zeros(size(refel.Q));
        RQ2 = zeros(size(refel.Q));
        RQ3 = zeros(size(refel.Q));
        RD1 = zeros(size(refel.Q));
        RD2 = zeros(size(refel.Q));
        RD3 = zeros(size(refel.Q));
        for i=1:length(J.rx)
            for j=1:length(J.rx)
                RQ1[i,j] = J.rx[i]*refel.Qr[i,j] + J.sx[i]*refel.Qs[i,j] + J.tx[i]*refel.Qt[i,j];
                RQ2[i,j] = J.ry[i]*refel.Qr[i,j] + J.sy[i]*refel.Qs[i,j] + J.ty[i]*refel.Qt[i,j];
                RQ3[i,j] = J.rz[i]*refel.Qr[i,j] + J.sz[i]*refel.Qs[i,j] + J.tz[i]*refel.Qt[i,j];
                RD1[i,j] = J.rx[i]*refel.Ddr[i,j] + J.sx[i]*refel.Dds[i,j] + J.tx[i]*refel.Ddt[i,j];
                RD2[i,j] = J.ry[i]*refel.Ddr[i,j] + J.sy[i]*refel.Dds[i,j] + J.ty[i]*refel.Ddt[i,j];
                RD3[i,j] = J.rz[i]*refel.Ddr[i,j] + J.sz[i]*refel.Dds[i,j] + J.tz[i]*refel.Ddt[i,j];
            end
        end
        return (RQ1, RQ2, RQ3, RD1, RD2, RD3);
    end
end