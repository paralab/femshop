#=
# 2D NS eq. Dirichlet bc, CG
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("NS");

# Try making an optional log
@useLog("NSlog")

# Set up the configuration (order doesn't matter)
@domain(2, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 2)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(EULER_IMPLICIT)            # time stepper (optional second arg is CFL#)

# Specify the problem
num_elem = 32;
@mesh(QUADMESH, num_elem, 4)                   # .msh file or generate our own

@variable(u)                        # same as @variable(u, SCALAR)
@variable(v)                        # same as @variable(u, SCALAR)
@variable(uold)                        # same as @variable(u, SCALAR)
@variable(vold)                        # same as @variable(u, SCALAR)
@variable(p)                        # same as @variable(u, SCALAR)
@variable(du)                        # same as @variable(u, SCALAR)
@variable(dv)                        # same as @variable(u, SCALAR)
@variable(dp)                        # same as @variable(u, SCALAR)

@testSymbol(w)                    # sets the symbol for a test function

T = 1
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
#@initial(u, "y > 0.9 ? 1 : 0")  # initial condition needed if time dependent
#@initial(uold, "y > 0.9 ? 1 : 0")  # initial condition needed if time dependent
@initial(u, "y")  # initial condition needed if time dependent
@initial(uold, "y")  # initial condition needed if time dependent
@initial(du, "0")  # initial condition needed if time dependent
@initial(v, "0")  # initial condition needed if time dependent
@initial(vold, "0")  # initial condition needed if time dependent
@initial(dv, "0")  # initial condition needed if time dependent
@initial(p, "0")  # initial condition needed if time dependent
@initial(dp, "0")  # initial condition needed if time dependent

#@boundary(u, 1, DIRICHLET, 0)
@boundary(du, 1, DIRICHLET, 0)
#@boundary(v, 1, DIRICHLET, 0)
@boundary(dv, 1, DIRICHLET, 0)
#@boundary(p, 1, DIRICHLET, 0)
#@boundary(dp, 1, DIRICHLET, 0)
#@boundary(u, 2, DIRICHLET, 0)
@boundary(du, 2, DIRICHLET, 0)
#@boundary(v, 2, DIRICHLET, 0)
@boundary(dv, 2, DIRICHLET, 0)
#@boundary(p, 2, DIRICHLET, 0)
#@boundary(dp, 2, DIRICHLET, 0)
#@boundary(u, 3, DIRICHLET, 0)
@boundary(du, 3, DIRICHLET, 0)
#@boundary(v, 3, DIRICHLET, 0)
@boundary(dv, 3, DIRICHLET, 0)
#@boundary(p, 3, DIRICHLET, 0)
#@boundary(dp, 3, DIRICHLET, 0)
#@boundary(u, 4, DIRICHLET, 1)
@boundary(du, 4, DIRICHLET, 0)
#@boundary(v, 4, DIRICHLET, 0)
@boundary(dv, 4, DIRICHLET, 0)
#@boundary(p, 4, DIRICHLET, 0)
#@boundary(dp, 4, DIRICHLET, 0)


# Write the weak form
@coefficient(mu, 0.01)
@coefficient(dtc, 0.01)
@coefficient(h, 1.0/32)
@coefficient(coe1, 4.0)
@coefficient(coe2, 36.0)

@parameter(tauM, "1.0 ./ (coe1 ./ dtc ./ dtc+ (u*u+v*v) .^ 0.5 ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")
@parameter(tauC, "h*h* (coe1 ./ dtc ./ dtc+ (u*u+v*v) .^ 0.5 ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")
#@parameter(tauC, "(h.^2 ./ tauM")

#@weakForm([du, dv], ["w*Dt(du)", "w*Dt(dv)"])

#@weakForm([du, dv, dp], ["w*(Dt(du) + (u*deriv(du,1)+v*deriv(du,2))) - deriv(w,1)*dp + mu*dot(grad(w), grad(du)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(Dt(du) + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1)) - (w*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2))) - deriv(w,1)*p + mu*dot(grad(w), grad(u)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1)))", 
@weakForm([du, dv, dp], ["w*(Dt(du)) - (w*u)",
                         "w*(Dt(dv) + (u*deriv(dv,1)+v*deriv(dv,2))) - deriv(w,2)*dp + mu*dot(grad(w), grad(dv)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*( Dt(dv) + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2) ) - (w*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2))) - deriv(w,2)*p + mu*dot(grad(w), grad(v)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2) ))", 
                         "w*(deriv(du,1)+deriv(dv,2)) + tauC*(deriv(w,1)*( (u-uold) ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1) ) + deriv(w,2)*( (v-vold) ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2) )) - (w*(deriv(u,1)+deriv(v,2)) + tauC*(deriv(w,1)*( (u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1) ) + deriv(w,2)*( (v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2) )))"])

solve([du, dv, dp], [u, v, p, uold, vold], nonlinear=true);

#analytical
#=
x = collect(1:1/num_elem/order:2);
u_a = zeros(length(x));
for i = 1:length(x)
    u_a[i] = x[i]*(1+1.0/3*x[i]^2)^(-1/2);
end

#check integral error
err_l2 = 0.0;
for i = 1:length(x)-1
    err_a = ((u.values[i]-u_a[i])^2+(u.values[i+1]-u_a[i+1])^2)/2.0*(1.0/num_elem/order);
    global err_l2 += err_a;
end

print("L2 error = ", err_l2, "\n")
=#
# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes[:,1], Femshop.grid_data.allnodes[:,2], u.values, st = :surface))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
