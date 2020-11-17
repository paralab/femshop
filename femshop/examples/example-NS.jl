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
@domain(2, SQUARE, UNIFORM_GRID)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 1)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(EULER_IMPLICIT)            # time stepper (optional second arg is CFL#)
dt = 0.05;
nsteps = 50;
@setSteps(dt, nsteps)                # manually set stepper dt and Nsteps, overriding defaults and interval

# Specify the problem
num_elem = 32;
@mesh(QUADMESH, num_elem, 4)

@variable(u)                        # same as @variable(u, SCALAR)
@variable(v)                        # same as @variable(u, SCALAR)
@variable(uold)                        # same as @variable(u, SCALAR)
@variable(vold)                        # same as @variable(u, SCALAR)
@variable(p)                        # same as @variable(u, SCALAR)
@variable(du)                        # same as @variable(u, SCALAR)
@variable(dv)                        # same as @variable(u, SCALAR)
@variable(dp)                        # same as @variable(u, SCALAR)

@testSymbol(w)                    # sets the symbol for a test function

#T = 2 # set manually using setSteps above
#@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "y > 0.99 ? 1 : 0")  # initial condition needed if time dependent
@initial(uold, "y > 0.99 ? 1 : 0")
@initial(du, "0")
@initial(v, "0")
@initial(vold, "0")
@initial(dv, "0")
@initial(p, "0")
@initial(dp, "0")

@boundary(du, 1, DIRICHLET, 0)
@boundary(du, 2, DIRICHLET, 0)
@boundary(du, 3, DIRICHLET, 0)
@boundary(du, 4, DIRICHLET, 0)

@boundary(dv, 1, DIRICHLET, 0)
@boundary(dv, 2, DIRICHLET, 0)
@boundary(dv, 3, DIRICHLET, 0)
@boundary(dv, 4, DIRICHLET, 0)

@boundary(dp, 1, NO_BC)
@boundary(dp, 2, NO_BC)
@boundary(dp, 3, NO_BC)
@boundary(dp, 4, NO_BC)
@referencePoint(dp, [0,0], 0)


# Write the weak form
@coefficient(mu, 0.01)
@coefficient(dtc, dt)
@coefficient(h, 1.0 / 32)
@coefficient(coe1, 4.0)
@coefficient(coe2, 6.0)

@parameter(tauM, "1.0 ./ (coe1 ./ dtc ./ dtc+ (u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")
@parameter(tauC, "0.1*h*h* ((u*u+v*v) ./ h ./ h+coe2*mu*mu ./ h ./ h ./ h ./ h) .^ 0.5")

@weakForm([du, dv, dp], ["w*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2) + deriv(u,2)*dv)) - deriv(w,1)*dp + mu*dot(grad(w), grad(du)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1)) - (w*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2))) - deriv(w,1)*p + mu*dot(grad(w), grad(u)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1)))", 
                         "w*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2) + deriv(v,1)*du)) - deriv(w,2)*dp + mu*dot(grad(w), grad(dv)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*(dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2)) - (w*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2))) - deriv(w,2)*p + mu*dot(grad(w), grad(v)) + tauM*(u*deriv(w,1)+v*deriv(w,2))*((v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2)))", 
	                     "w*(deriv(du,1)+deriv(dv,2)) + tauC*(deriv(w,1)*( du ./ dtc + (u*deriv(du,1)+v*deriv(du,2)) + deriv(dp,1) ) + deriv(w,2)*( dv ./ dtc + (u*deriv(dv,1)+v*deriv(dv,2)) + deriv(dp,2) )) - (w*(deriv(u,1)+deriv(v,2)) + tauC*(deriv(w,1)*( (u-uold) ./ dtc + (u*deriv(u,1)+v*deriv(u,2)) + deriv(p,1) ) + deriv(w,2)*( (v-vold) ./ dtc + (u*deriv(v,1)+v*deriv(v,2)) + deriv(p,2) )))"])
						  
solve([du, dv, dp], [u, v, p, uold, vold], nonlinear=true);

# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes[1,:], Femshop.grid_data.allnodes[2,:], u.values[:], st = :surface, reuse=false))
#display(plot(Femshop.grid_data.allnodes[1,:], Femshop.grid_data.allnodes[2,:], v.values[:], st = :surface, reuse=false))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
