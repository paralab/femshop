#=
# 1D nonlinear heat eq. Dirichlet bc, CG
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
init_femshop("nonlinear_heat");

# Try making an optional log
@useLog("nonlinear_heatlog")

num_elem = 32;
order = 2;

# Set up the configuration (order doesn't matter)
@domain(1, SQUARE, UNSTRUCTURED)    # dimension, geometry, decomposition
@solver(CG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, order)     # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement
@stepper(EULER_IMPLICIT)            # time stepper (optional second arg is CFL#)

# Specify the problem
@mesh(LINEMESH, num_elem, 2, [1,2]) # type, num of elements, num of BIDs, interval

@variable(u)                        # same as @variable(u, SCALAR)
@variable(du)                       # same as @variable(du, SCALAR)

@testSymbol(v)                    # sets the symbol for a test function

T = 0.1;
@timeInterval(T)                    # (start, end) using this sets problem to time dependent
@initial(u, "(1.3093073414159544-0.8660254037844387)*x-1.3093073414159544+0.8660254037844387*2")  # initial condition needed if time dependent
@initial(du, "0")  # initial condition needed if time dependent

@boundary(du, 1, DIRICHLET, 0)
@boundary(du, 2, DIRICHLET, 0)
@boundary(u, 1, DIRICHLET, 0.866)
@boundary(u, 2, DIRICHLET, 1.31)

# Write the weak form
@coefficient(f, "x^(-4)")
@weakForm(du, "(grad(v)*grad(du) - v*5*u*u*u*u*f*du) - (grad(v)*grad(u) - v*u*u*u*u*u*f)")

solve(du,u, nonlinear=true);

#analytical
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
# solution is stored in the variable's "values"
#using Plots
#pyplot();
#display(plot(Femshop.grid_data.allnodes[1,:], Femshop.grid_data.allnodes[2,:], u.values[:], st = :surface))

# check
log_dump_config(Femshop.config);
log_dump_prob(Femshop.prob);

@finalize()
