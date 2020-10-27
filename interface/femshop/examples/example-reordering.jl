#=
# Test out different node orderings
=#
if !@isdefined(Femshop)
    include("../Femshop.jl");
    using .Femshop
end
n = 31;
ord = 1;
gd = n*ord + 1; # grid size in each dimension
griddim = [gd,gd,gd];
te = 1; # elemental loop order is also tiled with this width
tds = [ord*te+1]; # tile size in each dimension

init_femshop("reordering");
@useLog("reorderinglog")

#cachesim(true);

@domain(3)
@functionSpace(LEGENDRE, ord)

@matrixFree(200,1e-6)

@mesh(HEXMESH, n)

@variable(u)
@variable(q1)
@variable(q2)
@variable(q3)
@variable(q4)
@variable(q5)
@variable(q6)

@testSymbol(v)

@boundary(u, 1, DIRICHLET, 0)
@boundary(q1, 1, DIRICHLET, 0)
@boundary(q2, 1, DIRICHLET, 0)
@boundary(q3, 1, DIRICHLET, 0)
@boundary(q4, 1, DIRICHLET, 0)
@boundary(q5, 1, DIRICHLET, 0)
@boundary(q6, 1, DIRICHLET, 0)

@coefficient(f, "-3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm(u, "-dot(grad(u), grad(v)) - f*v")

@coefficient(g, "(-3*pi*pi + 5)*sin(pi*x)*sin(pi*y)*sin(pi*z)")
@weakForm([q1,q2,q3,q4,q5,q6], ["-dot(grad(q1), grad(v)) + q2*v + q3*v + q4*v + q5*v + q6*v - g*v",
                                "-dot(grad(q2), grad(v)) + q3*v + q4*v + q5*v + q6*v + q1*v - g*v",
                                "-dot(grad(q3), grad(v)) + q4*v + q5*v + q6*v + q1*v + q2*v - g*v",
                                "-dot(grad(q4), grad(v)) + q5*v + q6*v + q1*v + q2*v + q3*v - g*v",
                                "-dot(grad(q5), grad(v)) + q6*v + q1*v + q2*v + q3*v + q4*v - g*v",
                                "-dot(grad(q6), grad(v)) + q1*v + q2*v + q3*v + q4*v + q5*v - g*v"])

#
#using Plots
#pyplot();

# First do u (one DOF)
times = 1; # average times times
println("One DOF per node, average of "*string(times)*" runs")
#warm up
solve(u);

# Lex. ordering
regtime = Base.Libc.time();
for iter=1:times
    solve(u);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
println("Lex. time = "*string(regtime)*"sec.");
#display(spy(Femshop.CGSolver.Amat, legend=nothing, reuse=false));

# Tiled ordering
for ti=1:length(tds)
    td = tds[ti];
    @mesh(HEXMESH, n)
    tiled_nodes(griddim,(td,td,td));
    #tiled_elements((n,n,n),(te,te,te));
    
    #solve(u);
    tiletime = Base.Libc.time();
    for iter=1:times
        solve(u);
    end
    tiletime = Base.Libc.time() - tiletime;
    tiletime /= times;

    println("tiled("*string(td)*") = "*string(tiletime)*"sec.");
end
#display(spy(Femshop.CGSolver.Amat, legend=nothing, reuse=false));

# morton ordering
@mesh(HEXMESH, n)
morton_nodes(griddim);
morton_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;

println("morton = "*string(tiletime)*"sec.");
#display(spy(Femshop.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
@mesh(HEXMESH, n)
hilbert_nodes(griddim);
hilbert_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve(u);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;

println("hilbert = "*string(tiletime)*"sec.");

###############################################################################
# then do qn (6 DOF)
times = 1;
println("Six DOF per node, averaged "*string(times)*" times")

@mesh(HEXMESH, n)

#warm up
solve([q1,q2,q3,q4,q5,q6]);
times = 1; # average times times

# Lex. ordering
regtime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
regtime = Base.Libc.time() - regtime;
regtime /= times;
println("Lex. time = "*string(regtime)*"sec.");

# Tiled ordering
for ti=1:length(tds)
    td = tds[ti];
    @mesh(HEXMESH, n)
    tiled_nodes((gd,gd,gd),(td,td,td));
    #tiled_elements((n,n,n),(2,2,2));
    
    #solve([q1,q2,q3,q4,q5,q6]);
    tiletime = Base.Libc.time();
    for iter=1:times
        solve([q1,q2,q3,q4,q5,q6]);
    end
    tiletime = Base.Libc.time() - tiletime;
    tiletime /= times;

    println("tiled("*string(td)*") = "*string(tiletime)*"sec.");
end

# morton ordering
@mesh(HEXMESH, n)
morton_nodes(griddim);
morton_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;

println("morton = "*string(tiletime)*"sec.");
#display(spy(Femshop.CGSolver.Amat, legend=nothing, reuse=false));

# hilbert ordering
@mesh(HEXMESH, n)
hilbert_nodes(griddim);
hilbert_elements([n,n,n]);

#solve(u);
tiletime = Base.Libc.time();
for iter=1:times
    solve([q1,q2,q3,q4,q5,q6]);
end
tiletime = Base.Libc.time() - tiletime;
tiletime /= times;

println("hilbert = "*string(tiletime)*"sec.");

# exact solution is sin(pi*x)*sin(pi*y)*sin(pi*z)
# check error
maxerr = 0;
maxerrq = 0;
exact(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);
exactq(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z);

for i=1:size(Femshop.grid_data.allnodes,1)
    x = Femshop.grid_data.allnodes[i,1];
    y = Femshop.grid_data.allnodes[i,2];
    z = Femshop.grid_data.allnodes[i,3];
    err = abs(u.values[i] - exact(x,y,z));
    errq = abs(q1.values[i] - exactq(x,y,z)) + abs(q2.values[i] - exactq(x,y,z)) + abs(q3.values[i] - exactq(x,y,z)) + 
            abs(q4.values[i] - exactq(x,y,z)) + abs(q5.values[i] - exactq(x,y,z)) + abs(q6.values[i] - exactq(x,y,z));
    global maxerr;
    global maxerrq;
    maxerr = max(err,maxerr);
    maxerrq = max(errq,maxerrq);
end
println("max error for u = "*string(maxerr));
println("max error for q = "*string(maxerrq));

@finalize()
