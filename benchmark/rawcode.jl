# Find a way to run these example in isolation
init_femshop("poisson1d")

useLog("poisson1dlog")

domain(1)

functionSpace(order=3)
mesh(LINEMESH, elsperdim=20)

u = variable("u")
testSymbol("v")
boundary(u, 1, DIRICHLET, 0)

coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) -
             pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
weakForm(u, "-grad(u)*grad(v) - f*v")

SUITE["poisson1d-assemble"] = @benchmarkable Femshop.CGSolver.assemble($u,
                         $Femshop.bilinears[u.index], $Femshop.linears[u.index])