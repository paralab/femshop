using Femshop
using Test

@testset "Femshop.jl" begin
    init_femshop("poisson1d");

    useLog("poisson1dlog")

    domain(1)                      # dimension
    functionSpace(order=3)         # basis function polynomial order

    mesh(LINEMESH, elsperdim=20)   # build uniform LINEMESH with 20 elements

    u = variable("u")              # make a scalar variable with symbol u
    testSymbol("v")                # sets the symbol for a test function

    boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

    coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
    weakForm(u, "-grad(u)*grad(v) - f*v")

    solve(u);

    maxerr = 0;
    exact(x) = sin(10*pi*x)*sin(pi*x);

    for i=1:size(Femshop.grid_data.allnodes,2)
        x = Femshop.grid_data.allnodes[1,i];
        err = abs(u.values[i] - exact(x));
        global maxerr;
        maxerr = max(err,maxerr);
    end
    println("max error = "*string(maxerr));
    
    @test(maxerr < 0.01)
    
    # Write your tests here.
end
