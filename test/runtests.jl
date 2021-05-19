using Femshop
using Test

@testset "Femshop.jl" begin

    @useLog("poisson1dlog")
    @domain(1)                      # dimension
    @functionSpace(LEGENDRE, 3)     # basis function, order
    
    @mesh(LINEMESH, 20)             # build uniform LINEMESH. 2nd arg=# of elements, (optional)3rd arg=# of BIDs
    
    @variable(u)                    # same as @variable(u, SCALAR)
    @testSymbol(v)                  # sets the symbol for a test function
    
    @boundary(u, 1, DIRICHLET, "0")   # boundary condition for BID 1 is Dirichlet with value 0
    
    @coefficient(f, "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
    @weakForm(u, "-grad(u)*grad(v) - f*v")
    
    solve(u);
    
    # exact solution is sin(10*pi*x)*sin(pi*x)
    # check error
    maxerr = 0;
    exact(x) = sin(10*pi*x)*sin(pi*x);

    for i=1:size(Femshop.grid_data.allnodes,2)
        x = Femshop.grid_data.allnodes[1,i];
        err = abs(u.values[i] - exact(x));
        maxerr = max(err,maxerr);
    end
    println(maxerr)
    
    @test(maxerr < 0.01)
    
    # Write your tests here.
end
