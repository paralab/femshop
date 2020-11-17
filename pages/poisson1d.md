---
title: Poisson 1D
---

## 1D Poisson

<img src="images/poisson1d.png" alt="poisson1d" width="400">

The Julia script: <a href="https://github.com/paralab/femshop/blob/master/femshop/examples/example-poisson1d.jl">example-poisson1d.jl<\a>

The 1D Poisson equation with Dirichlet boundary and smooth functions is about as simple as it gets. This example demonstrates the basics of setting up a problem in Fenshop. A uniform discretization of the unit domain is used with p=4 polynomial basis function space.

Begin by importing and using the Femshop module. Then initialize. The name here is only used when generating code files.
```
include("../Femshop.jl");
using .Femshop
init_femshop("poisson1d");
```
Then set up the configuration. This example simply sets dimensionality of the domain and basis function space.
```
@domain(1)                  # dimension
@functionSpace(LEGENDRE, 4) # basis function, order
```
Use the built-in simple mesh generator to make the mesh and set up all node mappings.
```
@mesh(LINEMESH, 20)         # uniform 1D mesh. 2nd arg is # of elements
```
Define the variable, test function, and forcing function symbols.
```
@variable(u)                # same as @variable(u, SCALAR)
@testSymbol(v)              # sets the symbol for a SCALAR test function
@coefficient(f, "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")
```
Convert the PDE
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=\begin{align}\Delta%20u&=f(x)\\ u(0)&=u(1)=0\end{align}"> <\div>
into the weak form
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=(grad(u),grad(v))=(f,v)"> <\div>

The boundary condition is specified.
```
@boundary(u, 1, DIRICHLET, "0") # boundary condition for BID 1 is Dirichlet with value 0
```
Then write the weak form expression by setting the weak form equation above equal to zero. Finally, solve for u.
```
@weakForm(u, "-grad(u)*grad(v) - f*v")
solve(u);
```
End things with `@finalize()` to finish up any generated files and the log.