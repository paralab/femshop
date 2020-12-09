---
title: Heat
---

## Heat

<img src="images/heat.png" alt="heat" width="400">

The Julia script: <a href="https://github.com/paralab/femshop/blob/master/femshop/examples/example-heat2d.jl">example-heat2d.jl</a>

A 2D heat equation demonstrates support for time dependent problems.

Begin by importing and using the Femshop module. Then initialize. The name here is only used when generating code files.
```
include("../Femshop.jl");
using .Femshop
init_femshop("heat");
```
Then set up the configuration. This example simply sets dimensionality of the domain and basis function space.
```
@domain(2)                  # dimension
@functionSpace(LEGENDRE, 4) # basis function, order
```
Use the built-in simple mesh generator to make the mesh and set up all node mappings.
```
@mesh(QUADMESH, 10) # has 10*10 uniform, square elements
```
Define the variable, test function, and coefficient symbols.
```
@variable(u)
@testSymbol(v)

@coefficient(f, "0.5*sin(6*pi*x)*sin(6*pi*y)")
```
Set up the time stepper and initial conditions. This example uses a low-storage RK4. Other explicit or implicit methods are available.
```
@stepper(LSRK4)  # Low-storage RK4
@timeInterval(1) # The end time
@initial(u, "abs(x-0.5)+abs(y-0.5) < 0.2 ? 1 : 0") # initial condition
```
Convert the PDE
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=\frac{d}{dt}u+D\Delta%20u=f"> </div>
into the weak form
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=\frac{d}{dt}(u,v)+D(\nabla%20u,\nabla%20v)=(f,v)"> </div>

The boundary conditions are specified.
```
@boundary(u, 1, DIRICHLET, 0)
```
Then write the weak form expression by setting the weak form equation above equal to zero. Finally, solve for u.
```
@weakForm(u, "Dt(u*v) + 0.01 * dot(grad(u),grad(v)) - f*v")
solve(u);
```
End things with `@finalize()` to finish up any generated files and the log.