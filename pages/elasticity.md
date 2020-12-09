---
title: Linear Elasticity
---

## Linear Elasticity

<img src="images/elasticity.png" alt="elasticity" width="400">

The Julia script: <a href="https://github.com/paralab/femshop/blob/master/femshop/examples/example-elasticity.jl">example-elasticity.jl</a>

A 3D linear elasticity equation that models the gravity induced deflection of a beam that is fixed at one end. This example demonstrates the use of vector valued variables, coefficients, and function spaces. It also introduces mixed Dirichlet and Neumann boundary conditions.

Begin by importing and using the Femshop module. Then initialize. The name here is only used when generating code files.
```
include("../Femshop.jl");
using .Femshop
init_femshop("elasticity");
```
Then set up the configuration. This example simply sets dimensionality of the domain and basis function space.
```
@domain(3)                  # dimension
@functionSpace(LEGENDRE, 4) # basis function, order
```
Use the built-in simple mesh generator to make the mesh and set up all node mappings.
```
n = [10,4,4]; # The numbers of elements in each dimension. 10x4x4 elements
interval = [0,1,0,0.2,0,0.2]; # The limits of the domain. A narrow beam.

# 3D mesh defined by n and interval above with two boundary regions.
@mesh(HEXMESH, n, 2, interval) 
```
Define the variable, test function, and coefficient symbols.
```
@variable(u, VECTOR)
@testSymbol(v, VECTOR)

@coefficient(mu, "x>0.5 ? 0.2 : 10") # discontinuous mu
@coefficient(lambda, 1.25)
@coefficient(f, VECTOR, ["0","0","-0.1"]) # gravitational force
```
Convert the PDE
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=u"> </div>
into the weak form
<div align="center"><img src="https://render.githubusercontent.com/render/math?math=((\lambda%20\Nabla%20\cdot%20uI+\mu%20(\Nabla%20u+\Nabla%20u^{T})):\Nabla%20v))=(f\cdot%20v)"> </div>

The boundary conditions are specified.
```
@boundary(u, 1, DIRICHLET, [0,0,0]) # x=0
@boundary(u, 2, NEUMANN, [0,0,0])   # other
```
Then write the weak form expression by setting the weak form equation above equal to zero. Finally, solve for u.
```
@weakForm(u, "inner( (lambda * div(u) .* [1 0 0; 0 1 0; 0 0 1] + mu .* (grad(u) + transpose(grad(u)))), grad(v)) - dot(f,v)")
solve(u);
```
End things with `@finalize()` to finish up any generated files and the log.