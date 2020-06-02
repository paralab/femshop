#=
# A struct containing configuration information
# This is the most general Femshop info that should be (almost)problem indepenent.
# Premade configs can be loaded for different classes of problems.
# Values can be set through various macros or automatically by problem specification.
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_config
    # Domain
    dimension::Int          # 1,2,3,4
    geometry::String        # square, irregular
    mesh_type::String       # unstructured, tree

    # FEM details
    solver_type::String     # cg, dg, hdg
    trial_function::String  # Legendre, nodal, modal
    test_function::String   # same as above
    elemental_nodes::String # uniform, gauss, lobatto (higher order node distribution within elements)
    quadrature::String      # uniform, gauss, lobatto (similar to above)
    p_adaptive::Bool        # Do adaptive p-refinement?
    basis_order_min::Int    # minimum order to use in p-refinement, or if p_adaptive is false
    basis_order_max::Int    # maximum order
    
    # Other solver details
    linear::Bool            # Is the equation linear?
    t_adaptive::Bool        # Do adaptive t_refinement?
    stepper::String         # Euler-explicit/implicit, RK4, LSRK4, etc. Type of time stepper to use
    linalg_matrixfree::Bool # Use matrix free methods?
    linalg_backend::String  # ours, petsc, ?? (What to use for linear algebra)
    output_format::String   # VTK, raw, custom (format for storing solutions)
    
    # Constructor builds a default config.
    Femshop_config() = new(
        1,
        SQUARE,
        TREE,
        DG,
        LEGENDRE,
        LEGENDRE,
        LOBATTO,
        LOBATTO,
        false,
        4,
        4,
        true,
        false,
        LSRK4,
        false,
        OURS,
        VTK
    );
end
