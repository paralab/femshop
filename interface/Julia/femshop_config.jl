#=
A struct containing configuration information
=#
if !@isdefined(IRREGULAR)
    include("femshop_constants.jl");
end

mutable struct Femshop_config
    # Domain
    dimension::Int
    geometry::String
    mesh_type::String

    # FEM details
    solver_type::String
    basis_type::String
    trial_function::String
    test_function::String
    elemental_nodes::String
    p_adaptive::Bool
    basis_order_min::Int
    basis_order_max::Int

    # Other solver details
    # These are just some examples, many more are needed.
    linear::Bool
    t_adaptive::Bool
    stepper::String
    linalg_matrixfree::Bool
    linalg_backend::String
    output_format::String

    # Variable and boundary things should be in the problem
    # specification rather than here.

    # Constructor builds a default config.
    Femshop_config() = new(
        1,
        SQUARE,
        TREE,
        DG,
        NODAL,
        LEGENDRE,
        LEGENDRE,
        LOBATTO,
        false,
        4,
        5,
        true,
        false,
        RK4,
        false,
        OURS,
        VTK
    );
end
