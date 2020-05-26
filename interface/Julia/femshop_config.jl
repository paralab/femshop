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
    solver_type::String #cg, dg, hdg
    # songzhe: I think basis_type and trial_function/test_function should be consistnent
    # if basis_type is nodal, then trial/test may not be legendre
    # actually I think with trial/test, we don't need basis_type
    basis_type::String
    # songzhe: add some options
    trial_function::String #nodal, modal
    test_function::String #same as above
    # songzhe: I assume this is higher order node distribution in one element
    elemental_nodes::String # uniform, GLL
    # songzhe: add gp distribution, which can be the same as above or different
    quadrature::String # GL, GLL
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
    # songzhe: I added a femshop_prob struct if this is what you mean
    
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
