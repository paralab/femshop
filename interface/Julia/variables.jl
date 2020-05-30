#=
# Variable struct and related functions
=#
export Variable, add_variable

#TODO move these to Femshop.jl
var_count = 0;
variables = [];

mutable struct Variable
    symbol                  # symbol used in expressions referring to this variable
    index::Int              # index in the Femshop list of variables
    type::String            # constants for SCALAR, VECTOR, etc.
    values::Array{Float64}  # an N x C array, C is number of components (SCALAR=1, etc.)
    dependson               # a list of variables that this depends on
    ready::Bool             # Is this variable's values ready? Can dependent variables use it?
end

function add_variable(var)
    global var_count += 1;
    # adjust values arrays
    if var.type == SCALAR
        var.values = zeros(prob.mesh_dofs);
    elseif var.type == VECTOR
        var.values = zeros(prob.mesh_dofs, config.dimension);
    elseif var.type == TENSOR
        var.values = zeros(prob.mesh_dofs, config.dimension*config.dimension);
    end
    global variables = [variables; var];
    
    log_entry("Added variable: "*string(var.symbol)*" of type: "*var.type);
end