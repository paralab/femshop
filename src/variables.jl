#=
# Variable struct and related functions
=#

mutable struct Variable
    symbol                  # symbol used in expressions referring to this variable
    symvar                  # SymType object
    index::Int              # index in the Femshop list of variables
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL, MODAL, CELL
    values::Array{Float64}  # a C x N array, C is number of components (SCALAR=1, etc.)
    
    # not currently used, may be removed
    dependson               # a list of variables that this depends on
    ready::Bool             # Is this variable's values ready? Can dependent variables use it?
end
