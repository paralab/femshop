#=
# A coefficient to be used in the weak form.
# Can be a number or a generated function.
# Values are in an array to allow vector/tensor valued coefficients
=#

struct Coefficient
    symbol::Symbol          # symbol used in expressions
    symvar                  # SymType object
    index::Int              # index in the Femshop list of coefficients
    type::String            # constants for SCALAR, VECTOR, etc.
    location::String        # constant for NODAL, MODAL, CELL
    value::Array            # An array of either constant values(numbers) or genfunctions
end