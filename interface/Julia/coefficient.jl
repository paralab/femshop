#=
# A coefficient to be used in the weak form.
# Can be a number or a generated function.
# Values are in an array to allow vector/tensor valued coefficients
=#

struct Coefficient
    symbol::Symbol
    value::Array
end