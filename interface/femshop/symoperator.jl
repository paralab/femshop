#=
Operators that work on SymType objects
=#

# An abstract operator type.
struct SymOperator
    symbol::Symbol      # The name used in the input expression
    op                  # Function handle for the operator
end
