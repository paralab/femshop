#=
# Represents bilinear and linear operators.
# Has a function that takes an element and refel and returns
# the elemental matrix or vector covering the degrees of freedom on that element
=#

mutable struct Bilinear
    func::GenFunction
end

mutable struct Linear
    func::GenFunction
end