
# External modules
using SparseArrays
using LinearAlgebra
# When we need SymPy
#=
######################
# NOTE: This is not a long term solution.
# Once the package is set up, we can put this dependency in the .toml
######################
try
    using SymPy
catch e
    println("Julia SymPy is not installed. Shall I install for you?[y/n]");
    if readline()[1] == 'y'
        println("Alright. This could take a minute.");
        using Pkg
        Pkg.update();
        Pkg.add("Conda")
        using Conda
        Conda.update()
        Pkg.add("SymPy")
        using SymPy
    else
        println("Proceeding without SymPy (Which is fine. We don't use it yet.)");
    end
end
=#

# include these first
include("femshop_constants.jl");
include("femshop_config.jl");
include("femshop_prob.jl");

# include these next
include("logging.jl");
include("mesh_read.jl");
include("mesh_write.jl")
include("function_utils.jl");
include("variables.jl");

# Femshop submodules
include("CodeGenerator.jl");
using .CodeGenerator
include("DGSolver.jl");
using .DGSolver

# include these last (depend on submodules)
include("rhs_params.jl");
