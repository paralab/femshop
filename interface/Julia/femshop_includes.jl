
# External modules
using SparseArrays
using LinearAlgebra

# ######################
# # NOTE: This is not a long term solution.
# # Once the package is set up, we can put this dependency in the .toml
# try
#     using SymEngine
# catch e
#     println("Julia SymEngine is not yet installed. Installing now.");
#     using Pkg
#     Pkg.add("SymEngine")
#     using SymEngine
# end
# ######################

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
include("coefficient.jl");
include("bilinear_linear.jl");

# Femshop submodules
include("SymbolicParser.jl")
using .SymbolicParser
include("CodeGenerator.jl");
using .CodeGenerator
include("DGSolver.jl");
using .DGSolver

# include these last (depend on submodules)
include("rhs_params.jl");
