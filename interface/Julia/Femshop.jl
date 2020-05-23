#=
The main module for Femshop.
We can reorganize things and make submodules as desired.
=#
module Femshop

include("femshop_includes.jl");

# Public macros and functions
export @language, @domain, @mesh, @solver, @functionSpace, @trialFunction,
        @testFunction, @nodes, @order, @useLog, @finalize
export init_femshop, set_language, add_mesh, finalize

# External modules
using SparseArrays
using LinearAlgebra
# Femshop submodules
using .CodeGenerator

# Module's global variables
config = Femshop_config();
project_name = "";
output_dir = pwd();
mesh_data = nothing;
gen_files = nothing;

include("macros.jl"); # included here after globals are defined

function set_language(lang, dirpath, name, head="")
    global output_dir = dirpath;
    global project_name = name;
    global gen_files = CodeGenerator.init_codegenerator(lang, dirpath, name, head);
end

function add_mesh(mesh)
    global mesh_data = mesh;
    log_entry("Added mesh with "*string(mesh.nx)*" nodes and "*string(mesh.nel)*" elements.");
end

function finalize()
    # Finalize generation
    if gen_files != nothing
        finalize_codegenerator();
    end
    # anything else
    close_log();
    println("Femshop has completed.");
end


end # module
