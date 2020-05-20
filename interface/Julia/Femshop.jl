#=
The main module for Femshop.
We can reorganize things and make submodules as desired.
=#

module Femshop

if !@isdefined(femshop_modules_loaded)
    femshop_modules_loaded = 1;
    # External modules
    using SparseArrays
    using LinearAlgebra

    # Femshop submodules
    include("CodeGenerator.jl");
    using .CodeGenerator

    # Directly included files
    include("femshop_constants.jl");
    include("femshop_config.jl");
end

# Public macros and functions
export @language, @domain, @solver, @functionSpace, @trialFunction,
        @testFunction, @nodes, @order, @finalize
export init_femshop, set_language, finalize

# Init global variables
# Optional name is the project name used for output
function init_femshop(name="")
    global config = Femshop_config();
    global project_name = name;
    global outputDir = pwd();
    global genfiles = nothing;
end

init_femshop("");

include("macros.jl"); # included here after globals are defined

function set_language(lang, dirpath, filename, head="")
    global outputDir = dirpath;
    global genfiles = CodeGenerator.init_codegenerator(lang, dirpath, filename, head);
end

function finalize()
    # Finalize generation
    if genfiles != nothing
        finalize_codegenerator();
    end
    # anything else
    println("Femshop has completed.");
end


end # module
