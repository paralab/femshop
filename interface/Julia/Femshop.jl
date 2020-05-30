#=
The main module for Femshop.
We can reorganize things and make submodules as desired.
=#
module Femshop

# Public macros and functions
export @language, @domain, @mesh, @solver, @functionSpace, @trialFunction,
        @testFunction, @nodes, @order, @boundary, @variable, @initial,
        @timeInterval,
        @outputMesh, @useLog, @finalize
export init_femshop, set_language, add_mesh, output_mesh, add_initial_condition, solve, finalize

#include("femshop_includes.jl");

# Module's global variables
config = nothing;
prob = nothing;
project_name = "";
output_dir = pwd();
mesh_data = nothing;
gen_files = nothing;

include("femshop_includes.jl");
include("macros.jl"); # included here after globals are defined

config = Femshop_config();
prob = Femshop_prob();

function set_language(lang, dirpath, name, head="")
    global output_dir = dirpath;
    global project_name = name;
    global gen_files = CodeGenerator.init_codegenerator(lang, dirpath, name, head);
end

function add_mesh(mesh)
    global mesh_data = mesh;
    log_entry("Added mesh with "*string(mesh.nx)*" nodes and "*string(mesh.nel)*" elements.");
end

function output_mesh(file, format)
    write_mesh(file, format, mesh_data);
    log_entry("Wrote mesh data to file: "*file*".msh");
end

function add_initial_condition(varindex)
    while length(prob.initial) < varindex
        prob.initial = [prob.initial; nothing];
    end
    prob.initial[varindex] = genfunctions[end];
    # hold off on initializing till solve or generate is determined.
end

function solve()
    # TEMPORARY ===============================
    init_dgsolver();
    DGSolver.solve();
    # =========================================
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
