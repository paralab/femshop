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

### Module's global variables ###
# config
config = nothing;
prob = nothing;
project_name = "unnamedProject";
output_dir = pwd();
gen_files = nothing;
#log
use_log = false;
log_file = nothing;
log_line_index = 1;
#mesh
mesh_data = nothing;
#problem variables
var_count = 0;
variables = [];
#generated functions
genfunc_count = 0;
genfunctions = [];

include("femshop_includes.jl");
include("macros.jl"); # included here after globals are defined

config = Femshop_config(); # These need to be initialized here
prob = Femshop_prob();

function init_femshop(name="unnamedProject")
    global project_name = name;
    global gen_files = nothing;
    if log_file != nothing
        close(log_file);
        use_log = false;
        log_line_index = 1;
    end
    global mesh_data = nothing;
    global var_count = 0;
    global variables = [];
    global genfun_count = 0;
    global genfunctions = [];
end

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
    lhs = nothing;
    lhspars = nothing;
    rhs = nothing;
    rhspars = nothing;
    if config.solver_type == DG
        t = @elapsed(init_dgsolver());
        log_entry("Set up DG solver.(took "*string(t)*" seconds)");
        if config.dimension == 1
            if prob.time_dependent
                rhs = DGSolver.rhs_dg_1d;
                rhspars = build_rhs_params();
            else
                # TODO
            end
        else
            # TODO
        end
        
        DGSolver.solve(lhs, lhspars, rhs, rhspars);
    else
        # TODO
    end
    
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
