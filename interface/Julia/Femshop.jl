#=
The main module for Femshop.
We can reorganize things and make submodules as desired.
=#
module Femshop

# Public macros and functions
export @language, @domain, @mesh, @solver, @stepper, @functionSpace, @trialFunction, @matrixFree,
        @testFunction, @nodes, @order, @boundary, @variable, @coefficient, @parameter, @testSymbol, @initial,
        @timeInterval, @weakForm, @LHS, @RHS, @customOperator, @customOperatorFile,
        @outputMesh, @useLog, @finalize
export init_femshop, set_language, dendro, set_solver, set_stepper, set_matrix_free, reformat_for_stepper, 
        add_mesh, output_mesh, add_test_function, 
        add_initial_condition, add_boundary_condition, set_rhs, set_lhs, solve, finalize, cachesim, cachesim_solve, 
        morton_nodes, tiled_nodes, tiled_elements
export build_cache_level, build_cache
export sp_parse
export generate_code_layer
export Variable, add_variable
export Coefficient, add_coefficient
export Parameter, add_parameter

### Module's global variables ###
# config
config = nothing;
prob = nothing;
project_name = "unnamedProject";
output_dir = pwd();
language = 0;
gen_files = nothing;
solver = nothing;
dendro_params = nothing;
#log
use_log = false;
log_file = "";
log_line_index = 1;
#mesh
mesh_data = nothing;
grid_data = nothing;
refel = nothing;
elemental_order = [];
#problem variables
var_count = 0;
variables = [];
coefficients = [];
parameters = [];
test_functions = [];
#generated functions
genfunc_count = 0;
genfunctions = [];
#rhs
linears = [];
#lhs
bilinears = [];
#time stepper
time_stepper = nothing;

use_cachesim = false;

include("femshop_includes.jl");
include("macros.jl"); # included here after globals are defined

config = Femshop_config(); # These need to be initialized here
prob = Femshop_prob();

function init_femshop(name="unnamedProject")
    global config = Femshop_config();
    global prob = Femshop_prob();
    global project_name = name;
    global language = JULIA;
    global gen_files = nothing;
    global dendro_params = nothing;
    global log_file = "";
    global use_log = false;
    global log_line_index = 1;
    global mesh_data = nothing;
    global grid_data = nothing;
    global refel = nothing;
    global elemental_order = [];
    global var_count = 0;
    global variables = [];
    global coefficients = [];
    global parameters = [];
    global test_functions = [];
    global genfunc_count = 0;
    global genfunctions = [];
    global linears = [];
    global bilinears = [];
    global time_stepper = nothing;
    global use_cachesim = false;
end

function set_language(lang, dirpath, name, head="")
    global language = lang;
    global output_dir = dirpath;
    global project_name = name;
    global gen_files = CodeGenerator.init_codegenerator(lang, dirpath, name, head);
end

function dendro(;max_depth=6, wavelet_tol=0.1, partition_tol=0.3, solve_tol=1e-6, max_iters=100)
    global dendro_params = (max_depth, wavelet_tol, partition_tol, solve_tol, max_iters);
end

function set_solver(s)
    global solver = s;
end

function set_stepper(type, cfl)
    global time_stepper = Stepper(type, cfl);
    global config.stepper = type;
    log_entry("Set time stepper to "*type);
end

function set_matrix_free(max, tol)
    config.linalg_matrixfree = true;
    config.linalg_matfree_max = max;
    config.linalg_matfree_tol = tol;
end

function add_mesh(mesh)
    if typeof(mesh) <: Tuple
        global mesh_data = mesh[1];
        global refel = mesh[2];
        global grid_data = mesh[3];
    else
        global mesh_data = mesh;
    end
    # set elemental loop ordering
    global elemental_order = 1:mesh_data.nel;

    log_entry("Added mesh with "*string(mesh_data.nx)*" vertices and "*string(mesh_data.nel)*" elements.");
    log_entry("Full grid has "*string(length(grid_data.allnodes))*" nodes.");
end

function output_mesh(file, format)
    write_mesh(file, format, mesh_data);
    log_entry("Wrote mesh data to file: "*file*".msh");
end

function add_test_function(v, type)
    varind = length(test_functions) + 1;
    # make SymType
    symvar = sym_var(string(v), type, config.dimension);

    push!(test_functions, Femshop.Coefficient(v, symvar, varind, type, []););
    log_entry("Set test function symbol: "*string(v)*" of type: "*type);
end

function add_variable(var)
    global var_count += 1;
    if language == JULIA || language == 0
        # adjust values arrays
        N = size(grid_data.allnodes)[1];
        if var.type == SCALAR
            var.values = zeros(N);
        elseif var.type == VECTOR
            var.values = zeros(N, config.dimension);
        elseif var.type == TENSOR
            var.values = zeros(N, config.dimension*config.dimension);
        elseif var.type == SYM_TENSOR
            var.values = zeros(N, Int((config.dimension*(config.dimension+1))/2));
        end
    end
    # make SymType
    symvar = sym_var(string(var.symbol), var.type, config.dimension);
    var.symvar = symvar;

    global variables = [variables; var];

    global linears = [linears; nothing];
    global bilinears = [bilinears; nothing];

    log_entry("Added variable: "*string(var.symbol)*" of type: "*var.type);
end

function add_coefficient(c, type, val, nfuns)
    global coefficients;
    vals = [];
    if nfuns == 0 # constant values
        vals = val;
        if length(vals) == 1 && !(typeof(vals) <: Array)
            vals = [val];
        end

    else # genfunction values
        if typeof(val) <: Array
            ind = length(genfunctions) - nfuns + 1;
            for i=1:length(val)
                if typeof(val[i]) == String
                    push!(vals, genfunctions[ind]);
                    ind += 1;
                else
                    push!(vals, val[i]);
                end
            end
        else
            push!(vals, genfunctions[end]);
        end
    end

    symvar = sym_var(string(c), type, config.dimension);

    index = length(coefficients);
    push!(coefficients, Coefficient(c, symvar, index, type, vals));

    log_entry("Added coefficient "*string(c)*" : "*string(val));

    return coefficients[end];
end

function add_parameter(p, type, val)
    
    index = length(parameters);
    push!(parameters, Parameter(p, index, type, val));

    log_entry("Added parameter "*string(p)*" : "*string(val));

    return parameters[end];
end

function swap_parameter_xyzt(ex)
    if typeof(ex) == Symbol
        if ex === :x
            return :parameterCoefficientForx ;
        elseif ex === :y
            return :parameterCoefficientFory ;
        elseif ex === :z
            return :parameterCoefficientForz ;
        elseif ex === :t
            return :parameterCoefficientFort ;
        end
    elseif typeof(ex) == Expr
        for i=1:length(ex.args)
            ex.args[i] = swap_parameter_xyzt(ex.args[i]); # Recursively swap
        end
    end
    return ex;
end

function add_initial_condition(varindex, ex, nfuns)
    global prob;
    while length(prob.initial) < varindex
        prob.initial = [prob.initial; nothing];
    end
    if typeof(ex) <: Array
        vals = [];
        ind = length(genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) == String
                push!(vals, genfunctions[ind]);
                ind += 1;
            else
                push!(vals, ex[i]);
            end
        end
        prob.initial[varindex] = vals;
    else
        if typeof(ex) == String
            prob.initial[varindex] = genfunctions[end];
        else
            prob.initial[varindex] = ex;
        end
    end

    log_entry("Initial condition for "*string(variables[varindex].symbol)*" : "*string(prob.initial[varindex]));
    # hold off on initializing till solve or generate is determined.
end

function add_boundary_condition(var, bid, type, ex, nfuns)
    global prob;
    # make sure the arrays are big enough
    if size(prob.bc_func)[1] < var_count || size(prob.bc_func)[2] < bid
        tmp1 = Array{String,2}(undef, (var_count, bid));
        tmp2 = Array{Any,2}(undef, (var_count, bid));
        tmp3 = zeros(Int, (var_count, bid));
        fill!(tmp1, "");
        fill!(tmp2, GenFunction("","","",0,0));
        tmp1[1:size(prob.bc_func)[1], 1:size(prob.bc_func)[2]] = Base.deepcopy(prob.bc_type);
        tmp2[1:size(prob.bc_func)[1], 1:size(prob.bc_func)[2]] = Base.deepcopy(prob.bc_func);
        tmp3[1:size(prob.bc_func)[1], 1:size(prob.bc_func)[2]] = prob.bid;
        prob.bc_type = tmp1;
        prob.bc_func = tmp2;
        prob.bid = tmp3;
    end
    if typeof(ex) <: Array
        vals = [];
        ind = length(genfunctions) - nfuns + 1;
        for i=1:length(ex)
            if typeof(ex[i]) == String
                push!(vals, genfunctions[ind]);
                ind += 1;
            else
                push!(vals, ex[i]);
            end
        end
        prob.bc_func[var.index, bid] = vals;
    else
        if typeof(ex) == String
            prob.bc_func[var.index, bid] = [genfunctions[end]];
        else
            prob.bc_func[var.index, bid] = [ex];
        end
    end
    prob.bc_type[var.index, bid] = type;
    prob.bid[var.index, bid] = bid;

    log_entry("Boundary condition: var="*string(var.symbol)*" bid="*string(bid)*" type="*type*" val="*string(ex));
end

function set_rhs(var, code="")
    global linears;
    if language == 0 || language == JULIA
        if typeof(var) <:Array
            for i=1:length(var)
                linears[var[i].index] = genfunctions[end];
            end
        else
            linears[var.index] = genfunctions[end];
        end
        
    else # external generation
        if typeof(var) <:Array
            for i=1:length(var)
                linears[var[i].index] = code;
            end
        else
            linears[var.index] = code;
        end
    end
end

function set_lhs(var, code="")
    global bilinears;
    if language == 0 || language == JULIA
        if typeof(var) <:Array
            for i=1:length(var)
                bilinears[var[i].index] = genfunctions[end];
            end
        else
            bilinears[var.index] = genfunctions[end];
        end
        
    else # external generation
        if typeof(var) <:Array
            for i=1:length(var)
                bilinears[var[i].index] = code;
            end
        else
            bilinears[var.index] = code;
        end
    end
end

function cachesim_solve(var, nlvar=nothing; nonlinear=false)
    if !(gen_files === nothing && (language == JULIA || language == 0))
        printerr("Cachesim solve is only ready for Julia direct solve");
    else
        if typeof(var) <: Array
            varnames = "["*string(var[1].symbol);
            for vi=2:length(var)
                varnames = varnames*", "*string(var[vi].symbol);
            end
            varnames = varnames*"]";
            varind = var[1].index;
        else
            varnames = string(var.symbol);
            varind = var.index;
        end

        lhs = bilinears[varind];
        rhs = linears[varind];
        
        t = @elapsed(result = CGSolver.linear_solve_cachesim(var, lhs, rhs));
        log_entry("Generated cachesim ouput for "*varnames*".(took "*string(t)*" seconds)");
    end
end

function solve(var, nlvar=nothing; nonlinear=false)
    if use_cachesim
        printerr("Use cachesim_solve(var) for generating cachesim output. Try again.");
        return nothing;
    end
	#@show(var[1].index)
    # Generate files or solve directly
    if prob.time_dependent
        global time_stepper = init_stepper(grid_data.allnodes, time_stepper);
    end
    if !(gen_files === nothing && (language == JULIA || language == 0))
        generate_main();
        if !(dendro_params === nothing)
            generate_config(dendro_params);
        else
            generate_config();
        end
        generate_prob();
        generate_mesh();
        generate_genfunction();
        generate_bilinear(bilinears[1]);
        generate_linear(linears[1]);
        #generate_stepper();
        generate_output();
    else
        if config.solver_type == CG
            init_cgsolver();
            if typeof(var) <: Array
                varnames = "["*string(var[1].symbol);
                for vi=2:length(var)
                    varnames = varnames*", "*string(var[vi].symbol);
                end
                varnames = varnames*"]";
                varind = var[1].index;
            else
                varnames = string(var.symbol);
                varind = var.index;
            end

            lhs = bilinears[varind];
            rhs = linears[varind];
            
            if prob.time_dependent
                global time_stepper = init_stepper(grid_data.allnodes, time_stepper);
				if (nonlinear)
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs, time_stepper));
				else
                	t = @elapsed(result = CGSolver.linear_solve(var, lhs, rhs, time_stepper));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	t = @elapsed(result = CGSolver.nonlinear_solve(var, nlvar, lhs, rhs));
                else
                    t = @elapsed(result = CGSolver.linear_solve(var, lhs, rhs));
				end

                # place the values in the variable value arrays
                if typeof(var) <: Array && length(result) > 1
                    tmp = 0;
                    totalcomponents = 0;
                    for vi=1:length(var)
                        totalcomponents = totalcomponents + length(var[vi].symvar.vals);
                    end
                    for vi=1:length(var)
                        components = length(var[vi].symvar.vals);
                        for compi=1:components
                            var[vi].values[:,compi] = result[(compi+tmp):totalcomponents:end];
                            tmp = tmp + 1;
                        end
                    end
                elseif length(result) > 1
                    components = length(var.symvar.vals);
                    for compi=1:components
                        var.values[:,compi] = result[compi:components:end];
                    end
                end
            end

            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)");
        end
    end

end

function finalize()
    # Finalize generation
    if gen_files != nothing
        finalize_codegenerator();
    end
    if use_cachesim
        CachesimOut.finalize();
    end
    # anything else
    close_log();
    println("Femshop has completed.");
end

function cachesim(use)
    log_entry("Using cachesim - Only cachesim output will be generated.");
    global use_cachesim = use;
end

function morton_nodes(n)
    t = @elapsed(global grid_data = reorder_grid_morton(grid_data, n));
    log_entry("Reordered nodes to Morton. Took "*string(t)*" sec.");
end

function tiled_nodes(griddim, tiledim)
    t = @elapsed(global grid_data = reorder_grid_tiled(grid_data, griddim, tiledim));
    log_entry("Reordered nodes to tiled. Took "*string(t)*" sec.");
end

function tiled_elements(griddim, tiledim)
    global elemental_order = get_tiled_order(config.dimension, griddim, tiledim, false);
    log_entry("Reordered elements to tiled("*string(tiledim)*").");
end

end # module
