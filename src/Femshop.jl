#=
The main module for Femshop.
=#
module Femshop

# Public macros and functions
export @generateFor, @domain, @mesh, @solver, @stepper, @setSteps, @functionSpace, @trialFunction, @matrixFree,
        @testFunction, @nodes, @order, @boundary, @addBoundaryID, @referencePoint, @variable, @coefficient, @parameter, @testSymbol, @initial,
        @timeInterval, @weakForm, @fluxAndSource, @LHS, @RHS, @exportCode, @importCode,
        @customOperator, @customOperatorFile,
        @outputMesh, @useLog, @finalize
export init_femshop, set_language, set_custom_gen_target, dendro, set_solver, set_stepper, set_specified_steps, set_matrix_free,
        reformat_for_stepper, reformat_for_stepper_fv,
        add_mesh, output_mesh, add_boundary_ID, add_test_function, 
        add_initial_condition, add_boundary_condition, add_reference_point, set_rhs, set_lhs, set_lhs_surface, set_rhs_surface,
        set_symexpressions
export build_cache_level, build_cache, build_cache_auto
export sp_parse
export generate_code_layer, generate_code_layer_surface, generate_code_layer_fv
export Variable, add_variable
export Coefficient, add_coefficient
export Parameter, add_parameter
export FVInfo

### Module's global variables ###
# configuration, output files, etc.
config = nothing;
prob = nothing;
project_name = "unnamedProject";
output_dir = pwd();
language = 0;
gen_framework = 0;
generate_external = false;
solver = nothing;
dendro_params = nothing;
#log
use_log = false;
log_file = "";
log_line_index = 1;
#mesh
mesh_data = nothing;    # The basic element information as read from a MSH file or generated here.
grid_data = nothing;    # The full collection of nodes(including internal nodes) and other mesh info in the actual DOF ordering.
geo_factors = nothing;  # Geometric factors
fv_info = nothing;      # Finite volume info
refel = nothing;        # Reference element (will eventually be an array of refels)
elemental_order = [];   # Determines the order of the elemental loops
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
face_linears = [];
#lhs
bilinears = [];
face_bilinears = [];
#symbolic layer
symexpressions = [[],[],[],[]];
#time stepper
time_stepper = nothing;
specified_dt = 0;
specified_Nsteps = 0;
use_specified_steps = false;

use_cachesim = false;

#handles for custom code gen functions
custom_gen_funcs = [];

###############################################################
include("femshop_includes.jl");
include("macros.jl");
include("femshop_interface.jl"); # included here after globals are defined
###############################################################

config = Femshop_config(); # These need to be initialized here
prob = Femshop_prob();

# This has already been initialized on loading, so this is more of a reset function
function init_femshop(name="unnamedProject")
    global config = Femshop_config();
    global prob = Femshop_prob();
    global project_name = name;
    global language = JULIA;
    global framework = 0;
    global generate_external = false;
    global dendro_params = nothing;
    global log_file = "";
    global use_log = false;
    global log_line_index = 1;
    global mesh_data = nothing;
    global grid_data = nothing;
    global geo_factors = nothing;
    global fv_info = nothing;
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
    global face_linears = [];
    global face_bilinears = [];
    global symexpressions = [[],[],[],[]];
    global time_stepper = nothing;
    global specified_dt = 0;
    global specified_Nsteps = 0;
    global use_specified_steps = false;
    global use_cachesim = false;
end

# Sets the code generation target
function set_included_gen_target(lang, framework, dirpath, name; head="")
    #TODO
    println("Premade generation targets currently under construction. Sorry. Reverting to Julia target.");
    return;
    
    global language = lang;
    global gen_framework = framework;
    global output_dir = dirpath;
    global project_name = name;
    global generate_external = true;
    init_codegenerator(dirpath, name, head);
    
    # Need to set these three functions
    
    set_generation_target(get_external_language_elements, generate_external_code_layer, generate_external_files);
end

# Setting a custom target requires three functions
# 1. get_external_language_elements() - file extensions, comment chars etc.
# 2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
# 3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code
function set_custom_gen_target(lang_elements, code_layer, file_maker, dirpath, name; head="")
    global language = CUSTOM_GEN_TARGET;
    global gen_framework = CUSTOM_GEN_TARGET;
    global output_dir = dirpath;
    global project_name = name;
    global generate_external = true;
    init_code_generator(dirpath, name, head);
    set_generation_target(lang_elements, code_layer, file_maker);
end

# Specific to Dendro. This should go into the Dendro code gen module eventually.
function dendro(;max_depth=6, wavelet_tol=0.1, partition_tol=0.3, solve_tol=1e-6, max_iters=100)
    global dendro_params = (max_depth, wavelet_tol, partition_tol, solve_tol, max_iters);
end

# Set the solver module and type
function set_solver(stype)
    config.solver_type = stype;
    if stype == DG
        global solver = DGSolver;
    elseif stype == CG
        global solver = CGSolver;
    elseif stype == FV
        global solver = FVSolver;
        sourcedir = @__DIR__
        add_custom_op_file(string(sourcedir)*"/fv_ops.jl");
    end
end

# Sets the time stepper
function set_stepper(type, cfl)
    global time_stepper = Stepper(type, cfl);
    global config.stepper = type;
    log_entry("Set time stepper to "*type, 1);
end

# Time steps can be chosen automatically or specified. This is for manual specifying.
function set_specified_steps(dt, steps)
    global specified_dt = dt;
    global specified_Nsteps = steps;
    global use_specified_steps = true;
    prob.time_dependent = true;
    prob.end_time = dt*steps;
    log_entry("Set time stepper values to dt="*string(dt)*", Nsteps="*string(steps), 1);
end

# Uses matrix free iterative solver with the given max iterations and error tolerance.
function set_matrix_free(max, tol)
    config.linalg_matrixfree = true;
    config.linalg_matfree_max = max;
    config.linalg_matfree_tol = tol;
end

# Adds a mesh and builds the full grid and reference elements.
# Actually, the mesh object passed here could be either just a MeshData struct,
# or a tuple of MeshData, Refel and Grid already built.
function add_mesh(mesh)
    if typeof(mesh) <: Tuple
        global mesh_data = mesh[1];
        global refel = mesh[2];
        global grid_data = mesh[3];
        
    else
        global mesh_data = mesh;
        global refel;
        global grid_data;
        (refel, grid_data) = grid_from_mesh(mesh_data);
        
    end
    
    constantJ = true;
    do_faces = false;
    do_vol = false;
    if config.solver_type == DG || config.solver_type == FV
        do_faces = true;
    end
    if config.solver_type == FV
        do_vol = true;
    end
    
    global geo_factors = build_geometric_factors(refel, grid_data, do_face_detj=do_faces, do_vol_area=do_vol, constant_jacobian=constantJ);
    
    log_entry("Added mesh with "*string(mesh_data.nx)*" vertices and "*string(mesh_data.nel)*" elements.", 1);
    log_entry("Full grid has "*string(size(grid_data.allnodes,2))*" nodes.", 2);
    
    # If using FV, set up the FV_info
    if config.solver_type == FV
        global fv_info = build_FV_info(grid_data);
    end
    
    # Set elemental loop ordering to match the order from mesh for now.
    global elemental_order = 1:mesh_data.nel;
end

# Write the mesh to a MSH file
function output_mesh(file, format)
    write_mesh(file, format, mesh_data);
    log_entry("Wrote mesh data to file.", 1);
end

# Maybe remove. This is in grid.jl
function add_boundary_ID(bid, on_bdry)
    add_boundary_ID_to_grid(bid, on_bdry, grid_data);
end

# Defines a test function symbol by creating a special coefficient object.
function add_test_function(v, type)
    varind = length(test_functions) + 1;
    # make SymType
    symvar = sym_var(string(v), type, config.dimension);

    push!(test_functions, Femshop.Coefficient(v, symvar, varind, type, NODAL, []););
    log_entry("Set test function symbol: "*string(v)*" of type: "*type, 2);
end

# Adds a variable and allocates everything associated with it.
function add_variable(var)
    global var_count += 1;
    if language == JULIA || language == 0
        # adjust values arrays
        if var.location == CELL
            N = size(grid_data.loc2glb, 2);
        else
            N = size(grid_data.allnodes,2);
        end
        
        if var.type == SCALAR
            var.values = zeros(1, N);
        elseif var.type == VECTOR
            var.values = zeros(config.dimension, N);
        elseif var.type == TENSOR
            var.values = zeros(config.dimension*config.dimension, N);
        elseif var.type == SYM_TENSOR
            var.values = zeros(Int((config.dimension*(config.dimension+1))/2), N);
        end
    end
    # make SymType
    symvar = sym_var(string(var.symbol), var.type, config.dimension);
    var.symvar = symvar;

    global variables = [variables; var];

    global linears = [linears; nothing];
    global bilinears = [bilinears; nothing];
    global face_linears = [face_linears; nothing];
    global face_bilinears = [face_bilinears; nothing];
    global symexpressions[1] = [symexpressions[1]; nothing];
    global symexpressions[2] = [symexpressions[2]; nothing];
    global symexpressions[3] = [symexpressions[3]; nothing];
    global symexpressions[4] = [symexpressions[4]; nothing];

    log_entry("Added variable: "*string(var.symbol)*" of type: "*var.type*", location: "*var.location, 2);
end

# Adds a coefficient with either constant value or some generated function of (x,y,z,t)
function add_coefficient(c, type, location, val, nfuns)
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
    push!(coefficients, Coefficient(c, symvar, index, type, location, vals));

    log_entry("Added coefficient "*string(c)*" : "*string(val), 2);

    return coefficients[end];
end

# Adds a parameter entity.
function add_parameter(p, type, val)
    index = length(parameters);
    push!(parameters, Parameter(p, index, type, val));

    log_entry("Added parameter "*string(p)*" : "*string(val), 2);

    return parameters[end];
end

# Needed to insert parameter expressions into the weak form expressions.
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

# Sets initial condition for the given variable, but does not initialize.
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

    log_entry("Initial condition for "*string(variables[varindex].symbol)*" : "*string(prob.initial[varindex]), 2);
    # hold off on initializing till solve or generate is determined.
end

# Sets boundary condition for a variable and BID.
function add_boundary_condition(var, bid, type, ex, nfuns)
    global prob;
    # make sure the arrays are big enough
    if size(prob.bc_func)[1] < var_count || size(prob.bc_func,2) < bid
        nbid = length(grid_data.bids);
        tmp1 = Array{String,2}(undef, (var_count, nbid));
        tmp2 = Array{Any,2}(undef, (var_count, nbid));
        tmp3 = zeros(Int, (var_count, nbid));
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

    log_entry("Boundary condition: var="*string(var.symbol)*" bid="*string(bid)*" type="*type*" val="*string(ex), 2);
end

# A reference point is a single point with a defined value.
# It's treated like a Dirichlet boundary with one boundary node.
function add_reference_point(var, pos, val)
    global prob;
    # make sure the array is big enough
    if size(prob.ref_point,1) < var_count
        tmp = Array{Any,2}(undef, (var_count, 3));
        for i=1:size(prob.ref_point,1)
            tmp[i,:] = prob.ref_point[i,:];
        end
        for i=(size(prob.ref_point,1)+1):var_count
            tmp[i,1] = false;
            tmp[i,2] = [0,0];
            tmp[i,3] = [0];
        end
        prob.ref_point = tmp;
    end
    if typeof(pos) <: Number
        pos = pos*ones(config.dimension);
    end
    if typeof(val) <: Number
        val = val*ones(length(var.symvar.vals));
    end
    
    # Find the closest vertex to pos
    # The stored pos is actually the index into glbvertex pointing to the closest vertex
    ind = [1,1];
    mindist = 12345;
    for ei=1:mesh_data.nel
        for i=1:size(grid_data.glbvertex,1)
            d = 0;
            for comp=1:length(pos)
                d = d + abs(grid_data.allnodes[comp, grid_data.glbvertex[i, ei]] - pos[comp]);
            end
            if d<mindist
                ind = [i,ei];
                mindist = d;
            end
        end
    end
    
    prob.ref_point[var.index, 1] = true;
    prob.ref_point[var.index, 2] = ind;
    prob.ref_point[var.index, 3] = val;
    
    log_entry("Reference point: var="*string(var.symbol)*" position="*string(pos)*" value="*string(val), 2);
end

# Sets the RHS(linear) function for the given variable(s).
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

# Sets the LHS(bilinear) function for the given variable(s).
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

# Sets the surface RHS(linear) function for the given variable(s).
function set_rhs_surface(var, code="")
    global face_linears;
    if language == 0 || language == JULIA
        if typeof(var) <:Array
            for i=1:length(var)
                face_linears[var[i].index] = genfunctions[end];
            end
        else
            face_linears[var.index] = genfunctions[end];
        end
        
    else # external generation
        if typeof(var) <:Array
            for i=1:length(var)
                face_linears[var[i].index] = code;
            end
        else
            face_linears[var.index] = code;
        end
    end
end

# Sets the surface LHS(bilinear) function for the given variable(s).
function set_lhs_surface(var, code="")
    global face_bilinears;
    if language == 0 || language == JULIA
        if typeof(var) <:Array
            for i=1:length(var)
                face_bilinears[var[i].index] = genfunctions[end];
            end
        else
            face_bilinears[var.index] = genfunctions[end];
        end
        
    else # external generation
        if typeof(var) <:Array
            for i=1:length(var)
                face_bilinears[var[i].index] = code;
            end
        else
            face_bilinears[var.index] = code;
        end
    end
end

function set_symexpressions(var, ex, lorr, vors)
    # If the ex is empty, don't set anything
    if ex == [] || ex == [[]]
        return;
    end
    
    if typeof(var) <:Array
        for i=1:length(var)
            set_symexpressions(var[i], ex, lorr, vors);
        end
        
    else
        global symexpressions;
        if lorr == LHS
            if vors == "volume"
                ind = 1;
            else
                ind = 2;
            end
        else # rhs
            if vors == "volume"
                ind = 3;
            else
                ind = 4;
            end
        end
        symexpressions[ind][var.index] = ex;
    end
end

function finalize()
    finalize_femshop();
end

end # module
