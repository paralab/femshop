#=
Module for code generation
=#
module CodeGenerator

export init_code_generator, finalize_code_generator, set_generation_target,
        generate_all_files, add_generated_file,
        # generate_main, generate_config, generate_prob, generate_mesh, generate_genfunction, 
        # generate_bilinear, generate_linear, generate_stepper, generate_output,
        generate_code_layer, generate_code_layer_surface, generate_code_layer_fv

import ..Femshop: JULIA, CPP, MATLAB, DENDRO, HOMG, CUSTOM_GEN_TARGET,
            SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, 
            CG, DG, HDG, FV,
            NODAL, MODAL, CELL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, 
            NONLINEAR_NEWTON, NONLINEAR_SOMETHING, 
            EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, ABM4, 
            DEFAULT_SOLVER, PETSC, 
            VTK, RAW_OUTPUT, CUSTOM_OUTPUT, 
            DIRICHLET, NEUMANN, ROBIN, NO_BC, FLUX,
            MSH_V2, MSH_V4,
            SCALAR, VECTOR, TENSOR, SYM_TENSOR,
            LHS, RHS,
            LINEMESH, QUADMESH, HEXMESH
import ..Femshop: Femshop_config, Femshop_prob, GenFunction, Variable, Coefficient
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, refel, mesh_data, grid_data, genfunctions, variables, coefficients, 
        test_functions, linears, bilinears, time_stepper, language, gen_framework
import ..Femshop: SymExpression, SymEntity
import ..Femshop: CachesimOut, use_cachesim
import ..Femshop: custom_gen_funcs

genDir = "";
genFileName = "";
genFileExtension = "";
commentChar = "";
blockCommentChar = [""; ""];
headerText = "";
genfiles = [];
external_get_language_elements_function = nothing;
external_generate_code_layer_function = nothing;
external_generate_code_files_function = nothing;

# for custom targets
using_custom_target = false;
# Temporary placeholders for external code gen functions that must be provided.
# These are reassigned in set_custom_target()
function default_language_elements_function() return (".jl", "#", ["#=", "=#"]) end;
function default_code_layer_function(var, entities, terms, lorr, vors) return ("","") end;
function default_code_files_function(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) return 0 end;

# general code generator functions
include("code_generator_utils.jl");
include("generate_code_layer.jl");

# code gen functions for each solver type and target
include("generate_code_layer_cg_julia.jl");
include("generate_code_layer_dg_julia.jl");
include("generate_code_layer_fv_julia.jl");

# # target specific code gen functions
# include("generate_code_layer_dendro.jl");
# include("generate_code_layer_homg.jl");
# include("generate_code_layer_matlab.jl");
# include("generate_code_layer_cachesim.jl");

# Surface integrals should be handled in the same place TODO
#include("generate_code_layer_surface.jl");

#Matlab
include("generate_matlab_utils.jl");
include("generate_matlab_files.jl");
include("generate_homg_files.jl");
#C++
include("generate_cpp_utils.jl");
include("generate_dendro_files.jl");

function init_code_generator(dir, name, header)
    # if lang == JULIA
    #     global genFileExtension = ".jl";
    #     global commentChar = "#";
    #     global blockCommentChar = ["#="; "=#"];
    #     log_entry("Set code generation language to Julia.");
    # elseif lang == CPP
    #     global genFileExtension = ".cpp";
    #     global commentChar = "//";
    #     global blockCommentChar = ["/*"; "*/"];
    #     log_entry("Set code generation language to C++.");
    # elseif lang == MATLAB
    #     global genFileExtension = ".m";
    #     global commentChar = "%";
    #     global blockCommentChar = ["%{"; "%}"];
    #     log_entry("Set code generation language to Matlab.");
    # elseif using_custom_target
    #     (ext, com, blo) = external_get_language_elements_function();
    #     global genFileExtension = ext;
    #     global commentChar = com;
    #     global blockCommentChar = blo;
    #     log_entry("Set code generation language to custom.");
    # else
    #     printerr("Invalid language, use JULIA, CPP, MATLAB");
    #     return nothing;
    # end
    
    global genFileExtension = ".jl";
    global commentChar = "#";
    global blockCommentChar = ["#="; "=#"];
    global genDir = dir;
    global genFileName = name;
    global headerText = header;
    
    global external_get_language_elements_function = default_language_elements_function;
    global external_generate_code_layer_function = default_code_layer_function;
    global external_generate_code_files_function = default_code_files_function;

    # src_dir = genDir*"/src";
    # if !isdir(src_dir)
    #     mkdir(src_dir);
    # end
    # # inc_dir = genDir*"/include";
    # # if !isdir(inc_dir)
    # #     mkdir(inc_dir);
    # # end
    
    # m = open(src_dir*"/"*name*genFileExtension, "w");
    # c = open(src_dir*"/Config"*genFileExtension, "w");
    # p = open(src_dir*"/Problem"*genFileExtension, "w");
    # n = open(src_dir*"/Mesh"*genFileExtension, "w");
    # d = open(src_dir*"/MeshData", "w");
    # g = open(src_dir*"/Genfunction"*genFileExtension, "w");
    # b = open(src_dir*"/Bilinear"*genFileExtension, "w");
    # l = open(src_dir*"/Linear"*genFileExtension, "w");
    # s = open(src_dir*"/Stepper"*genFileExtension, "w");
    # o = open(src_dir*"/Output"*genFileExtension, "w");
    
    # global genfiles = Genfiles(m,c,p,n,d,g,b,l,s,o);
    
    # # write headers
    # generate_head(m,headerText);
    # generate_head(c,"Configuration info");
    # generate_head(p,"Problem info");
    # generate_head(n,"Mesh");
    # generate_head(g,"Generated functions");
    # generate_head(b,"Bilinear term");
    # generate_head(l,"Linear term");
    # generate_head(s,"Time stepper");
    # generate_head(o,"Output");
    
    # log_entry("Created code files for: "*name);
end

# Sets the functions to be used during external code generation
function set_generation_target(lang_elements, code_layer, file_maker)
    global external_get_language_elements_function = lang_elements;
    global external_generate_code_layer_function = code_layer;
    global external_generate_code_files_function = file_maker;
    global using_custom_target = true;
end

function add_generated_file(filename; dir="")
    if length(dir) > 0
        code_dir = genDir*"/"*dir;
        if !isdir(code_dir)
            mkdir(code_dir);
        end
    else
        code_dir = genDir;
    end
    newfile = open(code_dir*"/"*filename, "w");
    push!(genfiles, newfile);
    return newfile;
end

# function generate_all_files(lhs_vol, lhs_surf, rhs_vol, rhs_surf)
#     # This is supplied as one of the target functions
#     external_generate_code_files_function(lhs_vol, lhs_surf, rhs_vol, rhs_surf);
# end

function finalize_code_generator()
    for f in genfiles
        close(f);
    end
    log_entry("Closed generated code files.");
end

############# Don't delete yet ####################
# # This is the primary code gen function that calls functions in the appropriate files
# function generate_code_layer(ex, var, lorr, vors)
#     if config.solver_type == FV
#         return generate_code_layer_fv(ex, var, lorr, vors, language, framework);
#     elseif config.solver_type == CG
#         return generate_code_layer_cg(ex, var, lorr, vors, language, framework);
#     elseif config.solver_type == DG
#         return generate_code_layer_dg(ex, var, lorr, vors, language, framework);
#     end
    
#     if use_cachesim
#         if language == 0 || language == JULIA
#             return generate_code_layer_cachesim(ex, var, lorr);
#         else
#             printerr("Using cachesim. Not generating code to solve.")
#         end
#     else
#         if language == 0 || language == JULIA
#             return generate_code_layer_julia(ex, var, lorr, vors);
#         elseif language == CPP
#             if framework == DENDRO
#                 return generate_code_layer_dendro(ex, var, lorr);
#             else
#                 printerr("Plain C++ is not ready for code layer gen.")
#             end
            
#         elseif language == MATLAB
#             if framework == HOMG
#                 return generate_code_layer_homg(ex, var, lorr);
#             else
#                 return generate_code_layer_matlab(ex, var, lorr);
#             end
            
#         elseif framework == CUSTOM_GEN_TARGET
#             return custom_code_layer_fun(ex, var, lorr);
#         end
#     end
# end

function generate_all_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf; parameters=0)
    if language == CPP
        if gen_framework == DENDRO
            if parameters == 0
                parameters = (5, 1, 0.3, 0.000001, 100);#(maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
            end
            dendro_main_file();
            dendro_config_file(parameters);
            dendro_prob_file();
            # dendro_mesh_file();
            dendro_genfunction_file();
            dendro_bilinear_file(lhs_vol);
            dendro_linear_file(rhs_vol);
            # dendro_stepper_file();
            dendro_output_file();
        else
            printerr("Plain C++ not ready")
        end
        
    elseif language == MATLAB
        if gen_framework == HOMG
            homg_main_file();
            homg_config_file();
            homg_prob_file();
            homg_mesh_file();
            homg_genfunction_file();
            homg_bilinear_file(lhs_vol);
            homg_linear_file(rhs_vol);
            homg_stepper_file();
            # homg_output_file();
        else
            matlab_main_file();
            matlab_config_file();
            matlab_prob_file();
            matlab_mesh_file();
            matlab_genfunction_file();
            matlab_bilinear_file(lhs_vol);
            matlab_linear_file(rhs_vol);
            matlab_stepper_file();
            matlab_output_file();
        end
        
    elseif using_custom_target
        external_generate_code_files_function(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf);
        
    end
end

function comment(file,line)
    println(file, commentChar * line);
end

function commentBlock(file,text)
    print(file, "\n"*blockCommentChar[1]*"\n"*text*"\n"*blockCommentChar[2]*"\n");
end

function generate_head(file, text)
    comment(file,"This file was generated by Femshop.");
    commentBlock(file, text);
end

# for writing structs to binary files
# format is | number of structs[Int64] | sizes of structs[Int64*num] | structs |
function write_binary_head(f, num, szs)
    Nbytes = 0;
    write(f, num);
    Nbytes += sizeof(num)
    for i=1:length(szs)
        write(f, szs[i])
        Nbytes += sizeof(szs[i])
    end
    return Nbytes;
end

# Write an array to a binary file.
# Return number of bytes written.
function write_binary_array(f, a)
    Nbytes = 0;
    for i=1:length(a)
        if isbits(a[i])
            write(f, a[i]);
            Nbytes += sizeof(a[i]);
        else
            Nbytes += write_binary_array(f,a[i]);
        end
    end
    return Nbytes;
end

# Assumes that the struct only has isbits->true types or arrays.
# Returns number of bytes written.
function write_binary_struct(f, s)
    Nbytes = 0;
    for fn in fieldnames(typeof(s))
        comp = getfield(s, fn);
        if isbits(comp)
            write(f, comp);
            Nbytes += sizeof(comp)
        else
            Nbytes += write_binary_array(f,comp);
        end
    end
    return Nbytes;
end

end # module