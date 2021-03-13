#=
Module for code generation
=#
module CodeGenerator

export init_codegenerator, finalize_codegenerator, Genfiles, set_custom_target,
        generate_all_files,
        generate_main, generate_config, generate_prob, generate_mesh, generate_genfunction, 
        generate_bilinear, generate_linear, generate_stepper, generate_output,
        generate_code_layer, generate_code_layer_surface

import ..Femshop: JULIA, CPP, MATLAB, DENDRO, HOMG, CUSTOM_GEN_TARGET,
            SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, 
            CG, DG, HDG,
            NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, 
            NONLINEAR_NEWTON, NONLINEAR_SOMETHING, 
            EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, ABM4, 
            DEFAULT_SOLVER, PETSC, 
            VTK, RAW_OUTPUT, CUSTOM_OUTPUT, 
            DIRICHLET, NEUMANN, ROBIN, NO_BC,
            MSH_V2, MSH_V4,
            SCALAR, VECTOR, TENSOR, SYM_TENSOR,
            LHS, RHS,
            LINEMESH, QUADMESH, HEXMESH
import ..Femshop: Femshop_config, Femshop_prob, GenFunction, Variable, Coefficient
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, refel, mesh_data, grid_data, genfunctions, variables, coefficients, 
        test_functions, linears, bilinears, time_stepper
import ..Femshop: CachesimOut, use_cachesim
import ..Femshop: custom_gen_funcs

# Holds a set of file streams for generated code
mutable struct Genfiles
    main;       # Runs the computation
    config;     # global configuration
    problem;    # global problem specification
    mesh;       # contains the mesh CODE
    meshdata;   # contains the refel/mesh/grid DATA
    genfunction;# Generated functions
    bilinear;   # bilinear function: bilinear(args) returns elemental matrix
    linear;     # linear function: linear(args) returns elemental vector
    stepper;    # optional time stepper for time dependent problems
    output;     # output
    
    files;      # an iterable list of these files
    
    Genfiles(m,c,p,n,nd,g,b,l,s,o) = new(m,c,p,n,nd,g,b,l,s,o,[m,c,p,n,nd,g,b,l,s,o]);
end

language = 0;
framework = 0;
genDir = "";
genFileName = "";
genFileExtension = "";
commentChar = "";
blockCommentChar = [""; ""];
headerText = "";
genfiles = nothing;

# for custom targets
using_custom_target = false;
custom_language_element_fun = 0;
custom_code_layer_fun = 0;
custom_files_fun = 0;

#General
include("generate_code_layer.jl");
#Matlab
include("generate_matlab_utils.jl");
include("generate_matlab_files.jl");
include("generate_homg_files.jl");
#C++
include("generate_cpp_utils.jl");
include("generate_dendro_files.jl");

function set_custom_target(lang_elements, code_layer, file_maker)
    global custom_language_element_fun = lang_elements;
    global custom_code_layer_fun = code_layer;
    global custom_files_fun = file_maker;
    global using_custom_target = true;
end

function init_codegenerator(lang, frame, dir, name, header)
    if lang == JULIA
        global genFileExtension = ".jl";
        global commentChar = "#";
        global blockCommentChar = ["#="; "=#"];
        log_entry("Set code generation language to Julia.");
    elseif lang == CPP
        global genFileExtension = ".cpp";
        global commentChar = "//";
        global blockCommentChar = ["/*"; "*/"];
        log_entry("Set code generation language to C++.");
    elseif lang == MATLAB
        global genFileExtension = ".m";
        global commentChar = "%";
        global blockCommentChar = ["%{"; "%}"];
        log_entry("Set code generation language to Matlab.");
    elseif using_custom_target
        (ext, com, blo) = custom_language_element_fun();
        global genFileExtension = ext;
        global commentChar = com;
        global blockCommentChar = blo;
        log_entry("Set code generation language to custom.");
    else
        printerr("Invalid language, use JULIA, CPP, MATLAB");
        return nothing;
    end
    
    global language = lang;
    global framework = frame;
    global genDir = dir;
    global genFileName = name;
    global headerText = header;
    
    src_dir = genDir*"/src";
    if !isdir(src_dir)
        mkdir(src_dir);
    end
    # inc_dir = genDir*"/include";
    # if !isdir(inc_dir)
    #     mkdir(inc_dir);
    # end
    
    m = open(src_dir*"/"*name*genFileExtension, "w");
    c = open(src_dir*"/Config"*genFileExtension, "w");
    p = open(src_dir*"/Problem"*genFileExtension, "w");
    n = open(src_dir*"/Mesh"*genFileExtension, "w");
    d = open(src_dir*"/MeshData", "w");
    g = open(src_dir*"/Genfunction"*genFileExtension, "w");
    b = open(src_dir*"/Bilinear"*genFileExtension, "w");
    l = open(src_dir*"/Linear"*genFileExtension, "w");
    s = open(src_dir*"/Stepper"*genFileExtension, "w");
    o = open(src_dir*"/Output"*genFileExtension, "w");
    
    global genfiles = Genfiles(m,c,p,n,d,g,b,l,s,o);
    
    # write headers
    generate_head(m,headerText);
    generate_head(c,"Configuration info");
    generate_head(p,"Problem info");
    generate_head(n,"Mesh");
    generate_head(g,"Generated functions");
    generate_head(b,"Bilinear term");
    generate_head(l,"Linear term");
    generate_head(s,"Time stepper");
    generate_head(o,"Output");
    
    log_entry("Created code files for: "*name);
    
    return genfiles;
end

function finalize_codegenerator()
    for f in genfiles.files
        close(f);
    end
    log_entry("Closed generated code files.");
end

function generate_all_files(bilinex, linex; parameters=0)
    if language == CPP
        if framework == DENDRO
            if parameters == 0
                parameters = (5, 1, 0.3, 0.000001, 100);#(maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
            end
            dendro_main_file();
            dendro_config_file(parameters);
            dendro_prob_file();
            # dendro_mesh_file();
            dendro_genfunction_file();
            dendro_bilinear_file(bilinex);
            dendro_linear_file(linex);
            # dendro_stepper_file();
            dendro_output_file();
        else
            printerr("Plain C++ not ready")
        end
        
    elseif language == MATLAB
        if framework == HOMG
            homg_main_file();
            homg_config_file();
            homg_prob_file();
            homg_mesh_file();
            homg_genfunction_file();
            homg_bilinear_file(bilinex);
            homg_linear_file(linex);
            homg_stepper_file();
            # homg_output_file();
        else
            matlab_main_file();
            matlab_config_file();
            matlab_prob_file();
            matlab_mesh_file();
            matlab_genfunction_file();
            matlab_bilinear_file(bilinex);
            matlab_linear_file(linex);
            matlab_stepper_file();
            # matlab_output_file();
        end
        
    elseif using_custom_target
        custom_files_fun(genDir, genfiles, bilinex, linex);
        
    end
end

# # Select the generator functions based on language
# function generate_main()
#     if language == CPP
#         if framework == DENDRO
#             dendro_main_file();
#         else
#             printerr("Plain C++ not ready")
#         end
        
#     elseif language == MATLAB
#         if framework == HOMG
#             homg_main_file();
#         else
#             matlab_main_file();
#         end
        
#     end
# end
# function generate_config(params=(5, 1, 0.3, 0.000001, 100))
#     if language == CPP
#         dendro_config_file(params); #params has (maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
#     elseif language == DENDRO
#         dendro_config_file(params); #params has (maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
#     elseif language == MATLAB
#         matlab_config_file();
#     elseif language == HOMG
#         homg_config_file();
#     end
# end
# function generate_prob()
#     if language == CPP
#         dendro_prob_file();
#     elseif language == DENDRO
#         dendro_prob_file();
#     elseif language == MATLAB
#         matlab_prob_file();
#     elseif language == HOMG
#         homg_prob_file();
#     end
# end
# function generate_mesh()
#     if language == CPP
#         # dendro_mesh_file();
#     elseif language == DENDRO
#         # dendro_mesh_file();
#     elseif language == MATLAB
#         matlab_mesh_file();
#     elseif language == HOMG
#         homg_mesh_file();
#     end
# end
# function generate_genfunction()
#     if language == CPP
#         dendro_genfunction_file();
#     elseif language == DENDRO
#         dendro_genfunction_file();
#     elseif language == MATLAB
#         matlab_genfunction_file();
#     elseif language == HOMG
#         homg_genfunction_file();
#     end
# end
# function generate_bilinear(ex)
#     if language == CPP
#         dendro_bilinear_file(ex);
#     elseif language == DENDRO
#         dendro_bilinear_file(ex);
#     elseif language == MATLAB
#         matlab_bilinear_file(ex);
#     elseif language == HOMG
#         homg_bilinear_file(ex);
#     end
# end
# function generate_linear(ex)
#     if language == CPP
#         dendro_linear_file(ex);
#     elseif language == DENDRO
#         dendro_linear_file(ex);
#     elseif language == MATLAB
#         matlab_linear_file(ex);
#     elseif language == HOMG
#         homg_linear_file(ex);
#     end
# end
# function generate_stepper()
#     if language == CPP
#         # dendro_stepper_file();
#     elseif language == DENDRO
#         # dendro_stepper_file();
#     elseif language == MATLAB
#         matlab_stepper_file();
#     elseif language == HOMG
#         homg_stepper_file();
#     end
# end
# function generate_output()
#     if language == CPP
#         dendro_output_file();
#     elseif language == DENDRO
#         dendro_output_file();
#     elseif language == MATLAB
#         # matlab_output_file();
#     elseif language == HOMG
#         # homg_output_file();
#     end
# end
        
# private ###

macro comment(file,line)
    return esc(quote
        println($file, commentChar * $line);
    end)
end

macro commentBlock(file,text)
    return esc(quote
        print($file, "\n"*blockCommentChar[1]*"\n"*text*"\n"*blockCommentChar[2]*"\n");
    end)
end

function generate_head(file, text)
    @comment(file,"This file was generated by Femshop.");
    @commentBlock(file, text);
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
