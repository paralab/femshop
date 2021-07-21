#=
This file contains all of the common interface functions.
Many of them simply call corresponding functions in jl.
=#
export generateFor, useLog, domain, solverType, functionSpace, trialSpace, testSpace, 
        nodeType, timeStepper, setSteps, matrixFree, customOperator, customOperatorFile,
        mesh, exportMesh, variable, coefficient, parameter, testSymbol, boundary, addBoundaryID,
        referencePoint, timeInterval, initial, weakForm, fluxAndSource, flux, source, exportCode, importCode,
        solve, cachesimSolve, finalize_femshop, cachesim, output_values,
        morton_nodes, hilbert_nodes, tiled_nodes, morton_elements, hilbert_elements, 
        tiled_elements, ef_nodes, random_nodes, random_elements

# Begin configuration setting functions

function generateFor(lang; filename=project_name, header="")
    outputDirPath = pwd()*"/"*uppercasefirst(filename);
    if !isdir(outputDirPath)
        mkdir(outputDirPath);
    end
    framew = 0;
    if !in(lang, [CPP,MATLAB,DENDRO,HOMG])
        # lang should be a filename for a custom target
        # This file must include these three functions:
        # 1. get_external_language_elements() - file extensions, comment chars etc.
        # 2. generate_external_code_layer(var, entities, terms, lorr, vors) - Turns symbolic expressions into code
        # 3. generate_external_files(var, lhs_vol, lhs_surf, rhs_vol, rhs_surf) - Writes all files based on generated code
        include(lang);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    else
        if lang == DENDRO
            framew = DENDRO;
            lang = CPP;
        elseif lang == HOMG
            framew = HOMG;
            lang = MATLAB;
        else
            framew = 0;
        end
        set_included_gen_target(lang, framew, outputDirPath, filename, head=header);
    end
end

function useLog(name=project_name; dir=output_dir, level=2)
    init_log(name, dir, level);
end

function domain(dims; shape=SQUARE, grid=UNIFORM_GRID)
    config.dimension = dims;
    config.geometry = shape;
    config.mesh_type = grid;
end

function solverType(type)
    set_solver(type);
end

function functionSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    config.trial_function = space;
    config.test_function = space;
    if orderMax > orderMin && orderMin >= 0
        config.p_adaptive = true;
        config.basis_order_min = orderMin;
        config.basis_order_max = orderMax;
    else
        config.basis_order_min = max(order, orderMin);
        config.basis_order_max = max(order, orderMin);
    end
end

function trialSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

function testSpace(;space=LEGENDRE, order=0, orderMin=0, orderMax=0)
    #TODO
    functionSpace(space=space, order=order, orderMin=orderMin, orderMax=orderMax);
end

function nodeType(type)
    config.elemental_nodes = type;
end

function timeStepper(type; cfl=0)
    set_stepper(type, cfl);
end

function setSteps(dt, steps)
    set_specified_steps(dt, steps);
end

function matrixFree(shallwe=true; maxiters=100, tol=1e-6)
    config.linalg_matrixfree = shallwe;
    config.linalg_matfree_max = maxiters;
    config.linalg_matfree_tol = tol;
end

function customOperator(name, handle)
    s = Symbol(name);
    add_custom_op(s, handle);
end

function customOperatorFile(filename)
    log_entry("Adding custom operators from file: "*string(filename), 2);
    add_custom_op_file(filename);
end

# End configuration functions, begin problem definition functions

function mesh(msh; elsperdim=5, bids=1, interval=[0,1])
    if msh == LINEMESH
        log_entry("Building simple line mesh with nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(add_mesh(simple_line_mesh(elsperdim.+1, bids, interval)));
        log_entry("Grid building took "*string(meshtime)*" seconds");
        
    elseif msh == QUADMESH
        log_entry("Building simple quad mesh with nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(add_mesh(simple_quad_mesh(elsperdim.+1, bids, interval)));
        log_entry("Grid building took "*string(meshtime)*" seconds");
        
    elseif msh == HEXMESH
        log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string(elsperdim));
        meshtime = @elapsed(add_mesh(simple_hex_mesh(elsperdim.+1, bids, interval)));
        log_entry("Grid building took "*string(meshtime)*" seconds");
        
    else # msh should be a mesh file name
        # open the file and read the mesh data
        mfile = open(msh, "r");
        log_entry("Reading mesh file: "*msh);
        meshtime = @elapsed(mshdat=read_mesh(mfile));
        log_entry("Mesh reading took "*string(meshtime)*" seconds");
        add_mesh(mshdat);
        close(mfile);
    end
end

function exportMesh(filename, format=MSH_V2)
    # open the file to write to
    mfile = open(filename, "w");
    log_entry("Writing mesh file: "*filename);
    output_mesh(mfile, format);
    close(mfile);
end

function variable(name, type=SCALAR, location=NODAL)
    varind = var_count + 1;
    varsym = Symbol(name);
    var = Variable(varsym, nothing, varind, type, location, [], [], false);
    add_variable(var);
    return var;
end

function coefficient(name, val, type=SCALAR, location=NODAL)
    csym = Symbol(name);
    nfuns = @makeFunctions(val); # if val is constant, nfuns will be 0
    return add_coefficient(csym, type, location, val, nfuns);
end

function parameter(name, val, type=SCALAR)
    if length(parameters) == 0
        coefficient("parameterCoefficientForx", "x")
        coefficient("parameterCoefficientFory", "y")
        coefficient("parameterCoefficientForz", "z")
        coefficient("parameterCoefficientFort", "t")
    end
    if typeof(val) <: Number
        newval = [val];
    elseif typeof(val) == String
        # newval will be an array of expressions
        # search for x,y,z,t symbols and replace with special coefficients like parameterCoefficientForx
        newval = [swap_parameter_xyzt(Meta.parse(val))];
        
    elseif typeof(val) <: Array
        newval = Array{Expr,1}(undef,length(val));
        for i=1:length(val)
            newval[i] = swap_parameter_xyzt(Meta.parse(val[i]));
        end
        newval = reshape(newval,size(val));
    else
        println("Error: use strings to define parameters");
        newval = 0;
    end
    
    return add_parameter(Symbol(name), type, newval);
end

function testSymbol(symb, type=SCALAR)
    add_test_function(Symbol(symb), type);
end

function boundary(var, bid, bc_type, bc_exp=0)
    nfuns = @makeFunctions(bc_exp);
    add_boundary_condition(var, bid, bc_type, bc_exp, nfuns);
end

function addBoundaryID(bid, trueOnBdry)
    # trueOnBdry(x, y, z) = something # points with x,y,z on this bdry segment evaluate true here
    if typeof(trueOnBdry) == String
        @stringToFunction("trueOnBdry", "x,y=0,z=0", trueOnBdry);
    end
    add_boundary_ID_to_grid(bid, trueOnBdry, grid_data);
end

function referencePoint(var, pos, val)
    add_reference_point(var, pos, val);
end

function timeInterval(T)
    prob.time_dependent = true;
    if time_stepper === nothing
        timeStepper(EULER_IMPLICIT);
    end
    prob.end_time = T;
end

function initial(var, ics)
    nfuns = @makeFunctions(ics);
    add_initial_condition(var.index, ics, nfuns);
end

function weakForm(var, wf)
    if typeof(var) <: Array
        # multiple simultaneous variables
        wfvars = [];
        wfex = [];
        if !(length(var) == length(wf))
            printerr("Error in weak form: # of unknowns must equal # of equations. (example: @weakform([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(wfvars, var[vi].symbol);
            push!(wfex, Meta.parse((wf)[vi]));
        end
    else
        wfex = Meta.parse(wf);
        wfvars = var.symbol;
    end
    
    log_entry("Making weak form for variable(s): "*string(wfvars));
    log_entry("Weak form, input: "*string(wf));
    
    # This is the parsing step. It goes from an Expr to arrays of Basic
    result_exprs = sp_parse(wfex, wfvars);
    if length(result_exprs) == 4 # has surface terms
        (lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr) = result_exprs;
    else
        (lhs_expr, rhs_expr) = result_exprs;
    end
    
    if typeof(lhs_expr) <: Tuple && length(lhs_expr) == 2 # has time derivative
        if length(result_exprs) == 4 # has surface terms
            log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" + surface("*string(lhs_surf_expr)*") = "*string(rhs_expr)*" + surface("*string(rhs_surf_expr)*")");
            
            (newlhs, newrhs, newsurflhs, newsurfrhs) = reformat_for_stepper(lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr, config.stepper);
            
            log_entry("Weak form, modified for time stepping: "*string(newlhs)*" + surface("*string(newsurflhs)*") = "*string(newrhs)*" + surface("*string(newsurfrhs)*")");
            
            lhs_expr = newlhs;
            rhs_expr = newrhs;
            lhs_surf_expr = newsurflhs;
            rhs_surf_expr = newsurfrhs;
            
        else # time derivative, but no surface
            log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" = "*string(rhs_expr));
            
            (newlhs, newrhs) = reformat_for_stepper(lhs_expr, rhs_expr, config.stepper);
            
            log_entry("Weak form, modified for time stepping: "*string(newlhs)*" = "*string(newrhs));
            
            lhs_expr = newlhs;
            rhs_expr = newrhs;
        end
        
    end
    
    # Here we build a SymExpression for each of the pieces. 
    # This is an Expr tree that is passed to the code generator.
    if length(result_exprs) == 4 # has surface terms
        (lhs_symexpr, rhs_symexpr, lhs_surf_symexpr, rhs_surf_symexpr) = build_symexpressions(wfvars, lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr);
        log_entry("lhs volume symexpression:\n"*string(lhs_symexpr));
        log_entry("lhs surface symexpression:\n"*string(lhs_surf_symexpr));
        log_entry("rhs volume symexpression:\n"*string(rhs_symexpr));
        log_entry("rhs surface symexpression:\n"*string(rhs_surf_symexpr));
    else
        (lhs_symexpr, rhs_symexpr) = build_symexpressions(wfvars, lhs_expr, rhs_expr);
        log_entry("lhs symexpression:\n"*string(lhs_symexpr));
        log_entry("rhs symexpression:\n"*string(rhs_symexpr));
    end
    
    
    ########## This part simply makes a string for printing #############
    # make a string for the expression
    if typeof(lhs_expr[1]) <: Array
        lhsstring = "";
        rhsstring = "";
        for i=1:length(lhs_expr)
            lhsstring = lhsstring*"lhs"*string(i)*" = "*string(lhs_expr[i][1]);
            rhsstring = rhsstring*"rhs"*string(i)*" = "*string(rhs_expr[i][1]);
            for j=2:length(lhs_expr[i])
                lhsstring = lhsstring*" + "*string(lhs_expr[i][j]);
            end
            for j=2:length(rhs_expr[i])
                rhsstring = rhsstring*" + "*string(rhs_expr[i][j]);
            end
            if length(result_exprs) == 4
                for j=1:length(lhs_surf_expr)
                    lhsstring = lhsstring*" + surface("*string(lhs_surf_expr[j])*")";
                end
                for j=1:length(rhs_surf_expr)
                    rhsstring = rhsstring*" + surface("*string(rhs_surf_expr[j])*")";
                end
            end
            lhsstring = lhsstring*"\n";
            rhsstring = rhsstring*"\n";
        end
    else
        lhsstring = "lhs = "*string(lhs_expr[1]);
        rhsstring = "rhs = "*string(rhs_expr[1]);
        for j=2:length(lhs_expr)
            lhsstring = lhsstring*" + "*string(lhs_expr[j]);
        end
        for j=2:length(rhs_expr)
            rhsstring = rhsstring*" + "*string(rhs_expr[j]);
        end
        if length(result_exprs) == 4
            for j=1:length(lhs_surf_expr)
                lhsstring = lhsstring*" + surface("*string(lhs_surf_expr[j])*")";
            end
            for j=1:length(rhs_surf_expr)
                rhsstring = rhsstring*" + surface("*string(rhs_surf_expr[j])*")";
            end
        end
    end
    log_entry("Weak form, symbolic layer:\n"*string(lhsstring)*"\n"*string(rhsstring));
    ################ string ####################################
    
    # change symbolic layer into code layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "volume", config.solver_type, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "volume", config.solver_type, language, gen_framework);
    if length(result_exprs) == 4
        (lhs_surf_string, lhs_surf_code) = generate_code_layer(lhs_surf_symexpr, var, LHS, "surface", config.solver_type, language, gen_framework);
        (rhs_surf_string, rhs_surf_code) = generate_code_layer(rhs_surf_symexpr, var, RHS, "surface", config.solver_type, language, gen_framework);
        log_entry("Weak form, code layer: LHS = \n"*string(lhs_string)*"\nsurfaceLHS = \n"*string(lhs_surf_string)*" \nRHS = \n"*string(rhs_string)*"\nsurfaceRHS = \n"*string(rhs_surf_string));
    else
        log_entry("Weak form, code layer: LHS = \n"*string(lhs_string)*" \n  RHS = \n"*string(rhs_string));
    end
    
    #log_entry("Julia code Expr: LHS = \n"*string(lhs_code)*" \n  RHS = \n"*string(rhs_code), 3);
    
    if language == JULIA || language == 0
        args = "args";
        @makeFunction(args, lhs_code);
        set_lhs(var);
        
        @makeFunction(args, rhs_code);
        set_rhs(var);
        
        if length(result_exprs) == 4
            args = "args";
            @makeFunction(args, string(lhs_surf_code));
            set_lhs_surface(var);
            
            @makeFunction(args, string(rhs_surf_code));
            set_rhs_surface(var);
        end
        
    else
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

function fluxAndSource(var, fex, sex)
    flux(var, fex);
    source(var, sex);
end

function flux(var, fex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symfex = [];
        if !(length(var) == length(fex))
            printerr("Error in flux function: # of unknowns must equal # of equations. (example: flux([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symfex, Meta.parse((fex)[vi]));
        end
    else
        symfex = Meta.parse(fex);
        symvars = var.symbol;
    end
    
    log_entry("Making flux for variable(s): "*string(symvars));
    log_entry("flux, input: "*string(fex));
    
    # The parsing step
    (flhs_expr, frhs_expr) = sp_parse(symfex, symvars);
    # Note that this automatically puts a (-) on all rhs terms. These need to be reversed for FV.
    
    # Modify the expressions according to time stepper.
    # There is an assumed Dt(u) added which is the only time derivative.
    log_entry("flux, before modifying for time: "*string(flhs_expr)*" - "*string(frhs_expr));
    
    (newflhs, newfrhs) = reformat_for_stepper_fv_flux(flhs_expr, frhs_expr, config.stepper);
    
    log_entry("flux, modified for time stepping: "*string(newflhs)*" + "*string(newfrhs));
    
    flhs_expr = newflhs;
    frhs_expr = newfrhs;
    
    # Here we build a SymExpression for each of the pieces. 
    # This is passed to the code generator.
    (lhs_symexpr, rhs_symexpr) = build_symexpressions(symvars, flhs_expr, frhs_expr);
    
    log_entry("flux lhs symexpression:\n"*string(lhs_symexpr));
    log_entry("flux rhs symexpression:\n"*string(rhs_symexpr));
    
    # change symbolic layer into code layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "surface", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "surface", FV, language, gen_framework);
    log_entry("flux, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args";
        @makeFunction(args, string(lhs_code));
        set_lhs_surface(var);
        
        @makeFunction(args, string(rhs_code));
        set_rhs_surface(var);
        
    else
        set_lhs_surface(var, lhs_code);
        set_rhs_surface(var, rhs_code);
    end
end

function source(var, sex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symsex = [];
        if !(length(var) == length(sex))
            printerr("Error in source function: # of unknowns must equal # of equations. (example: source([a,b,c], [s1,s2,s3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symsex, Meta.parse((sex)[vi]));
        end
    else
        symsex = Meta.parse(sex);
        symvars = var.symbol;
    end
    
    log_entry("Making source for variable(s): "*string(symvars));
    log_entry("source, input: "*string(sex));
    
    # The parsing step
    (slhs_expr, srhs_expr) = sp_parse(symsex, symvars);
    # Note that this automatically puts a (-) on all rhs terms. These need to be reversed for FV.
    
    # Modify the expressions according to time stepper.
    # There is an assumed Dt(u) added which is the only time derivative.
    log_entry("source, before modifying for time: "*string(slhs_expr)*" - "*string(srhs_expr));
    
    (newslhs, newsrhs) = reformat_for_stepper_fv_source(slhs_expr, srhs_expr, config.stepper);
    
    log_entry("source, modified for time stepping: "*string(newslhs)*" + "*string(newsrhs));
    
    slhs_expr = newslhs;
    srhs_expr = newsrhs;
    
    # Here we build a SymExpression for each of the pieces. 
    # This is passed to the code generator.
    (lhs_symexpr, rhs_symexpr) = build_symexpressions(symvars, slhs_expr, srhs_expr);
    
    log_entry("source lhs symexpression:\n"*string(lhs_symexpr));
    log_entry("source rhs symexpression:\n"*string(rhs_symexpr));
    
    # change symbolic layer into code code_layer
    (lhs_string, lhs_code) = generate_code_layer(lhs_symexpr, var, LHS, "volume", FV, language, gen_framework);
    (rhs_string, rhs_code) = generate_code_layer(rhs_symexpr, var, RHS, "volume", FV, language, gen_framework);
    log_entry("source, code layer: \n  LHS = "*string(lhs_string)*" \n  RHS = "*string(rhs_string));
    
    if language == JULIA || language == 0
        args = "args";
        @makeFunction(args, string(lhs_code));
        set_lhs(var);
        
        @makeFunction(args, string(rhs_code));
        set_rhs(var);
        
    else
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

function exportCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if language == JULIA || language == 0
        file = open(filename*".jl", "w");
        println(file, "#=\nGenerated functions for "*project_name*"\n=#\n");
        for LorR in [LHS, RHS]
            if LorR == LHS
                codevol = bilinears;
                codesurf = face_bilinears;
            else
                codevol = linears;
                codesurf = face_linears;
            end
            for i=1:length(variables)
                var = string(variables[i].symbol);
                if !(codevol[i] === nothing)
                    func_name = LorR*"_volume_function_for_"*var;
                    println(file, "function "*func_name*"(args)");
                    println(file, codevol[i].str);
                    println(file, "end #"*func_name*"\n");
                else
                    println(file, "# No "*LorR*" volume set for "*var*"\n");
                end
                
                if !(codesurf[i] === nothing)
                    func_name = LorR*"_surface_function_for_"*var;
                    println(file, "function "*func_name*"(args)");
                    println(file, codesurf[i].str);
                    println(file, "end #"*func_name*"\n");
                else
                    println(file, "# No "*LorR*" surface set for "*var*"\n");
                end
            end
        end
        
        close(file);
        
    else
        # Should we export for other targets?
    end
end

function importCode(filename)
    # For now, only do this for Julia code because others are already output in code files.
    if language == JULIA || language == 0
        file = open(filename*".jl", "r");
        lines = readlines(file, keep=true);
        for LorR in [LHS, RHS]
            if LorR == LHS
                codevol = bilinears;
                codesurf = face_bilinears;
            else
                codevol = linears;
                codesurf = face_linears;
            end
            
            # Loop over variables and check to see if a matching function is present.
            for i=1:length(variables)
                # Scan the file for a pattern like
                #   function LHS_volume_function_for_u
                #       ...
                #   end #LHS_volume_function_for_u
                #
                # Set LHS/RHS, volume/surface, and the variable name
                var = string(variables[i].symbol);
                vfunc_name = LorR*"_volume_function_for_"*var;
                vfunc_string = "";
                sfunc_name = LorR*"_surface_function_for_"*var;
                sfunc_string = "";
                for st=1:length(lines)
                    if occursin("function "*vfunc_name, lines[st])
                        # s is the start of the function
                        for en=(st+1):length(lines)
                            if occursin("end #"*vfunc_name, lines[en])
                                # en is the end of the function
                                st = en; # update st
                                break;
                            else
                                vfunc_string *= lines[en];
                            end
                        end
                    elseif occursin("function "*sfunc_name, lines[st])
                        # s is the start of the function
                        for en=(st+1):length(lines)
                            if occursin("end #"*sfunc_name, lines[en])
                                # en is the end of the function
                                st = en; # update st
                                break;
                            else
                                sfunc_string *= lines[en];
                            end
                        end
                    end
                end # lines loop
                
                # Generate the functions and set them in the right places
                if vfunc_string == ""
                    println("Warning: While importing, no "*LorR*" volume function was found for "*var);
                else
                    @makeFunction("args", string(vfunc_string));
                    if LorR == LHS
                        set_lhs(variables[i]);
                    else
                        set_rhs(variables[i]);
                    end
                end
                if sfunc_string == ""
                    println("Warning: While importing, no "*LorR*" surface function was found for "*var);
                else
                    @makeFunction("args", string(sfunc_string));
                    if LorR == LHS
                        set_lhs_surface(variables[i]);
                    else
                        set_rhs_surface(variables[i]);
                    end
                end
                
            end # vars loop
        end
        
    else
        # TODO non-julia
    end
end

# This will either solve the problem or generate the code for an external target.
function solve(var, nlvar=nothing; nonlinear=false)
    if use_cachesim
        printerr("Use cachesim_solve(var) for generating cachesim output. Try again.");
        return nothing;
    end
    
    global time_stepper; # This should not be necessary. It will go away eventually
    
    # Generate files or solve directly
    if !(!generate_external && (language == JULIA || language == 0)) # if an external code gen target is ready
        if typeof(var) <: Array
            varind = var[1].index;
        else
            varind = var.index;
        end
        generate_all_files(var, bilinears[varind], face_bilinears[varind], linears[varind], face_linears[varind]);
        
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
        
        if config.solver_type == CG
            init_cgsolver();
            
            lhs = bilinears[varind];
            rhs = linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
                if (nonlinear)
                    if time_stepper.type == EULER_EXPLICIT || time_stepper.type == LSRK4
                        printerr("Warning: Use implicit stepper for nonlinear problem. Aborting.");
                        return;
                    end
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
                            #println("putting result "*string(compi+tmp)*":"*string(totalcomponents)*":end in var["*string(vi)*"].values["*string(compi)*"]")
                            var[vi].values[compi,:] = result[(compi+tmp):totalcomponents:end];
                        end
                        tmp = tmp + components;
                    end
                elseif length(result) > 1
                    components = length(var.symvar.vals);
                    for compi=1:components
                        var.values[compi,:] = result[compi:components:end];
                    end
                end
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == DG
            init_dgsolver();
            
            lhs = bilinears[varind];
            rhs = linears[varind];
            slhs = face_bilinears[varind];
            srhs = face_linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs, time_stepper));
				else
                	t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs, time_stepper));
				end
                # result is already stored in variables
            else
                # solve it!
				if (nonlinear)
                	t = @elapsed(result = DGSolver.nonlinear_solve(var, nlvar, lhs, rhs, slhs, srhs));
                else
                    t = @elapsed(result = DGSolver.linear_solve(var, lhs, rhs, slhs, srhs));
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
                            var[vi].values[compi,:] = result[(compi+tmp):totalcomponents:end];
                            tmp = tmp + 1;
                        end
                    end
                elseif length(result) > 1
                    components = length(var.symvar.vals);
                    for compi=1:components
                        var.values[compi,:] = result[compi:components:end];
                    end
                end
            end
            
            log_entry("Solved for "*varnames*".(took "*string(t)*" seconds)", 1);
            
        elseif config.solver_type == FV
            init_fvsolver();
            
            slhs = bilinears[varind];
            srhs = linears[varind];
            flhs = face_bilinears[varind];
            frhs = face_linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    time_stepper.dt = specified_dt;
				    time_stepper.Nsteps = specified_Nsteps;
                end
				if (nonlinear)
                	t = 0;
                    #TODO
				else
                	t = @elapsed(result = FVSolver.linear_solve(var, slhs, srhs, flhs, frhs, time_stepper));
				end
                # result is already stored in variables
            else
                # does this make sense?
            end
        end
    end

end

# When using cachesim, this will be used to simulate the solve.
function cachesimSolve(var, nlvar=nothing; nonlinear=false)
    if !(!generate_external && (language == JULIA || language == 0))
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
        log_entry("Generated cachesim ouput for "*varnames*".(took "*string(t)*" seconds)", 1);
    end
end

# Writes the values for each variable to filename in the given format.
function output_values(vars, filename; format="raw", ascii=false)
    available_formats = ["raw", "csv", "vtk"];
    if format == "vtk"
        output_values_vtk(vars, filename, ascii)
        log_entry("Writing values to file: "*filename*".vtu");
        
    else
        filename *= "."*format;
        file = open(filename, "w");
        log_entry("Writing values to file: "*filename);
        
        if format == "raw"
            output_values_raw(vars, file)
        elseif format == "csv"
            output_values_csv(vars, file)
        else
            println("Unknown output file format("*string(format)*"). Choose from: "*string(available_formats));
        end
        log_entry("Finished writing to "*filename);
        
        close(file);
    end
end

function finalize_femshop()
    # Finalize generation
    finalize_code_generator();
    if use_cachesim
        CachesimOut.finalize();
    end
    # anything else
    close_log();
    println("Femshop has completed.");
end

### Other specialized functions ###

function cachesim(use)
    log_entry("Using cachesim - Only cachesim output will be generated.", 1);
    global use_cachesim = use;
end

function morton_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive(grid_data, griddim, MORTON_ORDERING));
    log_entry("Reordered nodes to Morton. Took "*string(t)*" sec.", 2);
end

function morton_elements(griddim)
    global elemental_order = get_recursive_order(MORTON_ORDERING, config.dimension, griddim);
    log_entry("Reordered elements to Morton.", 2);
end

function hilbert_nodes(griddim)
    t = @elapsed(global grid_data = reorder_grid_recursive(grid_data, griddim, HILBERT_ORDERING));
    log_entry("Reordered nodes to Hilbert. Took "*string(t)*" sec.", 2);
end

function hilbert_elements(griddim)
    global elemental_order = get_recursive_order(HILBERT_ORDERING, config.dimension, griddim);
    log_entry("Reordered elements to Hilbert.", 2);
end

function tiled_nodes(griddim, tiledim)
    t = @elapsed(global grid_data = reorder_grid_tiled(grid_data, griddim, tiledim));
    log_entry("Reordered nodes to tiled. Took "*string(t)*" sec.", 2);
end

function tiled_elements(griddim, tiledim)
    global elemental_order = get_tiled_order(config.dimension, griddim, tiledim, true);
    log_entry("Reordered elements to tiled("*string(tiledim)*").", 2);
end

function ef_nodes()
    t = @elapsed(global grid_data = reorder_grid_element_first(grid_data, config.basis_order_min, elemental_order));
    log_entry("Reordered nodes to EF. Took "*string(t)*" sec.", 2);
end

function random_nodes(seed = 17)
    t = @elapsed(global grid_data = reorder_grid_random(grid_data, seed));
    log_entry("Reordered nodes to random. Took "*string(t)*" sec.", 2);
end

function random_elements(seed = 17)
    global elemental_order = random_order(mesh_data.nel, seed);
    log_entry("Reordered elements to random.", 2);
end