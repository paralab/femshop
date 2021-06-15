#=
This file contains all of the common interface functions.
Many of them simply call corresponding functions in Femshop.jl.
=#
export generateFor, useLog, domain, solverType, functionSpace, trialSpace, testSpace, 
        nodeType, timeStepper, setSteps, matrixFree, customOperator, customOperatorFile,
        mesh, exportMesh, variable, coefficient, parameter, testSymbol, boundary, addBoundaryID,
        refernecePoint, timeInterval, initial, weakForm, fluxAndSource, exportCode, importCode,
        solve, cachesimSolve, finalize_femshop, cachesim,
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
        include(lang);
        set_custom_gen_target(get_external_language_elements, generate_external_code_layer, generate_external_files, outputDirPath, filename, head=header);
    else
        if lang == DENDRO
            framew = DENDRO;
            lang = CPP;
        elseif lang == HOMG
            framew = HOMG;
            lang = MATLAB;
        end
        set_language(lang, outputDirPath, filename, framework=framew, head=header);
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
    if orderMax > orderMin && orderMax > 0
        config.p_adaptive = true;
        config.basis_order_min = orderMin;
        config.basis_order_max = orderMax;
    else
        config.basis_order_min = order;
        config.basis_order_max = order;
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

function matrixFree(;maxiters=100, tol=1e-6)
    config.linalg_matrixfree = true;
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
        mfile = open(m, "r");
        log_entry("Reading mesh file: "*m);
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
        newval = [swap_parameter_xyzt(val)];
        
    elseif typeof(val) <: Array
        newval = Array{Expr,1}(undef,length(val));
        for i=1:length(val)
            newval[i] = swap_parameter_xyzt(val[i]);
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
    
    result_exprs = sp_parse(wfex, wfvars);
    if length(result_exprs) == 4
        (lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr) = result_exprs;
    else
        (lhs_expr, rhs_expr) = result_exprs;
    end
    
    if typeof(lhs_expr) <: Tuple && length(lhs_expr) == 2 # has time derivative
        if length(result_exprs) == 4
            log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" + surface("*string(lhs_surf_expr)*") = "*string(rhs_expr)*" + surface("*string(rhs_surf_expr)*")");
            
            (newlhs, newrhs, newsurflhs, newsurfrhs) = reformat_for_stepper(lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr, Femshop.config.stepper);
            #TODO reformat surface terms
            
            log_entry("Weak form, modified for time stepping: "*string(newlhs)*" + surface("*string(newsurflhs)*") = "*string(newrhs)*" + surface("*string(newsurfrhs)*")");
            
            lhs_expr = newlhs;
            rhs_expr = newrhs;
            lhs_surf_expr = newsurflhs;
            rhs_surf_expr = newsurfrhs;
            
        else
            log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" = "*string(rhs_expr));
            
            (newlhs, newrhs) = reformat_for_stepper(lhs_expr, rhs_expr, Femshop.config.stepper);
            
            log_entry("Weak form, modified for time stepping: "*string(newlhs)*" = "*string(newrhs));
            
            lhs_expr = newlhs;
            rhs_expr = newrhs;
        end
        
    end
    
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
    
    # change symbolic layer into code layer
    lhs_code = generate_code_layer(lhs_expr, var, LHS);
    rhs_code = generate_code_layer(rhs_expr, var, RHS);
    if length(result_exprs) == 4
        lhs_surf_code = generate_code_layer_surface(lhs_surf_expr, var, LHS);
        rhs_surf_code = generate_code_layer_surface(rhs_surf_expr, var, RHS);
        log_entry("Weak form, code layer: LHS = "*string(lhs_code)*"\nsurfaceLHS = "*string(lhs_surf_code)*" \nRHS = "*string(rhs_code)*"\nsurfaceRHS = "*string(rhs_surf_code));
    else
        log_entry("Weak form, code layer: LHS = "*string(lhs_code)*" \n  RHS = "*string(rhs_code));
    end
    
    if language == JULIA || language == 0
        args = "args";
        @makeFunction(args, string(lhs_code));
        set_lhs(var);
        
        @makeFunction(args, string(rhs_code));
        set_rhs(var);
        
        if length(result_exprs) == 4
            args = "args";
            @makeFunction(args, string(lhs_surf_code));
            set_lhs_surface(var);
            
            @makeFunction(args, string(rhs_surf_code));
            set_rhs_surface(var);
        end
        
    elseif language == CPP
        # Don't need to generate any functions
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    elseif language == MATLAB
        # Don't need to generate any functions
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    elseif language == -1
        # Don't need to generate any functions
        set_lhs(var, lhs_code);
        set_rhs(var, rhs_code);
    end
end

function fluxAndSource(var, fex, sex)
    if typeof(var) <: Array
        # multiple simultaneous variables
        symvars = [];
        symfex = [];
        symsex = [];
        symdex = [];
        if !(length(var) == length(fex) && length(var) == length(sex))
            printerr("Error in flux function: # of unknowns must equal # of equations. (example: @flux([a,b,c], [f1,f2,f3]))");
        end
        for vi=1:length(var)
            push!(symvars, var[vi].symbol);
            push!(symfex, Meta.parse((fex)[vi]));
            push!(symsex, Meta.parse((sex)[vi]));
            push!(symdex, Meta.parse("Dt("*string(var[vi].symbol)*")"));
        end
    else
        symfex = Meta.parse(fex);
        symsex = Meta.parse(sex);
        symdex = Meta.parse("Dt("*string(var.symbol)*")");
        symvars = var.symbol;
    end
    
    log_entry("Making flux and source for variable(s): "*string(symvars));
    log_entry("flux, input: "*string(fex));
    log_entry("source, input: "*string(sex));
    
    # The parsing step
    (flhs_expr, frhs_expr) = sp_parse(symfex, symvars);
    (slhs_expr, srhs_expr) = sp_parse(symsex, symvars);
    (dlhs_expr, drhs_expr) = sp_parse(symdex, symvars);
    
    # Modify the expressions according to time stepper.
    # There is an assumed Dt(u) added which is the only time derivative.
    log_entry("time derivative, before modifying for time: "*string(dlhs_expr));
    log_entry("flux, before modifying for time: "*string(flhs_expr)*" - "*string(frhs_expr));
    log_entry("source, before modifying for time: "*string(slhs_expr)*" - "*string(srhs_expr));
    
    (newflhs, newfrhs, newslhs, newsrhs) = reformat_for_stepper_fv(dlhs_expr[1], flhs_expr, frhs_expr, slhs_expr, srhs_expr, Femshop.config.stepper);
    
    log_entry("flux, modified for time stepping: "*string(newflhs)*" - "*string(newfrhs));
    log_entry("source, modified for time stepping: "*string(newslhs)*" - "*string(newsrhs));
    
    flhs_expr = newflhs;
    frhs_expr = newfrhs;
    slhs_expr = newslhs;
    srhs_expr = newsrhs;
    
    # change symbolic layer into code layer
    flhs_code = generate_code_layer_fv(flhs_expr, var, LHS, "flux");
    frhs_code = generate_code_layer_fv(frhs_expr, var, RHS, "flux");
    slhs_code = generate_code_layer_fv(slhs_expr, var, LHS, "source");
    srhs_code = generate_code_layer_fv(srhs_expr, var, RHS, "source");
    Femshop.log_entry("flux, code layer: \n  LHS = "*string(flhs_code)*" \n  RHS = "*string(frhs_code));
    Femshop.log_entry("source, code layer: \n  LHS = "*string(slhs_code)*" \n  RHS = "*string(srhs_code));
    
    if Femshop.language == JULIA || Femshop.language == 0
        args = "args";
        @makeFunction(args, string(slhs_code));
        set_lhs(var);
        
        @makeFunction(args, string(srhs_code));
        set_rhs(var);
        
        @makeFunction(args, string(flhs_code));
        set_lhs_surface(var);
        
        @makeFunction(args, string(frhs_code));
        set_rhs_surface(var);
        
    else
        # not ready
        println("Not ready to generate FV code foe external target.")
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
    if !(gen_files === nothing && (language == JULIA || language == 0)) # if an external code gen target is ready
        # generate_main();
        if !(dendro_params === nothing)
            generate_all_files(bilinears[1], linears[1], parameters=dendro_params);
        else
            generate_all_files(bilinears[1], linears[1]);
        end
        
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
                    Femshop.time_stepper.dt = specified_dt;
				    Femshop.time_stepper.Nsteps = specified_Nsteps;
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
            
        elseif config.solver_type == DG
            init_dgsolver();
            
            lhs = bilinears[varind];
            rhs = linears[varind];
            slhs = face_bilinears[varind];
            srhs = face_linears[varind];
            
            if prob.time_dependent
                time_stepper = init_stepper(grid_data.allnodes, time_stepper);
                if use_specified_steps
                    Femshop.time_stepper.dt = specified_dt;
				    Femshop.time_stepper.Nsteps = specified_Nsteps;
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
                    Femshop.time_stepper.dt = specified_dt;
				    Femshop.time_stepper.Nsteps = specified_Nsteps;
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
        log_entry("Generated cachesim ouput for "*varnames*".(took "*string(t)*" seconds)", 1);
    end
end

function finalize_femshop()
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