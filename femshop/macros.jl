#=
Macros for the interface.
=#

# Initializes code generation module
macro language(lang, filename) return esc(:(@language($lang, $filename, ""))); end
macro language(lang, filename, header)
    return esc(quote
        language = $lang;
        outputDirPath = pwd()*"/"*uppercasefirst($filename);
        if !isdir(outputDirPath)
            mkdir(outputDirPath);
        end
        outputFileName = $filename;
        set_language($lang, outputDirPath, outputFileName, $header);
    end)
end

# Optionally create a log file
macro useLog() return :(@uselog(Femshop.project_name)) end
macro useLog(name) return esc(:(@useLog($name, Femshop.output_dir))) end
macro useLog(name, dir)
    return esc(quote
        init_log($name, $dir);
    end)
end

# Sets dimension, domain shape, discretization type
macro domain(dims) return :(@domain($dims, SQUARE, UNIFORM_GRID)) end
macro domain(dims, geometry, mesh)
    return esc(quote
        Femshop.config.dimension = $dims;
        Femshop.config.geometry = $geometry;
        Femshop.config.mesh_type = $mesh
    end)
end

# Set the solver type: CG, etc.
macro solver(type)
    return esc(quote
        Femshop.config.solver_type = $type;
        if $type == DG
            set_solver(Femshop.DGSolver);
        elseif $type == CG
            set_solver(Femshop.CGSolver);
        end
    end)
end

# Set basis function type and order
macro functionSpace(space) return esc(:(@functionSpace($space, 1, 1))) end
macro functionSpace(space, N) return esc(:(@functionSpace($space, $N, $N))) end
macro functionSpace(space, Nmin, Nmax)
    return esc(quote
        @trialFunction($space, $Nmin, $Nmax);
        @testFunction($space, $Nmin, $Nmax);
    end)
end

macro trialFunction(space, N) return esc(:(@trialFunction($space, $N, $N))) end
macro trialFunction(space, Nmin, Nmax)
    return esc(quote
        Femshop.config.trial_function = $space;
        if $Nmin == $Nmax
            Femshop.config.p_adaptive = false;
        else
            Femshop.config.p_adaptive = true;
        end
        Femshop.config.basis_order_min = $Nmin;
        Femshop.config.basis_order_max = $Nmax;
    end)
end

macro testFunction(space, N) return esc(:(@testFunction($space, $N, $N))) end
macro testFunction(space, Nmin, Nmax)
    return esc(quote
        Femshop.config.test_function = $space;
        if $Nmin == $Nmax
            Femshop.config.p_adaptive = false;
        else
            Femshop.config.p_adaptive = true;
        end
        Femshop.config.basis_order_min = $Nmin;
        Femshop.config.basis_order_max = $Nmax;
    end)
end

# Sets elemental node locations: Lobatto, etc.
macro nodes(nodetype)
    return esc(quote
        Femshop.config.elemental_nodes = $nodetype;
    end)
end

#### Not needed ####
# macro order(N)
#     return esc(:(@order($N, $N)));
# end
# macro order(Nmin, Nmax)
#     return esc(quote
#         if $Nmin == $Nmax
#             Femshop.config.p_adaptive = false;
#         else
#             Femshop.config.p_adaptive = true;
#         end
#         Femshop.config.basis_order_min = $Nmin;
#         Femshop.config.basis_order_max = $Nmax;
#     end)
# end

# Sets time stepper type 
macro stepper(s) return esc(:(@stepper($s, 0))) end
macro stepper(s, cfl)
    return esc(quote
        set_stepper($s, $cfl);
    end)
end

# Sets time stepper dt and Nsteps to specified values
macro setSteps(dt, steps)
    return esc(quote
        set_specified_steps($dt, $steps);
    end)
end

# Selects matrix free
macro matrixFree() return esc(:(@matrixFree(1000, 1e-6))) end
macro matrixFree(max, tol)
    return esc(quote
        set_matrix_free($max, $tol);
    end)
end

# Adds a custom operator
macro customOperator(s, handle)
    symb = string(s);
    return esc(quote
        $s = Symbol($symb);
        Femshop.add_custom_op($s, $handle);
        log_entry("Added custom operator: "*string($s));
    end)
end

# Adds a set of custom operators defined in a file
macro customOperatorFile(file)
    return esc(quote
        Femshop.add_custom_op_file($file);
        log_entry("Added custom operators from file: "*string($file));
    end)
end

# ============= end config , begin prob ==============

macro mesh(m)
    return esc(quote
        if $m == LINEMESH || $m == QUADMESH || $m == HEXMESH
            @mesh($m, 5, 1)
        else
            # open the file and read the mesh data
            mfile = open($m, "r");
            log_entry("Reading mesh file: "*$m);
            add_mesh(read_mesh(mfile));
            close(mfile);
        end
        
    end)
end
macro mesh(m,N) return esc(quote @mesh($m,$N,1,[0,1]); end) end
macro mesh(m,N,bids) return esc(quote @mesh($m,$N,$bids,[0,1]); end) end
macro mesh(m, N, bids, interval)
    return esc(quote
        if $m == LINEMESH
            log_entry("Building simple line mesh with nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_line_mesh($N.+1, $bids, $interval)));
            log_entry("Grid building took "*string(meshtime)*" seconds");
        elseif $m == QUADMESH
            log_entry("Building simple quad mesh with nx*nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_quad_mesh($N.+1, $bids, $interval)));
            log_entry("Grid building took "*string(meshtime)*" seconds");
        elseif $m == HEXMESH
            log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_hex_mesh($N.+1, $bids, $interval)));
            log_entry("Grid building took "*string(meshtime)*" seconds");
        end
    end)
end

macro outputMesh(file) return esc(:(@outputMesh($file, MSH_V2))); end
macro outputMesh(file, format)
    return esc(quote
        # open the file to write to
        mfile = open($file, "w");
        log_entry("Writing mesh file: "*$file);
        output_mesh(mfile, $format);
        close(mfile);
    end)
end

macro variable(var) return esc(:(@variable($var, SCALAR))); end
macro variable(var, type)
    varsym = string(var);
    return esc(quote
        varind = Femshop.var_count + 1;
        $var = Symbol($varsym);
        $var = Femshop.Variable($var, nothing, varind, $type, [], [], false);
        add_variable($var);
    end)
end

macro coefficient(c, val) return esc(:(@coefficient($c, SCALAR, $val))); end
macro coefficient(c, type, val)
    csym = string(c);
    return esc(quote
        $c = Symbol($csym);
        nfuns = @makeFunctions($val); # if val is constant, nfuns will be 0
        $c = add_coefficient($c, $type, $val, nfuns);
    end)
end

macro parameter(p,val) return esc(:(@parameter($p, SCALAR, $val))); end
macro parameter(p, type, val)
    psym = string(p);
    return esc(quote
        $p = Symbol($psym);
        if !(@isdefined(parameter_coefficients_defined) && length(Femshop.parameters)>0)
            @coefficient(parameterCoefficientForx, "x")
            @coefficient(parameterCoefficientFory, "y")
            @coefficient(parameterCoefficientForz, "z")
            @coefficient(parameterCoefficientFort, "t")
            parameter_coefficients_defined = true;
        end
        if typeof($val) <: Number
            newval = [$val];
        elseif typeof($val) == String
            # newval will be an array of expressions
            # search for x,y,z,t symbols and replace with special coefficients like parameterCoefficientForx
            newval = [Femshop.swap_parameter_xyzt(Meta.parse($val))];
            
        elseif typeof($val) <: Array
            newval = Array{Expr,1}(undef,length(val));
            for i=1:length($val)
                newval[i] = Femshop.swap_parameter_xyzt(Meta.parse($val[i]));
            end
            newval = reshape(newval,size($val));
        else
            println("Error: use strings to define parameters");
            newval = 0;
        end
        
        $p = add_parameter($p, $type, newval);
    end)
end

macro boundary(var, bid, bc_type) return esc(:(@boundary($var,$bid,$bc_type, 0))); end
macro boundary(var, bid, bc_type, bc_exp)
    return esc(quote
        nfuns = @makeFunctions($bc_exp);
        add_boundary_condition($var, $bid, $bc_type, $bc_exp, nfuns);
    end)
end

macro referencePoint(var, pos, val)
    return esc(quote
        add_reference_point($var, $pos, $val);
    end)
end

macro timeInterval(t)
    return esc(quote
        Femshop.prob.time_dependent = true;
        if Femshop.time_stepper === nothing
            @stepper(EULER_IMPLICIT);
        end
        Femshop.prob.end_time = $t;
    end)
end

macro initial(var, ic)
    return esc(quote
        nfuns = @makeFunctions($ic);
        add_initial_condition($var.index, $ic, nfuns);
    end)
end

macro testSymbol(var) return esc(:(@testSymbol($var, SCALAR))); end
macro testSymbol(var, type)
    varsym = string(var);
    return esc(quote
        $var = Symbol($varsym);
        add_test_function($var, $type);
    end)
end

macro weakForm(var, ex)
    return esc(quote
        using LinearAlgebra
        
        if typeof($var) <: Array
            # multiple simultaneous variables
            wfvars = [];
            wfex = [];
            if !(length($var) == length($ex))
                printerr("Error in weak form: # of unknowns must equal # of equations. (example: @weakform([a,b,c], [f1,f2,f3]))");
            end
            for vi=1:length($var)
                push!(wfvars, $var[vi].symbol);
                push!(wfex, Meta.parse(($ex)[vi]));
            end
        else
            wfex = Meta.parse($ex);
            wfvars = $var.symbol;
        end
        
        log_entry("Making weak form for variable(s): "*string(wfvars));
        log_entry("Weak form, input: "*string($ex));
        
        result_exprs = sp_parse(wfex, wfvars);
        if length(result_exprs) == 4
            (lhs_expr, rhs_expr, lhs_surf_expr, rhs_surf_expr) = result_exprs;
        else
            (lhs_expr, rhs_expr) = result_exprs;
        end
        
        if typeof(lhs_expr) <: Tuple && length(lhs_expr) == 2 # has time derivative
            if length(result_exprs) == 4
                log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" + surface("*string(lhs_surf_expr)*") = "*string(rhs_expr)*" + surface("*string(rhs_surf_expr)*")");
                
                (newlhs, newrhs) = reformat_for_stepper(lhs_expr, rhs_expr, Femshop.config.stepper);
                #TODO reformat surface terms
                
                log_entry("Weak form, modified for time stepping: "*string(newlhs)*" + surface("*string(lhs_surf_expr)*") = "*string(newrhs)*" + surface("*string(rhs_surf_expr)*")");
                
                lhs_expr = newlhs;
                rhs_expr = newrhs;
                
            else
                log_entry("Weak form, before modifying for time: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" = "*string(rhs_expr));
                
                (newlhs, newrhs) = reformat_for_stepper(lhs_expr, rhs_expr, Femshop.config.stepper);
                
                log_entry("Weak form, modified for time stepping: "*string(newlhs)*" = "*string(newrhs));
                
                lhs_expr = newlhs;
                rhs_expr = newrhs;
            end
            
        end
        
        # make a string for the expression
        if length(lhs_expr) > 1
            global lhsstring = "";
            global rhsstring = "";
            for i=1:length(lhs_expr)
                global lhsstring = lhsstring*"lhs"*string(i)*" = "*string(lhs_expr[i][1]);
                global rhsstring = rhsstring*"rhs"*string(i)*" = "*string(rhs_expr[i][1]);
                for j=2:length(lhs_expr[i])
                    global lhsstring = lhsstring*" + "*string(lhs_expr[i][j]);
                end
                for j=2:length(rhs_expr[i])
                    global rhsstring = rhsstring*" + "*string(rhs_expr[i][j]);
                end
                if length(result_exprs) == 4
                    for j=2:length(lhs_surf_expr[i])
                        global lhsstring = lhsstring*" + surface("*string(lhs_surf_expr[i][j])*")";
                    end
                    for j=2:length(rhs_surf_expr[i])
                        global rhsstring = rhsstring*" + surface("*string(rhs_surf_expr[i][j])*")";
                    end
                end
                lhsstring = lhsstring*"\n";
                rhsstring = rhsstring*"\n";
            end
        else
            global lhsstring = "lhs = "*string(lhs_expr[1][1]);
            global rhsstring = "rhs = "*string(rhs_expr[1][1]);
            for j=2:length(lhs_expr[1])
                global lhsstring = lhsstring*" + "*string(lhs_expr[1][j]);
            end
            for j=2:length(rhs_expr[1])
                global rhsstring = rhsstring*" + "*string(rhs_expr[1][j]);
            end
            if length(result_exprs) == 4
                for j=2:length(lhs_surf_expr[1])
                    global lhsstring = lhsstring*" + surface("*string(lhs_surf_expr[1][j])*")";
                end
                for j=2:length(rhs_surf_expr[1])
                    global rhsstring = rhsstring*" + surface("*string(rhs_surf_expr[1][j])*")";
                end
            end
        end
        log_entry("Weak form, symbolic layer:\n"*string(lhsstring)*"\n"*string(rhsstring));
        
        # change symbolic layer into code layer
        lhs_code = generate_code_layer(lhs_expr, $var, LHS);
        rhs_code = generate_code_layer(rhs_expr, $var, RHS);
        if length(result_exprs) == 4
            lhs_surf_code = generate_code_layer_surface(lhs_surf_expr, $var, LHS);
            rhs_surf_code = generate_code_layer_surface(rhs_surf_expr, $var, RHS);
            Femshop.log_entry("Weak form, code layer: LHS = "*string(lhs_code)*"\n  + surface("*string(lhs_surf_code)*") \n  RHS = "*string(rhs_code)*"\n  + surface("*string(rhs_surf_code)*")");
        else
            Femshop.log_entry("Weak form, code layer: LHS = "*string(lhs_code)*" \n  RHS = "*string(rhs_code));
        end
        
        if Femshop.language == JULIA || Femshop.language == 0
            args = "args";
            @makeFunction(args, string(lhs_code));
            set_lhs($var);
            
            @makeFunction(args, string(rhs_code));
            set_rhs($var);
            
            if length(result_exprs) == 4
                args = "args";
                @makeFunction(args, string(lhs_surf_code));
                set_lhs_surface($var);
                
                @makeFunction(args, string(rhs_surf_code));
                set_rhs_surface($var);
            end
            
        elseif Femshop.language == CPP
            # Don't need to generate any functions
            set_lhs($var, lhs_code);
            set_rhs($var, rhs_code);
        elseif Femshop.language == MATLAB
            # Don't need to generate any functions
            set_lhs($var, lhs_code);
            set_rhs($var, rhs_code);
        end
        
    end)
end

# macro PDE(var, ex)
#     return esc(quote
#         # do the thing
#     end)
# end

macro finalize()
    return esc(:(Femshop.finalize()));
end
