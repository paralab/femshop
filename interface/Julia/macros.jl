#=
Macros for the interface.
Some of these may not be necessary.
Some may be better as regular functions.
We'll change as needed.
=#

# Initializes code generation module
macro language(lang, filename)
    return esc(:(@language($lang, $filename, "")));
end
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

# These set the configuration states
macro domain(dims) return :(@domain($dims, SQUARE, GRID)) end
macro domain(dims, geometry, mesh)
    return esc(quote
        Femshop.config.dimension = $dims;
        Femshop.config.geometry = $geometry;
        Femshop.config.mesh_type = $mesh
    end)
end

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

macro functionSpace(space)
    return esc(:(@functionSpace($space, 1, 1)))
end
macro functionSpace(space, N)
    return esc(:(@functionSpace($space, $N, $N)))
end
macro functionSpace(space, Nmin, Nmax)
    return esc(quote
        @trialFunction($space, $Nmin, $Nmax);
        @testFunction($space, $Nmin, $Nmax);
    end)
end

macro trialFunction(space, N)
    return esc(:(@trialFunction($space, $N, $N)))
end
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

macro testFunction(space, N)
    return esc(:(@testFunction($space, $N, $N)))
end
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

macro nodes(nodetype)
    return esc(quote
        Femshop.config.elemental_nodes = $nodetype;
    end)
end

macro order(N)
    return esc(:(@order($N, $N)));
end
macro order(Nmin, Nmax)
    return esc(quote
        if $Nmin == $Nmax
            Femshop.config.p_adaptive = false;
        else
            Femshop.config.p_adaptive = true;
        end
        Femshop.config.basis_order_min = $Nmin;
        Femshop.config.basis_order_max = $Nmax;
    end)
end

macro stepper(s) return esc(:(@stepper($s, 0))) end
macro stepper(s, cfl)
    return esc(quote
        set_stepper($s, $cfl);
    end)
end

macro matrixFree() return esc(:(@matrixFree(1000, 1e-6))) end
macro matrixFree(max, tol)
    return esc(quote
        set_matrix_free($max, $tol);
    end)
end

macro customOperator(s, handle)
    symb = string(s);
    return esc(quote
        varind = Femshop.var_count + 1;
        $s = Symbol($symb);
        add_custom_op($s, $handle);
        log_entry("Added custom operator: "*string($s));
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
macro mesh(m,N) return esc(quote @mesh($m,$N,1); end) end
macro mesh(m, N, bids)
    return esc(quote
        if $m == LINEMESH
            log_entry("Building simple line mesh with nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_line_mesh($N+1, $bids)));
            log_entry("Grid building took "*string(meshtime)*" seconds");
        elseif $m == QUADMESH
            log_entry("Building simple quad mesh with nx*nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_quad_mesh($N+1, $bids)));
            log_entry("Grid building took "*string(meshtime)*" seconds");
        elseif $m == HEXMESH
            log_entry("Building simple hex mesh with nx*nx*nx elements, nx="*string($N));
            meshtime = @elapsed(add_mesh(simple_hex_mesh($N+1, $bids)));
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
        global $var = Symbol($varsym);
        $var = Femshop.Variable($var, nothing, varind, $type, [], [], false);
        add_variable($var);
    end)
end

macro coefficient(c, val) return esc(:(@coefficient($c, SCALAR, $val))); end
macro coefficient(c, type, val)
    csym = string(c);
    return esc(quote
        global $c = Symbol($csym);
        nfuns = @makeFunctions($val); # if val is constant, nfuns will be 0
        $c = add_coefficient($c, $type, $val, nfuns);
    end)
end

macro boundary(var, bid, bc_type, bc_exp)
    return esc(quote
        nfuns = @makeFunctions($bc_exp);
        add_boundary_condition($var, $bid, $bc_type, $bc_exp, nfuns);
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
        global $var = Symbol($varsym);
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
        
        (lhs_expr, rhs_expr) = sp_parse(wfex, wfvars);
        
        if typeof(lhs_expr) <: Tuple && length(lhs_expr) == 2 # has time derivative
            log_entry("Weak form, symbolic layer: Dt("*string(lhs_expr[1])*") + "*string(lhs_expr[2])*" = "*string(rhs_expr));
            
            (newlhs, newrhs) = reformat_for_stepper(lhs_expr, rhs_expr, Femshop.config.stepper);
            
            log_entry("Weak form, modified for time stepping: "*string(newlhs)*" = "*string(newrhs));
            
            lhs_expr = newlhs;
            rhs_expr = newrhs;
        end
        
        # make a string for the expression
        if length(lhs_expr) > 1
            lhsstring = "";
            rhsstring = "";
            for i=1:length(lhs_expr)
                global lhsstring = lhsstring*"lhs"*string(i)*" = "*string(lhs_expr[i][1]);
                global rhsstring = rhsstring*"rhs"*string(i)*" = "*string(rhs_expr[i][1]);
                for j=2:length(lhs_expr[i])
                    global lhsstring = lhsstring*" + "*string(lhs_expr[i][j]);
                end
                for j=2:length(rhs_expr[i])
                    global rhsstring = rhsstring*" + "*string(rhs_expr[i][j]);
                end
                lhsstring = lhsstring*"\n";
                rhsstring = rhsstring*"\n";
            end
        else
            lhsstring = "lhs = "*string(lhs_expr[1][1]);
            rhsstring = "rhs = "*string(rhs_expr[1][1]);
            for j=2:length(lhs_expr[1])
                global lhsstring = lhsstring*" + "*string(lhs_expr[1][j]);
            end
            for j=2:length(rhs_expr[1])
                global rhsstring = rhsstring*" + "*string(rhs_expr[1][j]);
            end
        end
        log_entry("Weak form, symbolic layer:\n"*string(lhsstring)*"\n"*string(rhsstring));
        
        # change symbolic layer into code layer
        lhs_code = generate_code_layer(lhs_expr, $var, LHS);
        rhs_code = generate_code_layer(rhs_expr, $var, RHS);
        Femshop.log_entry("Weak form, code layer: LHS = "*string(lhs_code)*" \n  RHS = "*string(rhs_code));
        
        if Femshop.language == JULIA || Femshop.language == 0
            args = "args";
            @makeFunction(args, string(lhs_code));
            set_lhs($var);
            
            @makeFunction(args, string(rhs_code));
            set_rhs($var);
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

macro finalize()
    return esc(:(Femshop.finalize()));
end
