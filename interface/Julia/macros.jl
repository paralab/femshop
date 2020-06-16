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

# ============= end config , begin prob ==============

macro mesh(file)
    return esc(quote
        # open the file and read the mesh data
        mfile = open($file, "r");
        log_entry("Reading mesh file: "*$file);
        add_mesh(read_mesh(mfile));
        close(mfile);
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
        $var = Femshop.Variable($var, varind, $type, [], [], false);
        add_variable($var);
    end)
end

macro boundary(var, bid, bc_type, bc_exp)
    return esc(quote
        if Femshop.config.dimension == 1
            args = "x,t";
        elseif Femshop.config.dimension == 2
            args = "x,y,t";
        elseif Femshop.config.dimension == 3
            args = "x,y,z,t";
        end
        @makeFunction(args, $bc_exp);
        add_boundary_condition($var, $bid, $bc_type);
    end)
end

macro timeInterval(t)
    return esc(quote
        Femshop.prob.time_dependent = true;
        Femshop.prob.end_time = $t;
    end)
end

macro initial(var, ic)
    return esc(quote
        if Femshop.config.dimension == 1
            args = "x";
        elseif Femshop.config.dimension == 2
            args = "x,y";
        elseif Femshop.config.dimension == 3
            args = "x,y,z";
        end
        @initial($var, args, $ic);
    end)
end
macro initial(var, args, ic)
    return esc(quote
        @makeFunction($args, $ic);
        add_initial_condition($var.index);
    end)
end

macro testFunction(var)
    varsym = string(var);
    return esc(quote
        $var = Symbol($varsym);
        add_test_function($var);
    end)
end

macro weakForm(var, ex)
    return esc(quote
        wfex = Meta.parse($ex);
        args = "; time=0";
        for v in Femshop.variables
            global args = args*", "*string(v.symbol)*"=0";
        end
        (lhs_expr, rhs_expr) = sp_parse(wfex, $var.symbol, Femshop.test_function_symbol);
        @makeFunction(args, string(lhs_expr));
        set_lhs($var);
        
        @makeFunction(args, string(rhs_expr));
        set_rhs($var);
        
        Femshop.log_entry("Set weak form for variable "*string($var.symbol)*" to: "*string(lhs_expr)*" = "*string(rhs_expr));
    end)
end

################# temporary
macro RHS(var, ex)
    return esc(quote
        args = "v, flux, time";
        if $var == u
            eq = "Femshop.solver.mass_inv_advective(v[:,:,2]) .- Femshop.solver.surface_int(flux[:,:,2])";
        else
            eq = "Femshop.solver.mass_inv_advective(v[:,:,1]) .- Femshop.solver.surface_int(flux[:,:,1])";
        end
        @makeFunction(args, eq);
        set_rhs($var);
    end)
end

macro LHS(var, ex)
    if typeof(ex) == Symbol
        return esc(quote set_lhs($var, 0) end);
    elseif ex.head === :call && ex.args[1] === :Dt
        #arg = ex.args[2];
        return esc(quote set_lhs($var, 0, lhs_time_deriv=true) end);
    end
    return esc(quote
        set_lhs($var, 0)
    end)
end

############################

macro finalize()
    return esc(:(Femshop.finalize()));
end
