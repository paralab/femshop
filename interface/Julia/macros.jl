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
macro useLog() return :(@uselog("logFile")) end
macro useLog(name)
    return esc(quote
        if @isdefined(outputDirPath)
            @useLog($name, outputDirPath);
        else
            @useLog($name, pwd());
        end
    end)
end
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

macro mesh(file)
    return esc(quote
        # open the file and read the mesh data
        mfile = open($file, "r");
        log_entry("Reading mesh file: "*$file);
        add_mesh(read_mesh(mfile));
        close(mfile);
        # if mfile != nothing
        #     log_entry("Reading mesh file: "*$file);
        #     add_mesh(read_mesh(mfile));
        #     close(mfile);
        # else
        #     printerr("Error: couldn't open mesh file: "*$file);
        # end
        
    end)
end

macro solver(type)
    return esc(:(@solver($type, NODAL)));
end
macro solver(type, nodalORmodal)
    return esc(quote
        Femshop.config.solver_type = $type;
        Femshop.config.basis_type = $nodalORmodal;
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


macro finalize()
    return esc(:(Femshop.finalize()));
end
