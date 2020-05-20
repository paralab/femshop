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

# These set the configuration states
macro domain(dims, geometry, mesh)
    return esc(quote
        Femshop.config.dimension = $dims;
        Femshop.config.geometry = $geometry;
        Femshop.config.mesh_type = $mesh
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
    return esc(quote
        @trialFunction($space);
        @testFunction($space);
    end)
end

macro trialFunction(space)
    return esc(quote
        Femshop.config.trial_function = $space;
    end)
end

macro testFunction(space)
    return esc(quote
        Femshop.config.test_function = $space;
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
