#=
# function utils
=#
export GenFunction, @stringToFunction, add_genfunction, @makeFunction, @makeFunctions

# Stores everything we could need about a generated function.
struct GenFunction      # example:
    name::String        # "genfunction_7"
    args::String        # "x"
    str::String         # "sin(x)"
    expr                # expression: sin(x) NOTE: could be an Expr, Symbol, Number,
    func                # handle for: genfunction_7(x) = sin(x)
end

# Generates a function
macro stringToFunction(name, args, fun)
    return esc(quote
        eval(Meta.parse($name*"("*$args*")="*$fun));
    end)
end

# Adds the GenFunction to a global array of generated functions
function add_genfunction(genfun)
    global genfunc_count += 1;
    global genfunctions = [genfunctions ; genfun];
    log_entry("Generated function: "*genfun.name);
end

# Makes a GenFunction and adds it to the
# args is a string like "x,y,z"
# fun is a string like "sin(x)*y + 3*z"
macro makeFunction(args, fun)
    return esc(quote
        name = "genfunction_"*string(Femshop.genfunc_count);
        ex = Meta.parse($fun);
        nf = GenFunction(name, $args, $fun, ex, @stringToFunction(name, $args, $fun));
        add_genfunction(nf);
    end)
end

# Makes either: a constant number, a genfunction, or an array of genfunctions
macro makeFunctions(ex)
    return esc(quote
        nfuns = 0;
        args = "x=0,y=0,z=0,t=0";
        if typeof($ex) <: Array
            if typeof($ex[1]) <: Number
                nfuns = 0;
            else
                for i=1:length($ex)
                    if typeof($ex[i]) == String
                        @makeFunction(args, $ex[i]);
                        nfuns = nfuns + 1;
                    end
                end
            end
        else
            if typeof($ex) == String
                @makeFunction(args, $ex);
                nfuns = 1;
            elseif typeof($ex) <: Number
                nfuns = 0;
            end
        end
        nfuns
    end)
end
