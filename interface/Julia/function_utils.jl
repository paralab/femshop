#=
# function utils
=#
export GenFunction, @stringToFunction, add_genfunction, @makeFunction

#TODO: Move these global variables into Femshop.jl
genfunc_count = 0;
genfunctions = [];

# Stores everything we could need about a generated function.
struct GenFunction      # example:
    name::String        # "genfunction_7"
    args::String        # "x"
    str::String         # "sin(x)"
    expr::Expr          # expression: sin(x)
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
