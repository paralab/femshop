#=
# C++ specific code generation 
# In general they will return arrays of strings representing lines.
# Line ending chars like ";" are included, but newline chars are not.
=#

# numbers: 2 -> "2"
# strings: "thing" -> "\"thing\""
# arrays: [1 2; 3 4] -> "{{1,2},{3,4}}"
function cpp_gen_string(v)
    if typeof(v) == String
        return "\""*v*"\"";
        
    elseif typeof(v) <: Number
        return string(v);
        
    elseif typeof(v) <: Array
        if ndims(v) == 1 # "{a, b, c}"
            n = length(v);
            str = "{";
            for i=1:n
                str = str*cpp_gen_string(v[i]);
                if i < n
                    str = str*", ";
                end
            end
            str = str*"}";
        elseif ndims(v) == 2 # "{{a, b}, {c, d}}
            (n,m) = size(v);
            str = "{";
            for i=1:n
                for j=1:m
                    str = str*cpp_gen_string(v[i,j]);
                    if i*j < n*m
                        str = str*", ";
                    end
                end
            end
            str = str*"}";
        end
        return str;
        
    elseif typeof(v) == GenFunction
        return v.name;
    else
        return string(v);
    end
end

#= Generates:
for(iterator=range[1]; iterator<range[2]; iterator+=step){
    (content)
}
# Note: content should be an array of lines to allow indentation
=#
function cpp_for_loop(indent, iterator, range, step, content)
    lines = [];
    inner_indent = indent*"    ";
    push!(lines, indent*"for("*iterator*" = "*string(range[1])*";"*iterator*" < "*string(range[2])*";"*iterator*" += "*string(step)*"){");
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"}")
    
    return lines;
end

#= Generates:
std::function<returntype(argtypes)> name = [captures](args){
    (content)
};
# Note: content should be an array of lines to allow indentation
=#
function cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content)
    lines = [];
    inner_indent = indent*"    ";
    
    arg = string(argtypes[1])*" "*string(args[1]);
    argtype = string(argtypes[1]);
    for i=2:length(args)
        arg = arg*", "*string(argtypes[i])*" "*string(args[i]);
        argtype = argtype*", "*string(argtypes[i]);
    end
    ret = string(ret);
    
    push!(lines, indent*"std::function<"*rettype*"("*argtype*")> "*name*" = ["*captures*"]("*arg*"){");
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"};")
    
    return lines;
end

# #= Generates:
# name = val;
# # Note, if val is empty, nothing is needed, returns ""
# # Type is ignored
# =#
# function cpp_declare_var(indent, name, type=nothing, val=nothing)
#     if val === nothing
#         return "";
#     else
#         return indent*name*" = "*string(val)*";";
#     end
# end

# #= Generates:
# name = zeros(size);
# # Note, size should be a tuple for multidimensional arrays, integer for vectors
# # Type is ignored
# =#
# function cpp_alloc_array(indent, name, size, type=nothing)
#     if typeof(size) <: Tuple
#         s = string(size[1]);
#         for i=2:length(size)
#             s = s*", "*size[i];
#         end
#     elseif typeof(size) <: Int
#         s = string(size)*", 1";
#     else
#         println("Error generating cpp code! Tried allocating array with bad size arg");
#         s = "1,1";
#     end
    
# end

# #= Generates:
# name(args)
# =#
# function cpp_function_call(indent, name, args)
#     arg = string(args[1]);
#     for i=2:length(args)
#         arg = arg*", "*args[i];
#     end
#     return indent*name*"("*arg*")";
# end
