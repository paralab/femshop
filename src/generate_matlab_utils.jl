#=
# Matlab specific code generation 
# In general they will return arrays of strings representing lines.
# Line ending chars like ";" are included, but newline chars are not.
=#

# numbers: 2 -> "2"
# strings: "thing" -> "'thing'"
# arrays: [1 2; 3 4] -> "[1 2; 3 4]"
function matlab_gen_string(v)
    if typeof(v) == String
        return "'"*v*"'";
        
    elseif typeof(v) <: Number
        return string(v);
        
    elseif typeof(v) <: Array
        if ndims(v) == 1 # "[a; b; c]"
            n = length(v);
            str = "[";
            for i=1:n
                str = str*matlab_gen_string(v[i]);
                if i < n
                    str = str*"; ";
                end
            end
            str = str*"]";
        elseif ndims(v) == 2 # "[a b ; c d]
            (n,m) = size(v);
            str = "[";
            for i=1:n
                for j=1:m
                    str = str*matlab_gen_string(v[i,j])*" ";
                end
                if i < n
                    str = str*"; ";
                end
            end
            str = str*"]";
        end
        return str;
        
    elseif typeof(v) == GenFunction
        return v.name;
    else
        return string(v);
    end
end

#= Generates:
for iterator = range[1]:step:range[2]
    (content)
end
# Note: content should be an array of lines to allow indentation
=#
function matlab_for_loop(indent, iterator, range, step, content)
    lines = [];
    inner_indent = indent*"    ";
    push!(lines, indent*"for "*iterator*" = "*string(range[1])*":"*string(step)*":"*string(range[2]));
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"end")
    
    return lines;
end

#= Generates:
function [rets] = name(args)
    (content)
end
# Note: content should be an array of lines to allow indentation
=#
function matlab_function_def(indent, name, args, rets, content)
    lines = [];
    inner_indent = indent*"    ";
    arg = string(args[1]);
    for i=2:length(args)
        arg = arg*", "*string(args[i]);
    end
    ret = string(rets[1]);
    for i=2:length(rets)
        ret = ret*", "*string(rets[i]);
    end
    push!(lines, indent*"function ["*ret*"] = "*name*"("*arg*")");
    for i=1:length(content)
        push!(lines, inner_indent*content[i]);
    end
    push!(lines, indent*"end")
    
    return lines;
end

#= Generates:
name = val;
# Note, if val is empty, nothing is needed, returns ""
# Type is ignored
=#
function matlab_declare_var(indent, name, type=nothing, val=nothing)
    if val === nothing
        return "";
    else
        return indent*name*" = "*string(val)*";";
    end
end

#= Generates:
name = zeros(size);
# Note, size should be a tuple for multidimensional arrays, integer for vectors
# Type is ignored
=#
function matlab_alloc_array(indent, name, size, type=nothing)
    if typeof(size) <: Tuple
        s = string(size[1]);
        for i=2:length(size)
            s = s*", "*size[i];
        end
    elseif typeof(size) <: Int
        s = string(size)*", 1";
    else
        println("Error generating matlab code! Tried allocating array with bad size arg");
        s = "1,1";
    end
    
end

#= Generates:
name(args)
=#
function matlab_function_call(indent, name, args)
    arg = string(args[1]);
    for i=2:length(args)
        arg = arg*", "*args[i];
    end
    return indent*name*"("*arg*")";
end

# Returns a string like "fread(f, [3, 15], 'double')" for an array with size (3,15) and type Float64
function matlab_fread(A)
    if typeof(A) <: Array
        if length(A) == 0
            return "0";
        end
        if typeof(A[1]) == Int
            typ = "\'int64\'";
        else
            typ = "\'double\'";
        end
        if length(size(A)) == 1
            sz = "["*string(length(A))*",1]";
        else
            sz = "["*string(size(A,1))*","*string(size(A,2))*"]";
        end
        return "fread(f, "*sz*", "*typ*")"
    else
        if typeof(A) == Int64
            typ = "\'int64\'";
        else
            typ = "\'double\'";
        end
        return "fread(f, [1], "*typ*")"
    end
end

# produces code to read a binary struct into matlab
# The struct is labeled with the name and has the same fieldnames as s
function matlab_struct_reader(name, s)
    code = "";
    for fn in fieldnames(typeof(s))
        f = getfield(s,fn);
        if typeof(f) <: Array && length(f) > 0 && typeof(f[1]) <: Array
            code = code * name * "." * string(fn) * " = cell([" * string(length(f)) * ", 1]);\n";
            for i=1:length(f)
                code = code * name * "." * string(fn) * "{" * string(i) * "} = " * matlab_fread(f[i]) * ";\n";
            end
        else
            code = code * name * "." * string(fn) * " = " * matlab_fread(f) * ";\n";
        end
    end
    return code;
end