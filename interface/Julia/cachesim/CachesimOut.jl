module CachesimOut

export init_cachesimout, take_available_addr_range, add_cachesim_array, add_cachesim_tmp_array, remove_all_tmp_arrays,
        cachesim_load_range, cachesim_store_range, cachesim_load, cachesim_store, build_cache_level, build_cache
export CachesimArray

use_lib = false;
if use_lib
    include("Jpycachesim.jl");
    using .Jpycachesim
    
else
    println("Not using the cachesim library. Change use_lib in CachesimOut.jl to use.");
    println("Instead, writing access sequence to cachesim_output.out");
    println("Analyze with cachesim_read.py");
    output = open("cachesim/cachesim_output.out","w");
end

addr_offset = 0;
geofac_set = false;
arrays = [];
tmparrays = [];
tmpflags = [];
N=0; # number of global nodes
Np=0;# number of elemental nodes
Nn=0;# number of global dofs
Nd=0;# number of elemental dofs
Nv=0;# number of global vertices

# A special nonexistant array that also has an identifier, fake address range, size
mutable struct CachesimArray
    id::Int             # identifying unique number
    addr::Array{Int,1}  # fake address range
    size::Array{Int,1}  # size of the array(in elements)
    elsize::Int         # size of the elements in the array(in bytes)
    bytesize::Int       # total number of bytes in array
    
    CachesimArray(i, sz, elsz) = new(
        i,
        take_available_addr_range(sz, elsz),
        size_tuple_to_array(sz),
        elsz
    )
end

function init_cachesimout(n, refel, dof, vars)
    global N = n; # number of global nodes
    global Nn = N*dof; # number of global dof
    global Nd = refel.Np*dof; # elemental dof
    global Np = refel.Np; # number of elemental nodes
    
    #=
    The arrays are indexed as such:
    1 A
    2 b
    3 u
    4 Ak
    5 bk
    6 Q
    7 Q1
    8 Q2
    9 Q3
    10 D1
    11 D2
    12 D3
    13 Jx
    14 Jy
    15 Jz
    16 detJ
    17 w
    18+ variable values, then everything else
    =#
    push!(arrays, CachesimArray(1, [1], 8));
    push!(arrays, CachesimArray(2, [Nn], 8));
    push!(arrays, CachesimArray(3, [Nn], 8));
    push!(arrays, CachesimArray(4, [Nd, Nd], 8));
    push!(arrays, CachesimArray(5, [Nd], 8));
    push!(arrays, CachesimArray(6, size(refel.Q), 8));
    push!(arrays, CachesimArray(7, size(refel.Q), 8));
    push!(arrays, CachesimArray(8, size(refel.Q), 8));
    push!(arrays, CachesimArray(9, size(refel.Q), 8));
    push!(arrays, CachesimArray(10, size(refel.Q), 8));
    push!(arrays, CachesimArray(11, size(refel.Q), 8));
    push!(arrays, CachesimArray(12, size(refel.Q), 8));
    push!(arrays, CachesimArray(13, size(refel.Q), 8));
    push!(arrays, CachesimArray(14, size(refel.Q), 8));
    push!(arrays, CachesimArray(15, size(refel.Q), 8));
    push!(arrays, CachesimArray(16, [Nd], 8));
    push!(arrays, CachesimArray(17, [Nd], 8));
    
    for i=1:length(vars)
        push!(arrays, CachesimArray(17+i, size(vars[i].values), 8));
    end
    
    if use_lib
        pcs_get_cachesim_from_file("cachesim/cachedef"); # Just to test things. Should be set later by user.
    else
        
    end
    
end

function build_cache_level(level, sets, ways, cl_size, policy)
    if use_lib
        label = "L"*string(level);
        return Cache(label, sets, ways, cl_size, policy);
    else
        println("Can't build cache without library.");
    end
end

function build_cache(levels)
    if use_lib
        return pcs_build_cache(levels);
    else
        println("Can't build cache without library.");
    end
end

function add_cachesim_array(sz, elsz)
    id = length(arrays)+1;
    push!(arrays, CachesimArray(id, sz, elsz));
    return id;
end

function add_cachesim_tmp_array()
    # search the flags for an open tmp array
    id = -1;
    for i=1:length(tmpflags)
        if tmpflags[i] == false
            stupidjulia = id;
            id = tmparrays[i];
            tmpflags[i] = true;
            break;
        end
    end
    # get a new array if needed
    if id < 0
        id = add_cachesim_array([Np],8);
        push!(tmparrays, id);
        push!(tmpflags, true);
    end
    return id;
end

function remove_tmp_array(id)
    for i=1:length(tmparrays)
        if tmparrays[i] == id
            tmpflags[i] = false;
            break;
        end
    end
end

function remove_all_tmp_arrays()
    for i=1:length(tmparrays)
        tmpflags[i] = false;
    end
end

function cachesim_load_range(id, rows = [], cols = [])
    if length(rows) == 0
        rows = 1:arrays[id].size[1];
    end
    if length(cols) == 0 && length(arrays[id].size) > 1
        cols = 1:arrays[id].size[2];
    end
    if length(cols) > 0
        for j=1:length(cols)
            for i=1:length(rows)
                cachesim_load(id, rows[i], cols[j]);
            end
        end
    else
        for i=1:length(rows)
            cachesim_load(id, rows[i]);
        end
    end
end
function cachesim_store_range(id, rows = [], cols = [])
    if length(rows) == 0
        rows = 1:arrays[id].size[1];
    end
    if length(cols) == 0 && length(arrays[id].size) > 1
        cols = 1:arrays[id].size[2];
    end
    if length(cols) > 0
        for j=1:length(cols)
            for i=1:length(rows)
                cachesim_store(id, rows[i], cols[j]);
            end
        end
    else
        for i=1:length(rows)
            cachesim_store(id, rows[i]);
        end
    end
end

function cachesim_load(id, row, col=0)
    addr = ind_to_addr(id, row, col);
    if use_lib
        pcs_load(addr,8);
    else
        write(output,Int8(0));
        write(output,addr);
    end
end

function cachesim_store(id, row, col=0)
    addr = ind_to_addr(id, row, col);
    if use_lib
        pcs_store(addr,8);
    else
        write(output,Int8(1));
        write(output,addr);
    end
end

function cachesim_load_sparse(id, range=-1, rowoffset=0)
    a = arrays[id];
    alen = a.addr[2] - a.addr[1] + 1;
    if range == -1
        range = 1:alen;
    end
    
    for i=1:length(range)
        ad = a.addr[1] + rowoffset + (range[i]-1)*a.elsize;
        if use_lib
            pcs_load(ad,8);
        else
            write(output,Int8(0));
            write(output,ad);
        end
    end
    
end

function cachesim_store_sparse(id, range=-1, rowoffset=0)
    a = arrays[id];
    alen = a.addr[2] - a.addr[1] + 1;
    if range == -1
        range = 1:alen;
    end
    
    for i=1:length(range)
        ad = a.addr[1] + rowoffset + (range[i]-1)*a.elsize;
        if use_lib
            pcs_store(ad,8);
        else
            write(output,Int8(1));
            write(output,ad);
        end
    end
end

function take_available_addr_range(sz, elsz)
    range = [0,0];
    totalsz = sz[1];
    for i=2:length(sz)
        totalsz = totalsz*sz[i];
    end
    bytesz = totalsz*elsz;
    
    global addr_offset;
    range[1] = addr_offset;
    range[2] = range[1] + bytesz - 1;
    addr_offset += bytesz;
    
    return range;
end

function size_tuple_to_array(sz)
    if typeof(sz) <: Array
        return sz;
    elseif typeof(sz) <: Tuple
        ar = [];
        for i=1:length(sz)
            push!(ar,sz[i]);
        end
        return ar;
    end
end

function ind_to_addr(id, row, col=0)
    start = arrays[id].addr[1];
    elsize = arrays[id].elsize;
    if col == 0
        return start + (row-1)*elsize;
    else
        colsize = arrays[id].size[1];
        return start + (col-1)*colsize + (row-1)*elsize;
    end
end

function finalize()
    if use_lib
        pcs_print_stats();
    else
        close(output);
    end
end

end# module