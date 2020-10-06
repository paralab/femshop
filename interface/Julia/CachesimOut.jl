module CachesimOut

export init_cachesimout, take_available_addr_range, add_cachesim_array, cachesim_load, cachesim_store
export CachesimArray

#=
Want

l 2342
s 512 8
l 512 8

to be turned into

cs.load(2342)  # Loads one byte from address 2342
cs.store(512, length=8)  # Stores 8 bytes to addresses 512-519
cs.load(512, length=8)  # Loads from address 512 until (exclusive) 520 (eight bytes)
...
cs.force_write_back()
cs.print_stats()

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
7 J
8 detJ
9 w
10+ any other needed arrays
=#

addr_offset = 0;
geofac_set = false;
arrays = [];
output = open("cachesim_output.out", "w");

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

function init_cachesimout(N, refel, dof)
    # 1 A
    # 2 b
    # 3 u
    # 4 Ak
    # 5 bk
    # 6 Q
    # 7 RQ1
    # 8 RQ2
    # 9 RQ3
    # 10 TRQ1
    # 11 TRQ2
    # 12 TRQ3
    # 13 RD1
    # 14 RD2
    # 15 RD3
    # 17 wdetJ
    Nn = N*dof;
    Nd = refel.Np*dof;
    
    push!(arrays, CachesimArray(1, [Nn,Nn], 8));
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
    push!(arrays, CachesimArray(17, [Nd], 8));
end

function add_cachesim_array(sz, elsz)
    id = length(arrays)+1;
    push!(arrays, CachesimArray(id, sz, elsz));
    return id;
end

function cachesim_load(id, range=-1, rowoffset=0)
    a = arrays[id];
    alen = a.addr[2] - a.addr[1] + 1;
    if range == -1
        range = 1:alen;
    end
    
    for i=1:length(range)
        ad = a.addr[1] + rowoffset + (range[i]-1)*a.elsize;
        write_output(0,ad);
        #println("wrote load:"*string(ad));
    end
    
end

function cachesim_store(id, range=-1, rowoffset=0)
    a = arrays[id];
    alen = a.addr[2] - a.addr[1] + 1;
    if range == -1
        range = 1:alen;
    end
    
    for i=1:length(range)
        ad = a.addr[1] + rowoffset + (range[i]-1)*a.elsize;
        write_output(1,ad);
        #println("wrote store:"*string(ad));
    end
end

function write_output(type,addr)
    write(output,Int8(type));
    write(output,addr);
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

function finalize()
    close(output);
end

end# module