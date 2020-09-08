#=
# Contains info about all nodes on the domain
=#
struct Grid
    allnodes::Array{Float64}    # All node coordinates
    bdry::Array{Array{Int,1},1} # Indices of boundary nodes for each BID (bdry[bid][nodes])
    bids::Array{Int,1}          # BID corresponding to rows of bdrynodes
    
    Grid(a,b,c) = new(a,b,c)
end