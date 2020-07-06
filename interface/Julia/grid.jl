#=
# Contains info about all nodes on the domain
=#
struct Grid
    allnodes::Array{Float64}    # All node coordinates
    bdry::Array{Int,2}     # Indices of boundary nodes for each BID
    bids::Array{Int,1}          # BID corresponding to rows of bdrynodes
    
    Grid(a,b,c) = new(a,b,c)
end