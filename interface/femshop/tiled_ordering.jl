#=
# Builds a tiled ordering of a 2D or 3D grid
# Tiles are tiledim in size (except edge tiles which could be smaller)
# Grid has griddim size
# Returns the global order in which objects shall be indexed.
# example: 2D, griddim=(4,4), tiledim=(3,3)
#     13  14  15  16
#      7   8   9  12
#      4   5   6  11
#      1   2   3  10
# Invert=true will give  [1, 2, 3, 5, 6, 7 ...]   (encode?)
# Invert=false will give [1, 2, 3, 10, 4, 5, ...] (decode?)
=#
export reorder_grid_tiled, get_tiled_order

function tiled_order_2d(griddim, tiledim, invert = true)
    # griddim and tiledim are tuples of grid dimensions
    gridx = griddim[1];
    gridy = griddim[2];
    tilex = tiledim[1];
    tiley = tiledim[2];
    
    nnodes = gridx*gridy;
    tiled = zeros(Int64, nnodes); # The ordering
    
    # Count full tiles and renaining edge nodes.
    fullx = Int(floor(gridx / tilex));
    partialx = gridx - fullx*tilex;
    fully = Int(floor(gridy / tiley));
    partialy = gridy - fully*tiley;
    
    # Fill in ordering
    for ni=1:nnodes
        # Find destination tile index
        tx = Int(floor(mod(ni-1, gridx)/tilex)) + 1;
        ty = Int(floor((ni-1)/(gridx*tiley))) + 1;
    end
    
    tind = 0;
    for j = 1:(fully+1)
        ytill = j<=fully ? tiley : partialy;
        for i = 1:(fullx+1)
            xtill = i<=fullx ? tilex : partialx;
            
            # Add that tile's nodes one at a time
            for tj=1:ytill
                for ti=1:xtill
                    gind = ti + (i-1)*tilex + gridx*(tj + (j-1)*tiley - 1);
                    tind = tind + 1;
                    if invert
                        tiled[gind] = tind;
                    else
                        tiled[tind] = gind;
                    end
                end
            end
            
        end
    end
    
    return tiled;
end

function tiled_order_3d(griddim, tiledim, invert = true, tileorder=[])
    # griddim and tiledim are tuples of grid dimensions
    gridx = griddim[1];
    gridy = griddim[2];
    gridz = griddim[3];
    tilex = tiledim[1];
    tiley = tiledim[2];
    tilez = tiledim[3];
    
    nnodes = gridx*gridy*gridz;
    tiled = zeros(Int64, nnodes); # The ordering
    
    # Count full tiles and remaining edge nodes.
    fullx = Int(floor(gridx / tilex));
    partialx = gridx - fullx*tilex;
    fully = Int(floor(gridy / tiley));
    partialy = gridy - fully*tiley;
    fullz = Int(floor(gridz / tilez));
    partialz = gridz - fullz*tilez;
    
    # Fill in ordering
    if length(tileorder) == 0
        Ntiles = (fullx+1)*(fully+1)*(fullz+1);
        tileorder = 1:Ntiles;
    end
    
    tind = 0;
    for i=1:length(tileorder)
        # Get the tile coordinates in (x,y,z)
        tileind = tileorder[i];
        tx = Int(mod(tileind-1, fullx+1)) + 1;
        ty = Int(floor(mod(tileind-1, (fullx+1)*(fully+1))/(fullx+1))) + 1;
        tz = Int(floor((tileind-1)/((fullx+1)*(fully+1)))) + 1;
        xtill = tx<=fullx ? tilex : partialx;
        ytill = ty<=fully ? tiley : partialy;
        ztill = tz<=fullz ? tilez : partialz;
        
        # Add that tile's nodes one at a time
        for tk=1:ztill
            for tj=1:ytill
                for ti=1:xtill
                    gind = ti + (tx-1)*tilex + gridx*((tj + (ty-1)*tiley - 1) + gridy*(tk + (tz-1)*tilez - 1));
                    tind = tind + 1;
                    if invert
                        tiled[gind] = tind;
                    else
                        tiled[tind] = gind;
                    end
                end
            end
        end
    end
    
    return tiled;
end

function reorder_grid_tiled(grid, griddim, tiledim)
    dim = size(grid.allnodes,2);
    nel = size(grid.loc2glb,1);
    
    tiled = get_tiled_order(dim, griddim, tiledim, true);
    
    newnodes = zeros(size(grid.allnodes));
    newbdry = copy(grid.bdry);
    newloc2glb = copy(grid.loc2glb);
    newglbvertex = copy(grid.glbvertex);
    
    for mi=1:length(tiled)
        for d=1:dim
            newnodes[tiled[mi],d] = grid.allnodes[mi,d];
        end
    end
    
    for bid=1:length(newbdry)
        for i=1:length(newbdry[bid])
            newbdry[bid][i] = tiled[grid.bdry[bid][i]];
        end
    end
    
    for ei=1:nel
        for ni=1:size(newloc2glb,2)
            newloc2glb[ei,ni] = tiled[grid.loc2glb[ei,ni]];
        end
        for ni=1:size(newglbvertex,2)
            newglbvertex[ei,ni] = tiled[grid.glbvertex[ei,ni]];
        end
    end
    
    return Femshop.Grid(newnodes, newbdry, grid.bdryelem, grid.bdrynorm, grid.bids, newloc2glb, newglbvertex);
end

function get_tiled_order(dim, griddim, tiledim, invert = true)
    if dim == 3
        return tiled_order_3d(griddim, tiledim, invert);
    elseif dim == 2
        return tiled_order_2d(griddim, tiledim, invert);
    else
        return 1:griddim[1]; # 1D doesn't make sense
    end
end