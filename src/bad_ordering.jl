#=
# A random node reordering.
=#
export reorder_grid_random, random_order

function random_order(N, seed)
    rng = Random.MersenneTwister(seed); # Use a specific RNG seed to make results consistent
    return randperm(rng, N); # Randomly permute node order
end

function reorder_grid_random(grid, seed = 17)
    dim = size(grid.allnodes,1);
    nnodes = size(grid.allnodes,2);
    nel = size(grid.loc2glb,2);
    nfc = size(grid.face2glb,2);
    
    norder= random_order(nnodes, seed);
    
    newnodes = zeros(size(grid.allnodes));
    newbdry = copy(grid.bdry);
    newloc2glb = copy(grid.loc2glb);
    newglbvertex = copy(grid.glbvertex);
    newface2glb = copy(grid.face2glb);
    newfaceVertex2glb = copy(grid.faceVertex2glb);
    
    for mi=1:length(norder)
        for d=1:dim
            newnodes[d,norder[mi]] = grid.allnodes[d,mi];
        end
    end
    
    for bid=1:length(newbdry)
        for i=1:length(newbdry[bid])
            newbdry[bid][i] = norder[grid.bdry[bid][i]];
        end
    end
    
    for ei=1:nel
        for ni=1:size(newloc2glb,1)
            newloc2glb[ni,ei] = norder[grid.loc2glb[ni,ei]];
        end
        for ni=1:size(newglbvertex,1)
            newglbvertex[ni,ei] = norder[grid.glbvertex[ni,ei]];
        end
    end
    
    for fi=1:nfc
        for ni=1:size(newface2glb,1)
            newface2glb[ni,fi] = norder[grid.face2glb[ni,fi]];
        end
        for ni=1:size(newfaceVertex2glb,1)
            newfaceVertex2glb[ni,fi] = norder[grid.faceVertex2glb[ni,fi]];
        end
    end
    
    return Femshop.Grid(newnodes, newbdry, grid.bdryface, grid.bdrynorm, grid.bids, newloc2glb, newglbvertex, newface2glb, newfaceVertex2glb);
end