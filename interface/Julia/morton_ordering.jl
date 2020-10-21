# Builds a 3D morton ordering of the nodes
# Assumes the grid is (2^N)^3
# Also assumes first node and last node are at opposing extremes to compute bounding box if not given.
# Returns the order in which nodes shall be indexed.
export reorder_grid_morton

function morton_order_nodes_3d(nodes, N, invert = false)
    nnodes = (2^N)^3;
    if !(size(nodes,1) == nnodes)
        println("Number of nodes must be like (2^N)^3 for morton");
        return 1:size(nodes,1);
    end
    
    # Since node locations could be non-uniform, assign each an integer location
    lexorder = zeros(Int,size(nodes));
    twoN = 2^N;
    for k=1:twoN
        for j=1:twoN
            for i=1:twoN
                ni = i + twoN*((j-1) + twoN*(k-1));
                lexorder[ni,1] = i;
                lexorder[ni,2] = j;
                lexorder[ni,3] = k;
            end
        end
    end
    
    bbox = zeros(6);
    bbox[1] = 1;
    bbox[3] = 1;
    bbox[5] = 1;
    bbox[2] = twoN;
    bbox[4] = twoN;
    bbox[6] = twoN;
    
    center = [(bbox[1]+bbox[2])/2, (bbox[3]+bbox[4])/2, (bbox[5]+bbox[6])/2];
    tmpbbox = zeros(6);
    tmpcenter = zeros(3);
    
    morton = zeros(Int64, nnodes);
    for ni=1:nnodes
        vertex = lexorder[ni,:];
        for bi=1:6
            tmpbbox[bi] = bbox[bi];
        end
        tmpcenter[1] = center[1];
        tmpcenter[2] = center[2];
        tmpcenter[3] = center[3];
        
        # Find the morton index for this element
        index = Int64(0);
        for i=1:N
            stupid = index;
            index = index<<3;
            if vertex[1] > tmpcenter[1]
                index = index+1;
            end
            if vertex[2] > tmpcenter[2]
                index = index+2;
            end
            if vertex[3] > tmpcenter[3]
                index = index+4;
            end
            
            for j=1:3
                if vertex[j] > tmpcenter[j]
                    tmpbbox[2*j-1] = tmpcenter[j];
                else
                    tmpbbox[2*j] = tmpcenter[j];
                end
            end
            tmpcenter = [(tmpbbox[1]+tmpbbox[2])/2, (tmpbbox[3]+tmpbbox[4])/2, (tmpbbox[5]+tmpbbox[6])/2];
        end
        
        # add one for 1-based index
        index += 1;
        if invert
            morton[index] = ni;
        else
            morton[ni] = index;
        end
    end
    
    return morton;
end

function reorder_grid_morton(grid, N)
    dim = size(grid.allnodes,2);
    nel = size(grid.loc2glb,1);
    morton = morton_order_nodes_3d(grid.allnodes, N, true);
    
    newnodes = zeros(size(grid.allnodes));
    newbdry = copy(grid.bdry);
    newloc2glb = copy(grid.loc2glb);
    newglbvertex = copy(grid.glbvertex);
    
    for mi=1:length(morton)
        for d=1:dim
            newnodes[morton[mi],d] = grid.allnodes[mi,d];
        end
    end
    
    for bid=1:length(newbdry)
        for i=1:length(newbdry[bid])
            newbdry[bid][i] = morton[grid.bdry[bid][i]];
        end
    end
    
    for ei=1:nel
        for ni=1:size(newloc2glb,2)
            newloc2glb[ei,ni] = morton[grid.loc2glb[ei,ni]];
        end
        for ni=1:size(newglbvertex,2)
            newglbvertex[ei,ni] = morton[grid.glbvertex[ei,ni]];
        end
    end
    
    return Femshop.Grid(newnodes, newbdry, grid.bdryelem, grid.bdrynorm, grid.bids, newloc2glb, newglbvertex);
end