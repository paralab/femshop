#=
# Uses the elemental order to determine the nodal order.
# Order within elements is:
#   interior -> empty sides(z=0,y=0,x=0,z=1,y=1,x=1) or (y=0,x=0,y=1,x=1)
=#
export reorder_grid_element_first, get_element_first_order

function element_first_order_2d(grid, porder, elorder, invert = true)
    nnodes = size(grid.allnodes,2);
    eforder = zeros(Int, nnodes); # The ordering
    node_done = zeros(Bool, nnodes);
    
    # Loop over the elements in elorder, filling in the eforder as we go
    l2g = grid.loc2glb;
    efind = 0;
    for ei=elorder
        # Add interior first
        for j = 2:porder
            for i = 2:porder
                lind = i + (porder+1)*(j-1);
                gind = l2g[lind, ei];
                # These should not be done yet.
                efind = efind+1;
                if invert
                    eforder[gind] = efind;
                else
                    eforder[efind] = gind;
                end
                node_done[gind] = true;
            end
        end
        
        # y=0 side
        for i = 1:porder+1
            lind = i;
            gind = l2g[lind, ei];
            if !node_done[gind]
                efind = efind+1;
                if invert
                    eforder[gind] = efind;
                else
                    eforder[efind] = gind;
                end
                node_done[gind] = true;
            end
        end
        
        # x=0 side
        for i = 2:porder+1
            lind = (porder+1)*(i-1) + 1;
            gind = l2g[lind, ei];
            if !node_done[gind]
                efind = efind+1;
                if invert
                    eforder[gind] = efind;
                else
                    eforder[efind] = gind;
                end
                node_done[gind] = true;
            end
        end
        
        # y=1 side
        for i = 2:porder+1
            lind = i + (porder+1)*porder;
            gind = l2g[lind, ei];
            if !node_done[gind]
                efind = efind+1;
                if invert
                    eforder[gind] = efind;
                else
                    eforder[efind] = gind;
                end
                node_done[gind] = true;
            end
        end
        
        # x=1 side
        for i = 2:porder
            lind = (porder+1)*i;
            gind = l2g[lind, ei];
            if !node_done[gind]
                efind = efind+1;
                if invert
                    eforder[gind] = efind;
                else
                    eforder[efind] = gind;
                end
                node_done[gind] = true;
            end
        end
        
    end
    
    return eforder;
end

function element_first_order_3d(grid, porder, elorder, invert = true)
    nnodes = size(grid.allnodes,2);
    eforder = zeros(Int, nnodes); # The ordering
    node_done = zeros(Bool, nnodes);
    
    # Loop over the elements in elorder, filling in the eforder as we go
    l2g = grid.loc2glb;
    efind = 0;
    for ei=elorder
        # Add interior first
        for k=2:porder
            for j = 2:porder
                for i = 2:porder
                    lind = i + (porder+1)*((j-1) + (porder+1)*(k-1));
                    gind = l2g[lind, ei];
                    # These should not be done yet.
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # z=0 side
        for j = 1:porder+1
            for i = 1:porder+1
                lind = i + (porder+1)*(j-1);
                gind = l2g[lind, ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # y=0 side
        for j = 2:porder+1
            for i = 1:porder+1
                lind = i + (porder+1)*(porder+1)*(j-1);
                gind = l2g[lind, ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # x=0 side
        for j = 2:porder+1
            for i = 2:porder+1
                lind = (porder+1)*(i-1) + 1 + (porder+1)*(porder+1)*(j-1);
                gind = l2g[lind,ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # z=1 side
        for j = 2:porder+1
            for i = 2:porder+1
                lind = i + (porder+1)*(j-1) + (porder+1)*(porder+1)*porder;
                gind = l2g[lind, ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # y=1 side
        for j = 2:porder
            for i = 2:porder+1
                lind = i + (porder+1)*(porder+1)*(j-1) + (porder+1)*porder;
                gind = l2g[lind, ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
        # x=1 side
        for j = 2:porder
            for i = 2:porder
                lind = (porder+1)*(i-1) + 1 + (porder+1)*(porder+1)*(j-1) + porder;
                gind = l2g[lind, ei];
                if !node_done[gind]
                    efind = efind+1;
                    if invert
                        eforder[gind] = efind;
                    else
                        eforder[efind] = gind;
                    end
                    node_done[gind] = true;
                end
            end
        end
        
    end
    
    return eforder;
end

function reorder_grid_element_first(grid, porder, elorder)
    dim = size(grid.allnodes,1);
    nel = size(grid.loc2glb,2);
    nfc = size(grid.face2glb,2);
    
    norder= get_element_first_order(dim, grid, porder, elorder, true);
    
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

function get_element_first_order(dim, grid, porder, elorder, invert = true)
    if dim == 3
        return element_first_order_3d(grid, porder, elorder, invert);
    elseif dim == 2
        return eforder_order_2d(grid, porder, elorder, invert);
    else
        return 1:size(grid.allnodes,2); # 1D doesn't make sense
    end
end