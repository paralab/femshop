#=
# Format specific functions for writing data output files.
=#

function output_values_raw(vars, file)
    # do one variable/component at a time
    if typeof(vars) <: Array 
        for vi=1:length(vars) # each variable
            output_values_raw(vars[vi], file);
        end
    else
        # Here vars is just one variable
        N = size(vars.values, 2); # number of points
        for ci=1:size(vars.values,1) # each component for that variable
            file.write(vars.values[ci,:]); # write the values
        end
    end
    
end

function output_values_csv(vars, file)
    # columns will be like x, y, z, u1,...
    if typeof(vars) <: Array
        N = size(vars[1].values, 2); # number of points (must be same for all vars)
        if vars[1].location == CELL
            x = fv_info.cellCenters;
        else
            x = grid_data.allnodes;
        end
        
        # header
        line = "x";
        if dim > 1
            line *= ", y";
        end
        if dim > 2
            line *= ", z";
        end
        for vi=1:length(vars) # each variable
            if size(vars[vi].values,1) > 1
                for ci=1:size(vars[vi].values,1)
                    line *= ", "*string(vars[vi].symbol)*"_"*string(ci);
                end
            else
                line *= ", "*string(vars[vi].symbol);
            end
        end
        println(file, line);
        
        # data
        for i=1:N
            # x, y, z
            line = string(x[1,i]);
            if dim > 1
                line *= ", "*string(x[2,i]);
            end
            if dim > 2
                line *= ", "*string(x[3,i]);
            end
            
            # u1, u2, ...
            for vi=1:length(vars) # each variable
                for ci=1:size(vars[vi].values,1)
                    line *= ", "*string(vars[vi].values[ci,i]);
                end
            end
            println(file, line);
        end
        
    else
        N = size(vars.values, 2); # number of points
        if vars.location == CELL
            x = fv_info.cellCenters;
        else
            x = grid_data.allnodes;
        end
        dim = size(x,1);
        
        # header
        line = "x";
        if dim > 1
            line *= ", y";
        end
        if dim > 2
            line *= ", z";
        end
        if size(vars.values,1) > 1
            for ci=1:size(vars.values,1)
                line *= ", "*string(vars.symbol)*"_"*string(ci);
            end
        else
            line *= ", "*string(vars.symbol);
        end
        println(file, line);
        
        # data
        for i=1:N
            # x, y, z
            line = string(x[1,i]);
            if dim > 1
                line *= ", "*string(x[2,i]);
            end
            if dim > 2
                line *= ", "*string(x[3,i]);
            end
            
            # u1, u2, ...
            for ci=1:size(vars.values,1)
                line *= ", "*string(vars.values[ci,i]);
            end
            println(file, line);
        end
    end
end

function output_values_vtk(vars, file, format="BINARY")
    # This is not ready
    printerr("vtk output format not ready")
    
    # N = size(vars[1].values, 2); # number of points
    # ncells = size(grid_data.loc2glb, 2); # number of elements
    
    # file.println("# vtk DataFile Version 3.0");
    # file.println("Output from Femshop");
    # file.println(format);
    # file.println("DATASET STRUCTURED_GRID");
    # s
    # file.println("POINTS "*string(N)*" double");
    
    # file.println("CELLS "*string(ncells)*" double");
end
