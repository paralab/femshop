# cachsim output
import ..Femshop: init_cachesimout, add_cachesim_array, cachesim_load, cachesim_store

function linear_solve_cachesim(var, bilinear, linear, stepper=nothing)
    # if config.linalg_matrixfree
    #     return solve_matrix_free_sym(var, bilinear, linear, stepper);
    #     #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    # end
    N1 = size(grid_data.allnodes)[1];
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
    end
    init_cachesimout(N1, refel, dofs_per_node);
    
    if prob.time_dependent && !(stepper === nothing)
        #TODO time dependent coefficients
        assemble_t = @elapsed((A, b) = assemble_cachesim(var, bilinear, linear, 0, stepper.dt));
        log_entry("Assembly took "*string(assemble_t)*" seconds");

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = [];
        start_t = Base.Libc.time();
        for i=1:stepper.Nsteps
            b = assemble_rhs_only_cachesim(var, linear, t, stepper.dt);
            sol = A\b;
            # place the values in the variable value arrays
            if typeof(var) <: Array
                tmp = 0;
                totalcomponents = 0;
                for vi=1:length(var)
                    totalcomponents = totalcomponents + length(var[vi].symvar.vals);
                end
                for vi=1:length(var)
                    components = length(var[vi].symvar.vals);
                    for compi=1:components
                        var[vi].values[:,compi] = sol[(compi+tmp):totalcomponents:end];
                        tmp = tmp + 1;
                    end
                end
            else
                components = length(var.symvar.vals);
                for compi=1:components
                    var.values[:,compi] = sol[compi:components:end];
                end
            end

            t += stepper.dt;
        end
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        #display(sol);
		# outfile = "linear_sol.txt"
		# open(outfile, "w") do f
  		# 	for ii in sol
    	# 		println(f, ii)
  		# 	end
		# end # the file f is automatically closed after this block finishes
        return sol;

    else
        assemble_t = @elapsed((A, b) = assemble_cachesim(var, bilinear, linear));
        sol_t = @elapsed(sol = A\b);

        log_entry("Assembly took "*string(assemble_t)*" seconds");
        log_entry("Linear solve took "*string(sol_t)*" seconds");
        #display(A);
        #display(b);
        #display(sol);
        return sol;
    end
end

function assemble_cachesim(var, bilinear, linear, t=0.0, dt=0.0)
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes)[1];
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
    end
    Nn = dofs_per_node * N1;
    
    b = zeros(Nn);
    A = spzeros(Nn, Nn);
    
    allnodes_id = add_cachesim_array(size(grid_data.allnodes),8);
    
    # The elemental assembly loop
    for e=1:nel
        nv = mesh_data.nv[e];
        gis = zeros(Int, nv);
        for vi=1:nv
            gis[vi] = mesh_data.invind[mesh_data.elements[e,vi]]; # mesh indices of element's vertices
        end

        vx = mesh_data.nodes[gis,:];        # coordinates of element's vertices
        glb = grid_data.loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions
        cachesim_load_range(allnodes_id, 1:refel.dim, glb[:]);
        
        
        # The linear part. Compute the elemental linear part for each dof
        rhsargs = (var, xe, glb, refel, RHS, t, dt);
        lhsargs = (var, xe, glb, refel, LHS, t, dt);
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;
            cachesim_load_range(2, glb);
            cachesim_load_range(5);
            cachesim_store_range(2, glb);

            bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
            A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
            cachesim_load_range(1, glb, glb);
            cachesim_load_range(4);
            cachesim_store_range(1, glb, glb);
            
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear_cachesim!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear_cachesim!(A, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        else
            linchunk = linear.func(rhsargs);
            insert_linear_cachesim!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear_cachesim!(A, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        end
    end
    
    # Boundary conditions
    bidcount = length(grid_data.bids); # the number of BIDs
    if dofs_per_node > 1
        if multivar
            dofind = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar.vals)
                    dofind = dofind + 1;
                    for bid=1:bidcount
                        if prob.bc_type[var[vi].index, bid] == DIRICHLET
                            (A, b) = dirichlet_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
                        elseif prob.bc_type[var[vi].index, bid] == NEUMANN
                            (A, b) = neumann_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                        else
                            printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, bid]);
                        end
                    end
                end
            end
        else
            for d=1:dofs_per_node
                dofind = d;
                for bid=1:bidcount
                    if prob.bc_type[var.index, bid] == DIRICHLET
                        (A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, dofind, dofs_per_node);
                    elseif prob.bc_type[var.index, bid] == NEUMANN
                        (A, b) = neumann_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                    else
                        printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, bid]);
                    end
                end
            end
        end
    else
        for bid=1:bidcount
            if prob.bc_type[var.index, bid] == DIRICHLET
                (A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid], grid_data.bdry[bid], t);
            elseif prob.bc_type[var.index, bid] == NEUMANN
                (A, b) = neumann_bc(A, b, prob.bc_func[var.index, bid], grid_data.bdry[bid], bid, t);
            else
                printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, bid]);
            end
        end
    end

    return (A, b);
end

# assembles the A and b in Au=b
function assemble_rhs_only_cachesim(var, linear, t=0.0, dt=0.0)
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes)[1];
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
    end
    Nn = dofs_per_node * N1;

    b = zeros(Nn);

    for e=1:nel;
        glb = grid_data.loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions

        rhsargs = (var, xe, glb, refel, RHS, t, dt);

        #linchunk = linear.func(args);  # get the elemental linear part
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;

        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear_cachesim!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

        else
            linchunk = linear.func(rhsargs);
            insert_linear_cachesim!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

        end
    end

    # Just for testing. This should really be a loop over BIDs with corresponding calls
    #b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, 1], grid_data.bdry[1,:], t);
    # Boundary conditions
    bidcount = length(grid_data.bids); # the number of BIDs
    if dofs_per_node > 1
        if multivar
            dofind = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar.vals)
                    dofind = dofind + 1;
                    for bid=1:bidcount
                        b = dirichlet_bc_rhs_only(b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
                    end
                end
            end
        else
            for d=1:dofs_per_node
                #rows = ((d-1)*length(glb)+1):(d*length(glb));
                dofind = d;
                for bid=1:bidcount
                    b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, d, dofs_per_node);
                end
            end
        end
    else
        for bid=1:bidcount
            b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid], grid_data.bdry[bid], t);
        end
    end

    return b;
end

# Inset the single dof into the greater construct
function insert_linear_cachesim!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));

        b[ind] = b[ind] + bel[ind2];
        
        cachesim_load_range(2,ind);
        cachesim_load_range(5,ind2);
        cachesim_store_range(2,ind);
    end
end

function insert_bilinear_cachesim!(a, ael, glb, dof, Ndofs)
    # group nodal dofs
    for di=1:length(dof)
        indi = glb.*Ndofs .- (Ndofs-dof[di]);
        indi2 = ((di-1)*length(glb)+1):(di*length(glb));
        for dj=1:length(dof)
            indj = glb.*Ndofs .- (Ndofs-dof[dj]);
            indj2 = ((dj-1)*length(glb)+1):(dj*length(glb));

            a[indi, indj] = a[indi, indj] + ael[indi2, indj2];
            cachesim_load_range(1,indi, indj);
            cachesim_load_range(4,indi2, indj2);
            cachesim_store_range(1,indi, indj);
        end
    end
end