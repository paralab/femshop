#=
# CG solver
=#
module CGSolver

export init_cgsolver, solve, nonlinear_solve

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
            NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
            NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4,
            ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, NEUMANN, ROBIN,
            MSH_V2, MSH_V4,
            SCALAR, VECTOR, TENSOR, SYM_TENSOR,
            LHS, RHS,
            LINEMESH, QUADMESH, HEXMESH
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data, grid_data, refel, time_stepper, elemental_order
import ..Femshop: Variable, Coefficient, GenFunction
import ..Femshop: geometric_factors

using LinearAlgebra, SparseArrays

include("cg_boundary.jl");
include("nonlinear.jl")
include("cg_matrixfree.jl");
include("cachsim_solve.jl");

function init_cgsolver()
    dim = config.dimension;

    # build initial conditions
    for vind=1:length(variables)
        if vind <= length(prob.initial)
            if prob.initial[vind] != nothing
                variables[vind].values = zeros(size(variables[vind].values));
                # Set initial condition
                if typeof(prob.initial[vind]) <: Array
                    for ci=1:length(prob.initial[vind])
                        for ni=1:size(grid_data.allnodes)[1]
                            if dim == 1
                                variables[vind].values[ni,ci] = prob.initial[vind][ci].func(grid_data.allnodes[ni],0,0,0);
                            elseif dim == 2
                                variables[vind].values[ni,ci] = prob.initial[vind][ci].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],0,0);
                            elseif dim == 3
                                variables[vind].values[ni,ci] = prob.initial[vind][ci].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],grid_data.allnodes[ni,3],0);
                            end
                        end
                    end
                else
                    for ni=1:size(grid_data.allnodes)[1]
                        if dim == 1
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni],0,0,0);
                        elseif dim == 2
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],0,0);
                        elseif dim == 3
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],grid_data.allnodes[ni,3],0);
                        end
                    end
                end

                variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(variables[vind].symbol));
            end
        end
    end
end

function linear_solve(var, bilinear, linear, stepper=nothing)
    if config.linalg_matrixfree
        return solve_matrix_free_sym(var, bilinear, linear, stepper);
        #return solve_matrix_free_asym(var, bilinear, linear, stepper);
    end
    if prob.time_dependent && !(stepper === nothing)
        #TODO time dependent coefficients
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, 0, stepper.dt));
        log_entry("Assembly took "*string(assemble_t)*" seconds");

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = [];
        start_t = Base.Libc.time();
        for i=1:stepper.Nsteps
            b = assemble_rhs_only(var, linear, t, stepper.dt);
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
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear));
        sol_t = @elapsed(sol = A\b);

        log_entry("Assembly took "*string(assemble_t)*" seconds");
        log_entry("Linear solve took "*string(sol_t)*" seconds");
        #display(A);
        #display(b);
        #display(sol);
        return sol;
    end
end

function nonlinear_solve(var, nlvar, bilinear, linear, stepper=nothing)
    if prob.time_dependent && !(stepper === nothing)
        #TODO time dependent coefficients
        
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        start_t = Base.Libc.time();
        nl = nonlinear(100, 1e-12, 1e-12);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        for i=1:stepper.Nsteps
			newton(nl, assemble, assemble_rhs_only, nlvar, t, stepper.dt);
            t += stepper.dt;
			nlvar[2] = nlvar[1];
			nlvar[2] = nlvar[1];
        end
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        #display(sol);
		outfile = "nonlinear_sol.txt"
		open(outfile, "w") do f
  			for ii in nlvar.values
    			println(f, ii)
  			end
		end # the file f is automatically closed after this block finishes
        #return nlvar.values;
        return [];

    else
        start_t = Base.Libc.time();
        nl = nonlinear(100, 1e-12, 1e-12);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        newton(nl, assemble, assemble_rhs_only, nlvar);
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        
		outfile = "nonlinear_sol.txt"
		open(outfile, "w") do f
  			for ii in nlvar.values
    			println(f, ii)
  			end
		end # the file f is automatically closed after this block finishes
        #return nlvar.values;
        return [];
    end
end

# assembles the A and b in Au=b
function assemble(var, bilinear, linear, t=0.0, dt=0.0)
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
    #A = spzeros(Nn, Nn);
    AI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    AJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
    AV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
    
    loop_time = Base.Libc.time();
    # Elemental loop follows elemental ordering
    for e=elemental_order;
        gis = grid_data.glbvertex[e,:];
        vx = grid_data.allnodes[gis,:];         # coordinates of element's vertices
        glb = grid_data.loc2glb[e,:];           # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];      # coordinates of this element's nodes for evaluating coefficient functions
        
        Astart = (e-1)*Np*dofs_per_node*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element

        # The linear part. Compute the elemental linear part for each dof
        rhsargs = (var, xe, glb, refel, RHS, t, dt);
        lhsargs = (var, xe, glb, refel, LHS, t, dt);
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;

            bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
            #A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
            for jj=1:Np
                offset = Astart - 1 + (jj-1)*Np;
                for ii=1:Np
                    AI[offset + ii] = glb[ii];
                    AJ[offset + ii] = glb[jj];
                    AV[offset + ii] = bilinchunk[ii, jj];
                end
            end
            
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        else
			
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        end
    end
    loop_time = Base.Libc.time() - loop_time;
    
    # Build the sparse A. Uses default + to combine overlaps
    sparse_time = Base.Libc.time();
    A = sparse(AI, AJ, AV);
    sparse_time = Base.Libc.time() - sparse_time;
    
    # Boundary conditions
    bc_time = Base.Libc.time();
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
    bc_time = Base.Libc.time() - bc_time;
    
    log_entry("Elemental loop time:     "*string(loop_time));
    log_entry("Form sparse matrix time: "*string(sparse_time));
    log_entry("Boundary condition time: "*string(bc_time));
    
    return (A, b);
end

# assembles the A and b in Au=b
function assemble_rhs_only(var, linear, t=0.0, dt=0.0)
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
    
    # Elemental loop follows elemental ordering
    for e=elemental_order;
        glb = grid_data.loc2glb[e,:];       # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions

        rhsargs = (var, xe, glb, refel, RHS, t, dt);

        #linchunk = linear.func(args);  # get the elemental linear part
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;

        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

        else
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

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
function insert_linear!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));

        b[ind] = b[ind] + bel[ind2];
    end
end

function insert_bilinear!(AI, AJ, AV, Astart, ael, glb, dof, Ndofs)
    Np = length(glb);
    # group nodal dofs
    for dj=1:length(dof)
        indj = glb.*Ndofs .- (Ndofs-dof[dj]);
        indj2 = ((dj-1)*Np+1):(dj*Np);
        for di=1:length(dof)
            indi = glb.*Ndofs .- (Ndofs-dof[di]);
            indi2 = ((di-1)*Np+1):(di*Np);
            
            #a[indi, indj] = a[indi, indj] + ael[indi2, indj2];
            
            for jj=1:Np
                offset = Astart + (jj-1 + Np*(dj-1))*Np*Ndofs + Np*(di-1) - 1;
                for ii=1:Np
                    AI[offset + ii] = indi[ii];
                    AJ[offset + ii] = indj[jj];
                    AV[offset + ii] = ael[indi2[ii], indj2[jj]];
                end
            end
        end
    end
    
end

end #module
