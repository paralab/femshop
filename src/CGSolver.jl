#=
# CG solver
=#
module CGSolver

export init_cgsolver, solve, nonlinear_solve

import ..Femshop: JULIA, CPP, MATLAB, DENDRO, HOMG, CUSTOM_GEN_TARGET,
        SQUARE, IRREGULAR, UNIFORM_GRID, TREE, UNSTRUCTURED, 
        CG, DG, HDG, FV,
        NODAL, MODAL, CELL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, 
        NONLINEAR_NEWTON, NONLINEAR_SOMETHING, 
        EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4, ABM4, 
        DEFAULT_SOLVER, PETSC, 
        VTK, RAW_OUTPUT, CUSTOM_OUTPUT, 
        DIRICHLET, NEUMANN, ROBIN, NO_BC, FLUX,
        MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR, SYM_TENSOR,
        LHS, RHS,
        LINEMESH, QUADMESH, HEXMESH
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data, grid_data, refel, time_stepper, elemental_order
import ..Femshop: Variable, Coefficient, GenFunction
import ..Femshop: geometric_factors, build_deriv_matrix

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
                        for ni=1:size(grid_data.allnodes,2)
                            if dim == 1
                                variables[vind].values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[ni],0,0,0);
                            elseif dim == 2
                                variables[vind].values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],0,0);
                            elseif dim == 3
                                variables[vind].values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],grid_data.allnodes[3,ni],0);
                            end
                        end
                    end
                else
                    for ni=1:size(grid_data.allnodes,2)
                        if dim == 1
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni],0,0,0);
                        elseif dim == 2
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],0,0);
                        elseif dim == 3
                            variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],grid_data.allnodes[3,ni],0);
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
        # Shall we save geometric factors to reuse each time step?
        save_geometric_factors = true;
        
        if save_geometric_factors
            assemble_t = @elapsed((A, b, wdetJ, J) = assemble(var, bilinear, linear, 0, stepper.dt; keep_geometric_factors=true));
            saved_GF = (wdetJ, J);
        else
            assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, 0, stepper.dt));
            saved_GF = nothing;
        end
        
        log_entry("Initial assembly took "*string(assemble_t)*" seconds");

        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = [];
        start_t = Base.Libc.time();
        assemble_t = 0;
        linsolve_t = 0;
        last2update = 0;
        last10update = 0;
        print("Time stepping progress(%): 0");
        for i=1:stepper.Nsteps
            if stepper.stages > 1
                # LSRK4 is a special case, low storage
                if stepper.type == LSRK4
                    # Low storage RK4: 
                    # p0 = u
                    #   ki = ai*k(i-1) + dt*f(p(i-1), t+ci*dt)
                    #   pi = p(i-1) + bi*ki
                    # u = p5
                    
                    tmppi = get_var_vals(var);
                    tmpki = zeros(size(b));
                    for rki=1:stepper.stages
                        rktime = t + stepper.c[rki]*stepper.dt;
                        # p(i-1) is currently in u
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, rktime, stepper.dt; rhs_only = true, saved_geometric_factors = saved_GF));
                        
                        linsolve_t += @elapsed(sol = A\b);
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        # Need to apply boundary conditions now to tmppi
                        tmppi = apply_boundary_conditions_rhs_only(var, tmppi, rktime);
                        
                        place_sol_in_vars(var, tmppi, stepper);
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    result = get_var_vals(var);
                    # will be placed in var.values for each stage
                    tmpvals = result;
                    # Storage for each stage
                    tmpki = zeros(length(b), stepper.stages);
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        assemble_t += @elapsed(b = assemble(var, nothing, linear, stime, stepper.dt; rhs_only = true, saved_geometric_factors = saved_GF));
                        
                        linsolve_t += @elapsed(tmpki[:,stage] = A\b);
                        
                        tmpvals = result;
                        for j=1:(stage-1)
                            if stepper.a[stage, j] > 0
                                tmpvals += stepper.dt * stepper.a[stage, j] .* tmpki[:,j];
                            end
                        end
                        
                        # Need to apply boundary conditions now
                        tmpvals = apply_boundary_conditions_rhs_only(var, tmpvals, stime);
                        
                        place_sol_in_vars(var, tmpvals, stepper);
                    end
                    for stage=1:stepper.stages
                        result += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    result = apply_boundary_conditions_rhs_only(var, result, t + stepper.dt);
                    place_sol_in_vars(var, result, stepper);
                end
                
            else
                assemble_t += @elapsed(b = assemble(var, nothing, linear, t, stepper.dt; rhs_only = true, saved_geometric_factors = saved_GF));
                
                linsolve_t += @elapsed(sol = A\b);
                
                place_sol_in_vars(var, sol, stepper);
            end
            
            t += stepper.dt;
            
            progressPercent = Int(floor(i*100.0/stepper.Nsteps));
            if progressPercent - last2update >= 2
                last2update = progressPercent;
                if progressPercent - last10update >= 10
                    print(progressPercent);
                    last10update = progressPercent;
                else
                    print(".");
                end
            end
        end
        println("");
        end_t = Base.Libc.time();

        log_entry("Stepping took "*string(end_t-start_t)*" seconds. ("*string(assemble_t)*" for assembly, "*string(linsolve_t)*" for linear solve)");
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
        # uncomment to look at A
        #global Amat = A;
        
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
        nl = nonlinear(100, 1e-9, 1e-9);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        for i=1:stepper.Nsteps
			println("solve for time step ", i);
			newton(nl, assemble, assemble_rhs_only, nlvar, t, stepper.dt);
            t += stepper.dt;
			nlvar[4].values = copy(nlvar[1].values);
			nlvar[5].values = copy(nlvar[2].values);
			if (i % 10 ==0)
				outfile = string("LD_u_",i);
				open(outfile, "w") do f
					for ii=1:length(nlvar[1].values)
						println(f, nlvar[1].values[ii])
					end
					close(f)
				end
				outfile = string("LD_v_",i);
				open(outfile, "w") do f
					for ii=1:length(nlvar[2].values)
						println(f, nlvar[2].values[ii])
					end
					close(f)
				end
			end
        end
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        #display(sol);
        #return nlvar.values;
        return [];

    else
        start_t = Base.Libc.time();
        nl = nonlinear(100, 1e-12, 1e-12);
        init_nonlinear(nl, var, nlvar, bilinear, linear);
        newton(nl, assemble, assemble_rhs_only, nlvar);
        end_t = Base.Libc.time();

        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        
        #return nlvar.values;
        return [];
    end
end

# assembles the A and b in Au=b
function assemble(var, bilinear, linear, t=0.0, dt=0.0; rhs_only = false, keep_geometric_factors = false, saved_geometric_factors = nothing)
    Np = refel.Np;
    nel = mesh_data.nel;
    N1 = size(grid_data.allnodes,2);
    multivar = typeof(var) <: Array;
    maxvarindex = 0;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
            maxvarindex = max(maxvarindex,var[vi].index);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
        maxvarindex = var.index;
    end
    Nn = dofs_per_node * N1;

    b = zeros(Nn);
    
    if !rhs_only
        AI = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
        AJ = zeros(Int, nel*dofs_per_node*Np*dofs_per_node*Np);
        AV = zeros(nel*dofs_per_node*Np*dofs_per_node*Np);
    end
    
    # Save geometric factors if desired
    if keep_geometric_factors
        saved_wdetj = [];
        saved_J = [];
    end
    
    # Stiffness and mass are precomputed for uniform grid meshes
    precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
    if precomputed_mass_stiffness
        glb = grid_data.loc2glb[:,1];
        xe = grid_data.allnodes[:,glb[:]];
        if saved_geometric_factors === nothing
            (detJ, J) = geometric_factors(refel, xe);
            wdetj = refel.wg .* detJ;
            if keep_geometric_factors
                saved_wdetj = [wdetj];
                saved_J = [J];
            end
        else
            wdetj = saved_geometric_factors[1][1];
            J = saved_geometric_factors[2][1];
        end
        if config.dimension == 1
            (RQ1, RD1) = build_deriv_matrix(refel, J);
            TRQ1 = RQ1';
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1)];
        elseif config.dimension == 2
            (RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2) = (RQ1', RQ2');
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2)];
        else
            (RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);
            (TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');
            stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2) , (TRQ3 * diagm(wdetj) * RQ3)];
        end
        mass = (refel.Q)' * diagm(wdetj) * refel.Q;
    else
        stiffness = 0;
        mass = 0;
    end
    
    loop_time = Base.Libc.time();
    # Elemental loop follows elemental ordering
    for ei=1:nel
        e = elemental_order[ei];
        glb = grid_data.loc2glb[:,e];           # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[:,glb[:]];      # coordinates of this element's nodes for evaluating coefficient functions
        
        if !precomputed_mass_stiffness
            if saved_geometric_factors === nothing
                (detJ, J) = geometric_factors(refel, xe);
                wdetj = refel.wg .* detJ;
                
                if keep_geometric_factors
                    push!(saved_wdetj, wdetj);
                    push!(saved_J, J);
                end
            else
                wdetj = saved_geometric_factors[1][e];
                J = saved_geometric_factors[2][e];
            end
        end
        
        rhsargs = (var, xe, glb, refel, wdetj, J, t, dt, stiffness, mass);
        
        if !rhs_only
            Astart = (e-1)*Np*dofs_per_node*Np*dofs_per_node + 1; # The segment of AI, AJ, AV for this element
            lhsargs = (var, xe, glb, refel, wdetj, J, t, dt, stiffness, mass);
        end
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;
            
            if !rhs_only
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
            end
            
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);
            
            if !rhs_only
                bilinchunk = bilinear.func(lhsargs);
                insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
            end
        else
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);
            
            if !rhs_only
                bilinchunk = bilinear.func(lhsargs);
                insert_bilinear!(AI, AJ, AV, Astart, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
            end
        end
    end
    loop_time = Base.Libc.time() - loop_time;
    
    if !rhs_only
        # Build the sparse A. Uses default + to combine overlaps
        A = sparse(AI, AJ, AV);
    end
    
    # Boundary conditions
    bc_time = Base.Libc.time();
    bidcount = length(grid_data.bids); # the number of BIDs
    if !rhs_only
        dirichlet_rows = zeros(0);
        neumann_rows = zeros(0);
        neumann_Is = zeros(Int,0);
        neumann_Js = zeros(Int,0);
        neumann_Vs = zeros(0);
    end
    if dofs_per_node > 1
        if multivar
            dofind = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar.vals)
                    dofind = dofind + 1;
                    for bid=1:bidcount
                        if prob.bc_type[var[vi].index, bid] == NO_BC
                            # do nothing
                        elseif rhs_only
                            # When doing RHS only, values are simply placed in the vector.
                            b = dirichlet_bc_rhs_only(b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
                        elseif prob.bc_type[var[vi].index, bid] == DIRICHLET
                            #(A, b) = dirichlet_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
                            (tmprows, b) = dirichlet_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
                            append!(dirichlet_rows, tmprows);
                        elseif prob.bc_type[var[vi].index, bid] == NEUMANN
                            #(A, b) = neumann_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                            (tmprows, tmpIs, tmpJs, tmpVs, b) = neumann_bc(A, b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                            append!(neumann_rows, tmprows);
                            append!(neumann_Is, tmpIs);
                            append!(neumann_Js, tmpJs);
                            append!(neumann_Vs, tmpVs);
                        elseif prob.bc_type[var[vi].index, bid] == ROBIN
                            printerr("Robin BCs not ready.");
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
                    if prob.bc_type[var.index, bid] == NO_BC
                        # do nothing
                    elseif rhs_only
                        # When doing RHS only, values are simply placed in the vector.
                        b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, d, dofs_per_node);
                    elseif prob.bc_type[var.index, bid] == DIRICHLET
                        #(A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, dofind, dofs_per_node);
                        (tmprows, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, dofind, dofs_per_node);
                        append!(dirichlet_rows, tmprows);
                    elseif prob.bc_type[var.index, bid] == NEUMANN
                        #(A, b) = neumann_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                        (tmprows, tmpIs, tmpJs, tmpVs, b) = neumann_bc(A, b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], bid, t, dofind, dofs_per_node);
                        append!(neumann_rows, tmprows);
                        append!(neumann_Is, tmpIs);
                        append!(neumann_Js, tmpJs);
                        append!(neumann_Vs, tmpVs);
                    elseif prob.bc_type[var.index, bid] == ROBIN
                        printerr("Robin BCs not ready.");
                    else
                        printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, bid]);
                    end
                end
            end
        end
    else
        for bid=1:bidcount
            if prob.bc_type[var.index, bid] == NO_BC
                # do nothing
            elseif rhs_only
                # When doing RHS only, values are simply placed in the vector.
                b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], t);
            elseif prob.bc_type[var.index, bid] == DIRICHLET
                #(A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], t);
                (tmprows, b) = dirichlet_bc(A, b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], t);
                append!(dirichlet_rows, tmprows);
            elseif prob.bc_type[var.index, bid] == NEUMANN
                #(A, b) = neumann_bc(A, b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], bid, t);
                (tmprows, tmpIs, tmpJs, tmpVs, b) = neumann_bc(A, b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], bid, t);
                append!(neumann_rows, tmprows);
                append!(neumann_Is, tmpIs);
                append!(neumann_Js, tmpJs);
                append!(neumann_Vs, tmpVs);
            elseif prob.bc_type[var.index, bid] == ROBIN
                printerr("Robin BCs not ready.");
            else
                printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, bid]);
            end
        end
    end
    
    if !rhs_only
        if length(dirichlet_rows)>0
            A = identity_rows(A, dirichlet_rows, length(b));
        end
        if length(neumann_rows)>0
            A = insert_sparse_rows(A, neumann_Is, neumann_Js, neumann_Vs);
        end
    end
    
    # Reference points
    if size(prob.ref_point,1) >= maxvarindex
        if multivar
            posind = zeros(Int,0);
            vals = zeros(0);
            for vi=1:length(var)
                if prob.ref_point[var[vi].index,1]
                    eii = prob.ref_point[var[vi].index, 2];
                    tmp = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + var_to_dofs[vi][1];
                    if length(prob.ref_point[var[vi].index, 3]) > 1
                        tmp = tmp:(tmp+length(prob.ref_point[var[vi].index, 3])-1);
                    end
                    posind = [posind; tmp];
                    vals = [vals; prob.ref_point[var[vi].index, 3]];
                end
            end
            if length(vals) > 0
                if !rhs_only
                    A = identity_rows(A, posind, length(b));
                end
                b[posind] = vals;
            end
            
        else
            if prob.ref_point[var.index,1]
                eii = prob.ref_point[var.index, 2];
                posind = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + 1;
                if length(prob.ref_point[var.index, 3]) > 1
                    posind = posind:(posind+length(prob.ref_point[var[vi].index, 3])-1);
                else
                    posind = [posind];
                end
                if !rhs_only
                    A = identity_rows(A, posind, length(b));
                end
                b[posind] = prob.ref_point[var.index, 3];
            end
        end
    end
    
    bc_time = Base.Libc.time() - bc_time;
    
    if rhs_only
        if keep_geometric_factors
            return (b, saved_wdetj, saved_J)
        else
            return b;
        end
        
    else
        log_entry("Elemental loop time:     "*string(loop_time));
        log_entry("Boundary condition time: "*string(bc_time));
        if keep_geometric_factors
            return (A, b, saved_wdetj, saved_J);
        else
            return (A, b);
        end
    end
end

function assemble_rhs_only(var, linear, t=0.0, dt=0.0; keep_geometric_factors = false, saved_geometric_factors = nothing)
    return assemble(var, nothing, linear, t, dt; rhs_only=true, keep_geometric_factors = keep_geometric_factors, saved_geometric_factors = saved_geometric_factors)
end

# ######################################################
# # To be romoved. Fix nonlinear solve first.
# # assembles the A and b in Au=b
# function assemble_rhs_only(var, linear, t=0.0, dt=0.0; keep_geometric_factors = false, saved_geometric_factors = nothing)
#     Np = refel.Np;
#     nel = mesh_data.nel;
#     N1 = size(grid_data.allnodes,2);
#     multivar = typeof(var) <: Array;
#     maxvarindex = 0;
#     if multivar
#         # multiple variables being solved for simultaneously
#         dofs_per_node = 0;
#         var_to_dofs = [];
#         for vi=1:length(var)
#             tmp = dofs_per_node;
#             dofs_per_node += length(var[vi].symvar.vals);
#             push!(var_to_dofs, (tmp+1):dofs_per_node);
#             maxvarindex = max(maxvarindex,var[vi].index);
#         end
#     else
#         # one variable
#         dofs_per_node = length(var.symvar.vals);
#         maxvarindex = var.index;
#     end
#     Nn = dofs_per_node * N1;

#     b = zeros(Nn);
    
#     # Save geometric factors if desired
#     if keep_geometric_factors
#         saved_wdetj = [];
#         saved_J = [];
#     end
    
#     # Stiffness and mass are precomputed for uniform grid meshes
#     precomputed_mass_stiffness = config.mesh_type == UNIFORM_GRID && config.geometry == SQUARE
#     if precomputed_mass_stiffness
#         glb = grid_data.loc2glb[:,1];
#         xe = grid_data.allnodes[:,glb[:]];
#         (detJ, J) = geometric_factors(refel, xe);
#         wdetj = refel.wg .* detJ;
#         if keep_geometric_factors
#             saved_wdetj = wdetj;
#             saved_J = J;
#         end
#         if config.dimension == 1
#             (RQ1, RD1) = build_deriv_matrix(refel, J);
#             TRQ1 = RQ1';
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1)];
#         elseif config.dimension == 2
#             (RQ1, RQ2, RD1, RD2) = build_deriv_matrix(refel, J);
#             (TRQ1, TRQ2) = (RQ1', RQ2');
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2)];
#         else
#             (RQ1, RQ2, RQ3, RD1, RD2, RD3) = build_deriv_matrix(refel, J);
#             (TRQ1, TRQ2, TRQ3) = (RQ1', RQ2', RQ3');
#             stiffness = [(TRQ1 * diagm(wdetj) * RQ1) , (TRQ2 * diagm(wdetj) * RQ2) , (TRQ3 * diagm(wdetj) * RQ3)];
#         end
#         mass = (refel.Q)' * diagm(wdetj) * refel.Q;
#     else
#         stiffness = 0;
#         mass = 0;
#     end
    
#     # Elemental loop follows elemental ordering
#     for e=elemental_order;
#         glb = grid_data.loc2glb[:,e];       # global indices of this element's nodes for extracting values from var arrays
#         xe = grid_data.allnodes[:,glb[:]];  # coordinates of this element's nodes for evaluating coefficient functions
        
#         if !precomputed_mass_stiffness
#             if saved_geometric_factors === nothing
#                 (detJ, J) = geometric_factors(refel, xe);
#                 wdetj = refel.wg .* detJ;
                
#                 if keep_geometric_factors
#                     push!(saved_wdetj, wdetj);
#                     push!(saved_J, J);
#                 end
#             else
#                 wdetj = saved_geometric_factors[1][e];
#                 J = saved_geometric_factors[2][e];
#             end
#         end
        
#         rhsargs = (var, xe, glb, refel, wdetj, J, RHS, t, dt, stiffness, mass);

#         #linchunk = linear.func(args);  # get the elemental linear part
#         if dofs_per_node == 1
#             linchunk = linear.func(rhsargs);  # get the elemental linear part
#             b[glb] .+= linchunk;

#         elseif typeof(var) == Variable
#             # only one variable, but more than one dof
#             linchunk = linear.func(rhsargs);
#             insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

#         else
#             linchunk = linear.func(rhsargs);
#             insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);

#         end
#     end
    
#     # Boundary conditions
#     bidcount = length(grid_data.bids); # the number of BIDs
#     if dofs_per_node > 1
#         if multivar
#             dofind = 0;
#             for vi=1:length(var)
#                 for compo=1:length(var[vi].symvar.vals)
#                     dofind = dofind + 1;
#                     for bid=1:bidcount
#                         if prob.bc_type[var[vi].index, bid] == NO_BC
#                             # do nothing
#                         else
#                             b = dirichlet_bc_rhs_only(b, prob.bc_func[var[vi].index, bid][compo], grid_data.bdry[bid], t, dofind, dofs_per_node);
#                         end
#                     end
#                 end
#             end
#         else
#             for d=1:dofs_per_node
#                 #rows = ((d-1)*length(glb)+1):(d*length(glb));
#                 dofind = d;
#                 for bid=1:bidcount
#                     if prob.bc_type[var.index, bid] == NO_BC
#                         # do nothing
#                     else
#                         b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][d], grid_data.bdry[bid], t, d, dofs_per_node);
#                     end
#                 end
#             end
#         end
#     else
#         for bid=1:bidcount
#             if prob.bc_type[var.index, bid] == NO_BC
#                 # do nothing
#             else
#                 b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, bid][1], grid_data.bdry[bid], t);
#             end
#         end
#     end
    
#     # Reference points
#     if size(prob.ref_point,1) >= maxvarindex
#         if multivar
#             posind = zeros(Int,0);
#             vals = zeros(0);
#             for vi=1:length(var)
#                 if prob.ref_point[var[vi].index,1]
#                     eii = prob.ref_point[var[vi].index, 2];
#                     tmp = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + var_to_dofs[vi][1];
#                     if length(prob.ref_point[var[vi].index, 3]) > 1
#                         tmp = tmp:(tmp+length(prob.ref_point[var[vi].index, 3])-1);
#                     end
#                     posind = [posind; tmp];
#                     vals = [vals; prob.ref_point[var[vi].index, 3]];
#                 end
#             end
#             if length(vals) > 0
#                 b[posind] = vals;
#             end
            
#         else
#             if prob.ref_point[var.index,1]
#                 eii = prob.ref_point[var.index, 2];
#                 posind = (grid_data.glbvertex[eii[1], eii[2]] - 1)*dofs_per_node + 1;
#                 if length(prob.ref_point[var.index, 3]) > 1
#                     posind = posind:(posind+length(prob.ref_point[var[vi].index, 3])-1);
#                 end
#                 b[posind] = prob.ref_point[var.index, 3];
#             end
#         end
#     end

#     return b;
# end
# ##########################################################

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

function place_sol_in_vars(var, sol, stepper)
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
                if stepper.type == EULER_EXPLICIT
                    var[vi].values[compi,:] += sol[(compi+tmp):totalcomponents:end];
                else 
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                end
                tmp = tmp + 1;
            end
        end
    else
        components = length(var.symvar.vals);
        for compi=1:components
            if stepper.type == EULER_EXPLICIT
                var.values[compi,:] += sol[compi:components:end];
            else
                var.values[compi,:] = sol[compi:components:end];
            end
        end
    end
end

function get_var_vals(var)
    # place the variable values in a vector
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + length(var[vi].symvar.vals);
        end
        vect = zeros(totalcomponents * length(var[1].values[1,:]));
        for vi=1:length(var)
            components = length(var[vi].symvar.vals);
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
                tmp = tmp + 1;
            end
        end
    else
        components = length(var.symvar.vals);
        vect = zeros(components * length(var.values[1,:]));
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
end

end #module
