#=
# FV solver
=#
module FVSolver

export solve, nonlinear_solve

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
                SCALAR, VECTOR, TENSOR, SYM_TENSOR, VAR_ARRAY,
                LHS, RHS,
                LINEMESH, QUADMESH, HEXMESH
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data, grid_data, refel, time_stepper, elemental_order, genfunctions, indexers
import ..Femshop: Variable, Coefficient, GenFunction
import ..Femshop: GeometricFactors, geo_factors, geometric_factors, geometric_factors_face, build_deriv_matrix
import ..Femshop: FVInfo, fv_info, FV_cell_to_node, FV_node_to_cell

using LinearAlgebra, SparseArrays

include("fv_boundary.jl");

function linear_solve(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing, assemble_func=nothing)
    # If more than one variable
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        dofs_per_loop = 0;
        for vi=1:length(var)
            dofs_per_loop += length(var[vi].symvar);
            dofs_per_node += var[vi].total_components;
        end
    else
        # one variable
        dofs_per_loop = length(var.symvar);
        dofs_per_node = var.total_components;
    end
    nel = size(grid_data.loc2glb, 2);
    nfaces = size(grid_data.face2element, 2);
    Nn = dofs_per_node * nel;
    Nf = dofs_per_node * nfaces
    
    # Allocate arrays that will be used by assemble
    # These vectors will hold the integrated values(one per cell).
    # They will later be combined.
    sourcevec = zeros(Nn);
    fluxvec = zeros(Nn);
    facefluxvec = zeros(Nf);
    face_done = zeros(Int, nfaces); # Increment when the corresponding flux value is computed.
    allocated_vecs = [sourcevec, fluxvec, facefluxvec, face_done];
    
    if prob.time_dependent && !(stepper === nothing)
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        sol = get_var_vals(var);
        
        # allocate storage used by steppers
        if stepper.type == LSRK4
            tmppi = zeros(size(sol));
            tmpki = zeros(size(sol));
        elseif stepper.type == RK4
            tmpki = zeros(length(sol), stepper.stages);
        end
        
        start_t = Base.Libc.time();
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
                    
                    tmppi = get_var_vals(var, tmppi);
                    tmpki = zeros(size(sol));
                    for rki=1:stepper.stages
                        rktime = t + stepper.c[rki]*stepper.dt;
                        # p(i-1) is currently in u
                        
                        if assemble_func === nothing
                            sol = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rh, allocated_vecs, dofs_per_node, rktime, stepper.dt);
                        else
                            sol = assemble_using_loop_func(var, assemble_func, source_lhs, source_rhs, flux_lhs, flux_rh, allocated_vecs, dofs_per_node, dofs_per_loop, rktime, stepper.dt);
                        end
                        
                        
                        if rki == 1 # because a1 == 0
                            tmpki = stepper.dt .* sol;
                        else
                            tmpki = stepper.a[rki].*tmpki + stepper.dt.*sol;
                        end
                        tmppi = tmppi + stepper.b[rki].*tmpki
                        
                        place_sol_in_vars(var, tmppi, stepper);
                    end
                    
                else
                    # Explicit multi-stage methods: 
                    # x = x + dt*sum(bi*ki)
                    # ki = rhs(t+ci*dt, x+dt*sum(aij*kj)))   j < i
                    
                    # will hold the final result
                    sol = get_var_vals(var, sol);
                    # will be placed in var.values for each stage
                    tmpvals = sol;
                    for stage=1:stepper.stages
                        stime = t + stepper.c[stage]*stepper.dt;
                        
                        if assemble_func === nothing
                            tmpki[:,stage] = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, stime, stepper.dt);
                        else
                            tmpki[:,stage] = assemble_using_loop_function(var, assemble_func, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, stime, stepper.dt);
                        end
                        
                        tmpvals = sol;
                        for j=1:(stage-1)
                            if stepper.a[stage, j] > 0
                                tmpvals += stepper.dt * stepper.a[stage, j] .* tmpki[:,j];
                            end
                        end
                        
                        place_sol_in_vars(var, tmpvals, stepper);
                    end
                    for stage=1:stepper.stages
                        sol += stepper.dt * stepper.b[stage] .* tmpki[:, stage];
                    end
                    place_sol_in_vars(var, sol, stepper);
                end
                
            elseif stepper.type == EULER_EXPLICIT
                if assemble_func === nothing
                    sol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, t, stepper.dt);
                else
                    sol = sol .+ stepper.dt .* assemble_using_loop_function(var, assemble_func, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, stepper.dt);
                end
                
                place_sol_in_vars(var, sol, stepper);
                
            else
                printerr("Only explicit time steppers for FV are ready. TODO")
                return sol;
            end
            
            # ########## uncomment to return after one time step
            # return sol
            # ############
            
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
        
        log_entry("Stepping took "*string(end_t-start_t)*" seconds.");
        #log_entry("Stepping took "*string(end_t-start_t)*" seconds. ("*string(assemble_t)*" for assembly, "*string(linsolve_t)*" for linear solve)");
        #display(sol);
		# outfile = "linear_sol.txt"
		# open(outfile, "w") do f
  		# 	for ii in sol
    	# 		println(f, ii)
  		# 	end
		# end # the file f is automatically closed after this block finishes
        return sol;

    else
        # Does it make any sense to do this for time-independent problems?
        return sol;
    end
end

# macro loops(indices, ranges, content)
#     n = length(indices.args);
#     if n == 1
#         this_ind = indices.args[1];
#         this_range = ranges.args[1];
#         return esc(quote
#             for $this_ind in $this_range
#                 $content
#             end
#         end)
#     else
#         this_ind = indices.args[1];
#         this_range = ranges.args[1];
#         sub_ind = copy(indices);
#         sub_ind.args = sub_ind.args[2:end];
#         sub_range = copy(ranges);
#         sub_range.args = sub_range.args[2:end];
#         return esc(:(@loops([$this_ind], [$this_range], @loops($sub_ind, $sub_range, $content))));
#     end
# end

#

function assemble_using_loop_function(var, assemble_loops, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, dofs_per_loop = 1, t=0, dt=0)
    return assemble_loops.func(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, dofs_per_loop, t, dt);
end

function assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, t=0, dt=0)
    nel = size(grid_data.loc2glb, 2);
    
    # Label things that were allocated externally
    sourcevec = allocated_vecs[1];
    fluxvec = allocated_vecs[2];
    facefluxvec = allocated_vecs[3];
    face_done = allocated_vecs[4];
    
    face_done .= 0;
    
    # Elemental loop
    Threads.@threads for ei=1:nel
        eid = elemental_order[ei]; # The index of this element
        # Zero the result vectors for this element
        sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .= 0;
        fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .= 0;
        
        ##### Source integrated over the cell #####
        # Compute RHS volume integral
        if !(source_rhs === nothing)
            #sourceargs = prepare_args(var, eid, 0, RHS, "volume", t, dt); # (var, e, nodex, loc2glb, refel, detj, J, t, dt)
            sourceargs = (var, eid, 0, grid_data, geo_factors, fv_info, refel, t, dt);
            source = source_rhs.func(sourceargs) ./ geo_factors.volume[eid];
            # Add to global source vector
            sourcevec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] = source;
        end
        
        ##### Flux integrated over the faces #####
        # Loop over this element's faces.
        for i=1:refel.Nfaces
            fid = grid_data.element2face[i, eid];
            
            if !(flux_rhs === nothing)
                if face_done[fid] == 0
                    face_done[fid] = 1; # Possible race condition, but in the worst case it will be computed twice.
                    
                    #fluxargs = prepare_args(var, eid, fid, RHS, "surface", t, dt); #(var, (e, neighbor), refel, vol_loc2glb, nodex, cellx, frefelind, facex, face2glb, normal, fdetj, face_area, (J, vol_J_neighbor), t, dt);
                    fluxargs = (var, eid, fid, grid_data, geo_factors, fv_info, refel, t, dt);
                    flux = flux_rhs.func(fluxargs) .* geo_factors.area[fid];
                    # Add to global flux vector for faces
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] .= flux;
                    # Combine all flux for this element
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .+= flux ./ geo_factors.volume[eid];
                    
                else
                    # This flux has either been computed or is being computed by another thread.
                    # The state will need to be known before paralellizing, but for now assume it's complete.
                    fluxvec[((eid-1)*dofs_per_node + 1):(eid*dofs_per_node)] .-= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] ./ geo_factors.volume[eid];
                end
            end
            
            # Boundary conditions are applied to flux
            fbid = grid_data.facebid[fid]; # BID of this face
            if fbid > 0
                facex = grid_data.allnodes[:, grid_data.face2glb[:,1,fid]];  # face node coordinates
                
                if typeof(var) <: Array
                    dofind = 0;
                    for vi=1:length(var)
                        for compo=1:length(var[vi].symvar)
                            dofind = dofind + 1;
                            if prob.bc_type[var[vi].index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                                Qvec = Qvec ./ geo_factors.area[fid];
                                bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                                
                                fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                                facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                            else
                                printerr("Unsupported boundary condition type: "*prob.bc_type[var[vi].index, fbid]);
                            end
                        end
                    end
                else
                    for d=1:dofs_per_node
                        dofind = d;
                        if prob.bc_type[var.index, fbid] == NO_BC
                            # do nothing
                        elseif prob.bc_type[var.index, fbid] == FLUX
                            # compute the value and add it to the flux directly
                            Qvec = (refel.surf_wg[grid_data.faceRefelInd[1,fid]] .* geo_factors.face_detJ[fid])' * (refel.surf_Q[grid_data.faceRefelInd[1,fid]])[:, refel.face2local[grid_data.faceRefelInd[1,fid]]]
                            Qvec = Qvec ./ geo_factors.area[fid];
                            bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* geo_factors.area[fid];
                            
                            fluxvec[(eid-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) ./ geo_factors.volume[eid];
                            facefluxvec[(fid-1)*dofs_per_node + dofind] = bflux;
                        else
                            printerr("Unsupported boundary condition type: "*prob.bc_type[var.index, fbid]);
                        end
                    end
                end
            end# BCs
            
        end# face loop
    end# element loop
    
    return sourcevec + fluxvec;
end

# Gathers all info needed to be passed in args
function prepare_args(var, ei, fi, lorr, vors, t, dt)
    if vors == "volume"
        # nodes and maps
        loc2glb = grid_data.loc2glb[:,ei];               # elemental node global index
        nodex = grid_data.allnodes[:,glb[:]];       # elemental node coordinates
        
        # geometric factors
        detj = geo_factors.detJ[ei];
        J = geo_factors.J[ei];
        
        args = (var, ei, nodex, loc2glb, refel, detj, J, t, dt);
        
    else # surface
        if grid_data.face2element[1, fi] == ei
            # The normal points out of e
            normal = grid_data.facenormals[:, fi];
            neighbor = grid_data.face2element[2, fi];
            frefelind = [grid_data.faceRefelInd[1,fi], grid_data.faceRefelInd[2,fi]]; # refel based index of face in both elements
        else
            # The normal points into e, reverse it
            normal = -grid_data.facenormals[:, fi];
            neighbor = grid_data.face2element[1, fi];
            frefelind = [grid_data.faceRefelInd[2,fi], grid_data.faceRefelInd[1,fi]]; # refel based index of face in both elements
        end
        
        face2glb = grid_data.face2glb[:,:,fi];         # global index for face nodes for each side of each face
        facex = grid_data.allnodes[:, face2glb[:, 1]];  # face node coordinates
        
        # geometric factors
        fdetj = geo_factors.face_detJ[fi];
        face_area = geo_factors.area[fi];
        
        if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
            neighbor = ei;
        end
        els = (ei, neighbor);
        
        vol_J = (geo_factors.J[ei], geo_factors.J[neighbor]);
        
        vol_loc2glb = (grid_data.loc2glb[:,ei], grid_data.loc2glb[:, neighbor]); # volume local to global
        
        nodex = (grid_data.allnodes[:,vol_loc2glb[1][:]], grid_data.allnodes[:,vol_loc2glb[2][:]]); # volume node coordinates
        cellx = (fv_info.cellCenters[:, ei], fv_info.cellCenters[:, neighbor]); # cell center coordinates
        
        args = (var, els, refel, vol_loc2glb, nodex, cellx, frefelind, facex, face2glb, normal, fdetj, face_area, vol_J, t, dt);
        
    end
    
    return args
end

# Insert the single dof into the greater vector
function insert_linear!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));
        
        b[ind] = b[ind] + bel[ind2];
    end
end

function place_sol_in_vars(var, sol, stepper)
    # place the values in the variable value arrays
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                if stepper.type == EULER_EXPLICIT
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                else 
                    var[vi].values[compi,:] = sol[(compi+tmp):totalcomponents:end];
                end
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        for compi=1:components
            if stepper.type == EULER_EXPLICIT
                var.values[compi,:] = sol[compi:components:end];
            else
                var.values[compi,:] = sol[compi:components:end];
            end
        end
    end
end

function get_var_vals(var, vect=nothing)
    # place the variable values in a vector
    if typeof(var) <: Array
        tmp = 0;
        totalcomponents = 0;
        for vi=1:length(var)
            totalcomponents = totalcomponents + var[vi].total_components;
        end
        if vect === nothing
            vect = zeros(totalcomponents * length(var[1].values[1,:]));
        end
        
        for vi=1:length(var)
            components = var[vi].total_components;
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
                tmp = tmp + 1;
            end
        end
    else
        components = var.total_components;
        if vect === nothing
            vect = zeros(components * length(var.values[1,:]));
        end
        for compi=1:components
            vect[compi:components:end] = var.values[compi,:];
        end
    end
    
    return vect;
end

end #module
