#=
# FV solver
=#
module FVSolver

export init_fvsolver, solve, nonlinear_solve

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
import ..Femshop: GeometricFactors, geo_factors, geometric_factors, geometric_factors_face, build_deriv_matrix
import ..Femshop: FVInfo, fv_info, FV_cell_to_node, FV_node_to_cell

using LinearAlgebra, SparseArrays

include("fv_boundary.jl");

function init_fvsolver()
    dim = config.dimension;

    # build initial conditions
    for vind=1:length(variables)
        if vind <= length(prob.initial)
            if prob.initial[vind] != nothing
                #variables[vind].values = zeros(size(variables[vind].values));
                
                # Set initial condition
                if typeof(prob.initial[vind]) <: Array 
                    nodal_values = zeros(length(prob.initial[vind]), size(grid_data.allnodes,2));
                    for ci=1:length(prob.initial[vind])
                        for ni=1:size(grid_data.allnodes,2)
                            if dim == 1
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[ni],0,0,0);
                            elseif dim == 2
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],0,0);
                            elseif dim == 3
                                nodal_values[ci,ni] = prob.initial[vind][ci].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],grid_data.allnodes[3,ni],0);
                            end
                        end
                    end
                    # compute cell averages using nodal values
                    nel = size(grid_data.loc2glb, 2);
                    for ei=1:nel
                        e = elemental_order[ei];
                        glb = grid_data.loc2glb[:,e];
                        vol = geo_factors.volume[e];
                        detj = geo_factors.detJ[e];
                        
                        for ci=1:length(prob.initial[vind])
                            variables[vind].values[ci,e] = detj / vol * (refel.wg' * refel.Q * (nodal_values[ci,glb][:]))[1];
                        end
                    end
                    
                else
                    nodal_values = zeros(size(grid_data.allnodes,2));
                    for ni=1:size(grid_data.allnodes,2)
                        if dim == 1
                            nodal_values[ni] = prob.initial[vind].func(grid_data.allnodes[ni],0,0,0);
                        elseif dim == 2
                            nodal_values[ni] = prob.initial[vind].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],0,0);
                        elseif dim == 3
                            nodal_values[ni] = prob.initial[vind].func(grid_data.allnodes[1,ni],grid_data.allnodes[2,ni],grid_data.allnodes[3,ni],0);
                        end
                    end
                    # compute cell averages using nodal values
                    nel = size(grid_data.loc2glb, 2);
                    for ei=1:nel
                        e = elemental_order[ei];
                        glb = grid_data.loc2glb[:,e];
                        vol = geo_factors.volume[e];
                        detj = geo_factors.detJ[e];
                        
                        variables[vind].values[e] = detj / vol * (refel.wg' * refel.Q * nodal_values[glb])[1];
                    end
                end
                
                variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(variables[vind].symbol));
            end
        end
    end
end

function linear_solve(var, source_lhs, source_rhs, flux_lhs, flux_rhs, stepper=nothing)
    # If more than one variable
    if typeof(var) <: Array
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
        end
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
    end
    nel = size(grid_data.loc2glb, 2);
    Nn = dofs_per_node * nel;
    Nf = dofs_per_node * size(grid_data.face2element, 2);
    
    # Allocate arrays that will be used by assemble
    # These vectors will hold the integrated values(one per cell).
    # They will later be combined and interpolated to nodes.
    sourcevec = zeros(Nn);
    fluxvec = zeros(Nn);
    facefluxvec = zeros(Nf);
    face_done = zeros(Bool, Nf); # Set to true when the corresponding fluxvec value is computed.
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
                        
                        sol = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rh, allocated_vecs, dofs_per_node, rktime, stepper.dt);
                        
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
                        
                        tmpki[:,stage] = assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, stime, stepper.dt);
                        
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
                sol = sol .+ stepper.dt .* assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node, t, stepper.dt);
                
                place_sol_in_vars(var, sol, stepper);
                
            else
                printerr("Only explicit time steppers for FV. TODO")
                return sol;
            end
            
            ########### uncomment to return after one time step
            #return sol
            #############
            
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

function assemble(var, source_lhs, source_rhs, flux_lhs, flux_rhs, allocated_vecs, dofs_per_node=1, t=0, dt=0)
    nel = size(grid_data.loc2glb, 2);
    dim = size(grid_data.allnodes, 1);
    
    # name things that were allocated externally
    sourcevec = allocated_vecs[1];
    fluxvec = allocated_vecs[2];
    facefluxvec = allocated_vecs[3];
    face_done = allocated_vecs[4];
    
    face_done .= false;
    
    # Elemental loop
    for ei=1:nel
        e = elemental_order[ei];                    # This element index
        
        ##### Source integrated over the cell #####
        
        # nodes and maps
        glb = grid_data.loc2glb[:,e];               # elemental node global index
        nodex = grid_data.allnodes[:,glb[:]];       # elemental node coordinates
        
        # geometric factors
        detj = geo_factors.detJ[e];
        J = geo_factors.J[e];
        inv_vol = 1/geo_factors.volume[e];
        
        if source_rhs === nothing
            source = zeros(dofs_per_node);
        else
            sourceargs = (var, e, nodex, glb, refel, detj, J, t, dt);
            source = source_rhs.func(sourceargs) .* inv_vol;
        end
        
        if dofs_per_node > 1
            sourcevec[((e-1)*dofs_per_node + 1):(e*dofs_per_node)] = source;
        else
            sourcevec[e] = source[1];
        end
        
        ##### Flux integrated over the faces #####
        
        nfaces = refel.Nfaces;                      # number of faces
        faces = grid_data.element2face[1:nfaces, e];# face index
        fluxvec[((e-1)*dofs_per_node + 1):(e*dofs_per_node)] .= 0;
        # Loop over this element's faces.
        for i=1:nfaces
            fid = faces[i];
            if face_done[fid] == false
                face_done[fid] = true; # Possible race condition, but in the worst case it will be computed twice.
                
                if grid_data.face2element[1, fid] == e
                    # The normal points out of e
                    normal = grid_data.facenormals[:, fid];
                    neighbor = grid_data.face2element[2, fid];
                    frefelind = [grid_data.faceRefelInd[1,fid], grid_data.faceRefelInd[2,fid]]; # refel based index of face in both elements
                else
                    # The normal points into e, reverse it
                    normal = -grid_data.facenormals[:, fid];
                    neighbor = grid_data.face2element[1, fid];
                    frefelind = [grid_data.faceRefelInd[2,fid], grid_data.faceRefelInd[1,fid]]; # refel based index of face in both elements
                end
                
                face2glb = grid_data.face2glb[:,:,fid];         # global index for face nodes for each side of each face
                facex = grid_data.allnodes[:, face2glb[:, 1]];  # face node coordinates
                
                # geometric factors
                fdetj = geo_factors.face_detJ[fid];
                face_area = geo_factors.area[fid];
                
                if neighbor == 0 # This is a boundary face. For now, just compute as if neighbor is identical. BCs handled later.
                    neighbor = e;
                end
                
                vol_J_neighbor = geo_factors.J[neighbor];
                vol_loc2glb = (glb, grid_data.loc2glb[:, neighbor]); # volume local to global
                
                nodex = (grid_data.allnodes[:,glb[:]], grid_data.allnodes[:,vol_loc2glb[2][:]]); # volume node coordinates
                cellx = (fv_info.cellCenters[:, e], fv_info.cellCenters[:, neighbor]); # cell center coordinates
                
                fluxargs = (var, (e, neighbor), refel, vol_loc2glb, nodex, cellx, frefelind, facex, face2glb, normal, fdetj, face_area, (J, vol_J_neighbor), t, dt);
                
                flux = flux_rhs.func(fluxargs) .* face_area;
                if dofs_per_node > 1
                    facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] = flux;
                    # Combine with source
                    fluxvec[((e-1)*dofs_per_node + 1):(e*dofs_per_node)] += flux .* inv_vol;
                else
                    facefluxvec[fid] = flux;
                    # Combine with source
                    fluxvec[e] += flux .* inv_vol;
                end
                
            else
                # This flux has either been computed or is being computed by another thread.
                # The state will need to be known before paralellizing, but for now assume it's complete.
                if dofs_per_node > 1
                    fluxvec[((e-1)*dofs_per_node + 1):(e*dofs_per_node)] -= facefluxvec[((fid-1)*dofs_per_node + 1):(fid*dofs_per_node)] .* inv_vol;
                else
                    fluxvec[e] -= facefluxvec[fid] .* inv_vol;
                end
            end
            
            # Boundary conditions are applied to flux
            fbid = grid_data.facebid[fid]; # BID of this face
            if fbid > 0
                facex = grid_data.allnodes[:, grid_data.face2glb[:,1,fid]];  # face node coordinates
                
                if typeof(var) <: Array
                    dofind = 0;
                    for vi=1:length(var)
                        for compo=1:length(var[vi].symvar.vals)
                            dofind = dofind + 1;
                            if prob.bc_type[var[vi].index, fbid] == NO_BC
                                # do nothing
                            elseif prob.bc_type[var[vi].index, fbid] == FLUX
                                # compute the value and add it to the flux directly
                                Qvec = (refel.surf_wg[frefelind[1]] .* fdetj)' * (refel.surf_Q[frefelind[1]])[:, refel.face2local[frefelind[1]]]
                                Qvec = Qvec ./ face_area
                                bflux = FV_flux_bc_rhs_only(prob.bc_func[var[vi].index, fbid][compo], facex, Qvec, t, dofind, dofs_per_node) .* face_area;
                                
                                fluxvec[(e-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) .* inv_vol;
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
                            Qvec = (refel.surf_wg[frefelind[1]] .* fdetj)' * (refel.surf_Q[frefelind[1]])[:, refel.face2local[frefelind[1]]]
                            Qvec = Qvec ./ face_area
                            bflux = FV_flux_bc_rhs_only(prob.bc_func[var.index, fbid][d], facex, Qvec, t, dofind, dofs_per_node) .* face_area;
                            
                            fluxvec[(e-1)*dofs_per_node + dofind] += (bflux - facefluxvec[(fid-1)*dofs_per_node + dofind]) .* inv_vol;
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
            totalcomponents = totalcomponents + length(var[vi].symvar.vals);
        end
        for vi=1:length(var)
            components = length(var[vi].symvar.vals);
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
        components = length(var.symvar.vals);
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
            totalcomponents = totalcomponents + length(var[vi].symvar.vals);
        end
        if vect === nothing
            vect = zeros(totalcomponents * length(var[1].values[1,:]));
        end
        
        for vi=1:length(var)
            components = length(var[vi].symvar.vals);
            for compi=1:components
                vect[(compi+tmp):totalcomponents:end] = var[vi].values[compi,:];
                tmp = tmp + 1;
            end
        end
    else
        components = length(var.symvar.vals);
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
