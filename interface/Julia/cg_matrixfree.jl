#=
Functions used by the CGSolver for matrix free solutions.
=#

# Note: This uses the conjugate gradient iterative method,
# which assumes an SPD matrix.
function solve_matrix_free_sym(var, bilinear, linear, stepper=nothing)
    start_time = time_ns();
    tol = config.linalg_matfree_tol;
    maxiters = config.linalg_matfree_max;
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
    
    if dofs_per_node > 1
        println("Oops. Still working on multi DOF matrix-free");
        return zeros(Nn);
    end
    
    if prob.time_dependent && !(stepper === nothing)
        #TODO
    else
        # Use regular rhs assembly
        b = assemble_rhs_only(var, linear, 0, 0);
        normb = norm(b, Inf);
        
        x = zeros(Nn);
        
        # Do initial matvec
        #Ax = elem_matvec(x, bilinear, dofs_per_node, var);
        # nevermind, this will just be zeros
        
        r0 = copy(b); # = b - Ax;
        p = copy(r0);
        
        iter = 0;
        err = 1;
        while iter < maxiters && err > tol
            iter = iter+1;
            
            Ap = elem_matvec(p,bilinear, dofs_per_node, var);
            alpha = dot(r0,r0) / dot(p,Ap);
            
            x = x .+ alpha.*p;
            r1 = r0 .- alpha.*Ap;
            
            beta = dot(r1,r1) / dot(r0, r0);
            p = r1 .+ beta.*p;
            
            r0 = copy(r1);
            err = norm(r0, Inf)/normb;
            
            if iter%50 == 0
                println("iteration "*string(iter)*": res = "*string(err));
            end
        end
        
        println("Converged to "*string(err)*" in "*string(iter)*" iterations");
        log_entry("Converged to "*string(err)*" in "*string(iter)*" iterations");
        
        total_time = time_ns() - start_time;
        log_entry("Linear solve took "*string(total_time/1e9)*" seconds");
        
        # Rearrange the values to match the regular solve
        x = block_to_interlace(x,N1, dofs_per_node);
        
        return x;
    end
    
    
end

# Note: This uses the stabilized biconjugate gradient method which works for nonsymmetric matrices
function solve_matrix_free_asym(var, bilinear, linear, stepper=nothing)
    start_time = time_ns();
    tol = config.linalg_matfree_tol;
    maxiters = config.linalg_matfree_max;
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
    
    if prob.time_dependent && !(stepper === nothing)
        #TODO
    else
        # Use regular rhs assembly
        b = assemble_rhs_only(var, linear, 0, 0);
        normb = norm(b, Inf);
        
        x = zeros(Nn);
        
        r0 = copy(b); # = b - Ax; but x is initially zeros
        rs = copy(b);
        rho0 = 1;
        alpha = 1;
        w0 = 1;
        p = zeros(Nn);
        Ap = zeros(Nn);
        As = zeros(Nn);
        
        iter = 0;
        err = 1;
        while iter < maxiters && err > tol
            iter = iter+1;
            
            rho1 = dot(rs,r0);
            beta = (rho1*alpha) / (rho0*w0);
            p = r0 + beta*(p - w0*Ap);
            
            Ap = elem_matvec(p,bilinear, dofs_per_node, var);
            
            alpha = rho1/dot(rs,Ap);
            s = r0 - alpha*Ap;
            
            As = elem_matvec(s,bilinear, dofs_per_node, var);
            
            w1 = dot(s,As) / dot(As,As);
            x = x + alpha*p + w1*s;
            r0 = s - w1*As;
            
            err = norm(r0, Inf)/normb;
            
            if iter%50 == 0
                println("iteration "*string(iter)*": res = "*string(err));
            end
        end
        
        println("Converged to "*string(err)*" in "*string(iter)*" iterations");
        log_entry("Converged to "*string(err)*" in "*string(iter)*" iterations");
        
        total_time = time_ns() - start_time;
        log_entry("Linear solve took "*string(total_time/1e9)*" seconds");
        
        return x;
    end
end

# Does the elemental matvec Ax=b
# x is input, b is output, A is from bilinear
function elem_matvec(x, bilinear, dofs_per_node, var, t = 0.0, dt = 0.0)
    Np = refel.Np;
    nel = mesh_data.nel;
    multivar = typeof(var) <: Array;
    if multivar
        # multiple variables being solved for simultaneously
        dofs_per_node = 0;
        var_to_dofs = [];
        for vi=1:length(var)
            tmp = dofs_per_node;
            dofs_per_node += length(var[vi].symvar.vals);
            push!(var_to_dofs, (tmp+1):dofs_per_node);
            
            # # check for neumann bcs
            # for bi=1:length(prob.bc_type[vi,:])
            #     if prob.bc_type[vi,bi] == NEUMANN
                    
            #     end
            # end
        end
        
    else
        # one variable
        dofs_per_node = length(var.symvar.vals);
    end
    
    Ax = zeros(size(x));
    #Elemental loop follows elemental ordering
    for e=elemental_order;
        gis = grid_data.glbvertex[e,:];
        vx = grid_data.allnodes[gis,:];         # coordinates of element's vertices
        glb = grid_data.loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions
        
        subx = extract_linear(x, glb, dofs_per_node);
        
        lhsargs = (var, xe, glb, refel, LHS, t, dt);
        bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
        # Neumann bcs need to be applied to this
        
        
        if dofs_per_node == 1
            #bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
            Ax[glb] = Ax[glb] + bilinchunk * subx;
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            #bilinchunk = bilinear.func(lhsargs);
            insert_linear_matfree!(Ax, bilinchunk*subx, glb, 1:dofs_per_node, dofs_per_node);
        else
            #bilinchunk = bilinear.func(lhsargs);
            insert_linear_matfree!(Ax, bilinchunk*subx, glb, 1:dofs_per_node, dofs_per_node);
        end
    end
    
    # Apply boudary conditions
    bidcount = length(grid_data.bids); # the number of BIDs
    if dofs_per_node > 1
        if multivar
            rowoffset = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar.vals)
                    rowoffset = rowoffset + 1;
                    for bid=1:bidcount
                        Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid], rowoffset, dofs_per_node);
                    end
                end
            end
        else
            for d=1:dofs_per_node
                #rows = ((d-1)*length(glb)+1):(d*length(glb));
                rowoffset = (d-1)*Np;
                for bid=1:bidcount
                    Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid], d, dofs_per_node);
                end
            end
        end
    else
        for bid=1:bidcount
            Ax = dirichlet_bc_matfree(Ax, x, grid_data.bdry[bid]);
        end
    end
    
    return Ax;
end

# sets b[bdry] = x[bdry]
function dirichlet_bc_matfree(b, x, bdryind, dofind=1, totaldofs=1)
    if totaldofs > 1
        #bdry = copy(bdryind);
        #bdry = (bdry .- 1) .* totaldofs .+ dofind;
        #b[bdry] = x[bdry];
        N1 = length(grid_data.allnodes[:,1]);
        ind = bdryind .+ (dofind-1)*N1;
        
        b[ind] = x[ind];
    else
        b[bdryind] = x[bdryind];
    end
    
    return b;
end

function extract_linear(b, glb, dofs)
    if dofs == 1
        return b[glb];
    else
        N = length(glb);
        part = zeros(N*dofs);
        
        for d=1:dofs
            part[((d-1)*N+1):(d*N)] = b[glb.+(d-1)];
        end
        
        return part;
    end
end

# Inset the single dof into the greater construct
function insert_linear_matfree!(b, bel, glb, dof, Ndofs)
    # group nodal dofs
    for d=1:length(dof)
        #ind = glb.*Ndofs .- (Ndofs-dof[d]);
        ind2 = ((d-1)*length(glb)+1):(d*length(glb));
        
        b[ind2] = b[ind2] + bel[ind2];
    end
end

function block_to_interlace(x, Np, dofs)
    newx = similar(x);
    for d=1:dofs
        newx[d:dofs:end] = x[((d-1)*Np+1):(d*Np)];
    end
    return newx;
end
function interlace_to_block(x, Np, dofs)
    newx = similar(x);
    for d=1:dofs
        newx[((d-1)*Np+1):(d*Np)] = x[d:dofs:end];
    end
    return newx;
end