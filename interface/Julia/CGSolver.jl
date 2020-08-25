#=
# CG solver
=#
module CGSolver

export init_cgsolver, solve

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
            NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
            NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICHOLSON, RK4, LSRK4,
            ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, NEUMANN, ROBIN,
            MSH_V2, MSH_V4,
            SCALAR, VECTOR, TENSOR, SYM_TENSOR,
            LHS, RHS,
            LINEMESH, QUADMESH, HEXMESH
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data, grid_data, loc2glb, refel, time_stepper
import ..Femshop: Variable, Coefficient, GenFunction

using LinearAlgebra, SparseArrays

#include("cg_operators.jl");
include("cg_boundary.jl");
#include("elemental_matrix.jl");

function init_cgsolver()
    dim = config.dimension;
    if dim == 1
        setup1D();
    elseif dim == 2
        setup2D();
    elseif dim == 3
        setup3D();
    end
    # build initial conditions
    for vind=1:length(variables)
        if vind <= length(prob.initial)
            if prob.initial[vind] != nothing
                variables[vind].values = zeros(size(grid_data.allnodes)[1]);
                # Set initial condition
                for ni=1:size(grid_data.allnodes)[1]
                    if dim == 1
                        variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni],0,0,0);
                    elseif dim == 2
                        variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],0,0);
                    elseif dim == 3
                        variables[vind].values[ni] = prob.initial[vind].func(grid_data.allnodes[ni,1],grid_data.allnodes[ni,2],grid_data.allnodes[ni,3],0);
                    end
                end
                variables[vind].ready = true;
                log_entry("Built initial conditions for: "*string(variables[vind].symbol));
            end
        end
    end
end

function setup1D()
    # This may be unnecessary. I'll keep this for now in case I need it.
end
function setup2D()
    
end
function setup3D()
    
end

function solve(var, bilinear, linear, stepper=nothing)
    if prob.time_dependent && !(stepper === nothing)
        #TODO time dependent coefficients
        assemble_t = @elapsed((A, b) = assemble(var, bilinear, linear, 0, stepper.dt));
        log_entry("Assembly took "*string(assemble_t)*" seconds");
        
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        t = 0;
        start_t = Base.Libc.time();
        for i=1:stepper.Nsteps
            b = assemble_rhs_only(var, linear, t, stepper.dt);
            var.values = A\b;
            t += stepper.dt;
        end
        end_t = Base.Libc.time();
        
        log_entry("Solve took "*string(end_t-start_t)*" seconds");
        #display(sol);
        return var.values;
        
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
    A = spzeros(Nn, Nn);
    
    # The elemental assembly loop
    for e=1:nel;
        nv = mesh_data.nv[e];
        gis = zeros(Int, nv);
        for vi=1:nv
            gis[vi] = mesh_data.invind[mesh_data.elements[e,vi]]; # mesh indices of element's vertices
        end
        
        vx = mesh_data.nodes[gis,:];        # coordinates of element's vertices
        glb = loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions
        
        # The linear part. Compute the elemental linear part for each dof
        rhsargs = (var, xe, glb, refel, RHS, t, dt);
        lhsargs = (var, xe, glb, refel, LHS, t, dt);
        if dofs_per_node == 1
            linchunk = linear.func(rhsargs);  # get the elemental linear part
            b[glb] .+= linchunk;
            
            bilinchunk = bilinear.func(lhsargs); # the elemental bilinear part
            A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
        elseif typeof(var) == Variable
            # only one variable, but more than one dof
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);
            
            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear!(A, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        else
            linchunk = linear.func(rhsargs);
            insert_linear!(b, linchunk, glb, 1:dofs_per_node, dofs_per_node);
            
            bilinchunk = bilinear.func(lhsargs);
            insert_bilinear!(A, bilinchunk, glb, 1:dofs_per_node, dofs_per_node);
        end
    end
    
    # Just for testing. This should really be a loop over BIDs with corresponding calls
    if dofs_per_node > 1
        if multivar
            rowoffset = 0;
            for vi=1:length(var)
                for compo=1:length(var[vi].symvar.vals)
                    rowoffset = rowoffset + 1;
                    (A, b) = dirichlet_bc(A, b, prob.bc_func[var[vi].index, 1][compo], grid_data.bdry[1,:], t, rowoffset, dofs_per_node);
                end
            end
        else
            for d=1:dofs_per_node
                #rows = ((d-1)*length(glb)+1):(d*length(glb));
                rowoffset = (d-1)*Np;
                (A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, 1][d], grid_data.bdry[1,:], t, d, dofs_per_node);
            end
        end
    else
        (A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, 1], grid_data.bdry[1,:], t);
    end
    
    return (A, b);
end

# assembles the A and b in Au=b
function assemble_rhs_only(var, linear, t=0.0, dt=0.0)
    Np = refel.Np;
    nel = mesh_data.nel;
    Nn = size(grid_data.allnodes)[1];
    
    b = zeros(Nn);
    
    for e=1:nel;
        glb = loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions
        
        args = (var, xe, glb, refel, RHS, t, dt);
        linchunk = linear.func(args);  # get the elemental linear part
        b[glb] .+= linchunk;
    end
    
    # Just for testing. This should really be a loop over BIDs with corresponding calls
    b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, 1], grid_data.bdry[1,:], t);

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

function insert_bilinear!(a, ael, glb, dof, Ndofs)
    # group nodal dofs
    for di=1:length(dof)
        indi = glb.*Ndofs .- (Ndofs-dof[di]);
        indi2 = ((di-1)*length(glb)+1):(di*length(glb));
        for dj=1:length(dof)
            indj = glb.*Ndofs .- (Ndofs-dof[dj]);
            indj2 = ((dj-1)*length(glb)+1):(dj*length(glb));
            
            a[indi, indj] = a[indi, indj] + ael[indi2, indj2];
        end
    end
end

end #module