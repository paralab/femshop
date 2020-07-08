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
            LINEMESH, QUADMESH
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data, grid_data, loc2glb, refel, time_stepper
import ..Femshop: Variable, Coefficient, GenFunction

using LinearAlgebra, SparseArrays

include("cg_operators.jl");
include("cg_boundary.jl");
include("elemental_matrix.jl");

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
        variables[vind].values = zeros(size(grid_data.allnodes)[1]);
        if vind <= length(prob.initial)
            if prob.initial[vind] != nothing
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
    Nn = size(grid_data.allnodes)[1];
    
    b = zeros(Nn);
    A = spzeros(Nn, Nn);
    
    for e=1:nel;
        nv = mesh_data.nv[e];
        gis = zeros(Int, nv);
        for vi=1:nv
            gis[vi] = mesh_data.invind[mesh_data.elements[e,vi]]; # mesh indices of element's vertices
        end
        
        vx = mesh_data.nodes[gis,:];        # coordinates of element's vertices
        glb = loc2glb[e,:];                 # global indices of this element's nodes for extracting values from var arrays
        xe = grid_data.allnodes[glb[:],:];  # coordinates of this element's nodes for evaluating coefficient functions
        
        args = (var, xe, glb, refel, RHS, t, dt);
        linchunk = linear.func.func(args);  # get the elemental linear part
        b[glb] .+= linchunk;
        
        args = (var, xe, glb, refel, LHS, t, dt);
        bilinchunk = bilinear.func.func(args); # the elemental bilinear part
        A[glb, glb] .+= bilinchunk;         # This will be very inefficient for sparse A
    end
    
    # Just for testing. This should really be a loop over BIDs with corresponding calls
    (A, b) = dirichlet_bc(A, b, prob.bc_func[var.index, 1], grid_data.bdry[1,:], t);

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
        linchunk = linear.func.func(args);  # get the elemental linear part
        b[glb] .+= linchunk;
    end
    
    # Just for testing. This should really be a loop over BIDs with corresponding calls
    b = dirichlet_bc_rhs_only(b, prob.bc_func[var.index, 1], grid_data.bdry[1,:], t);

    return b;
end

end #module