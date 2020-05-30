#=
# DG solver
=#
module DGSolver

export init_dgsolver, solve

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data

using LinearAlgebra

include("refel.jl");
include("time_steppers.jl");

# Globals
refel = nothing;
allnodes = [];
rx = [];
J = [];
Fmask = [];
normals = [];
Fscale = [];
e2e = [];
e2f = [];
vmapm = [];
vmapp = [];
vmapb = [];
mapb = [];
vmapI = [];
vmapO = [];
mapI = [];
mapO = [];

function init_dgsolver()
    if config.dimension == 1
        setup1D();
        log_entry("Set up DG solver.");
        # build initial conditions
        for vind=1:length(variables)
            variables[vind].values = zeros(size(allnodes));
            if vind <= length(prob.initial)
                if prob.initial[vind] != nothing
                    # Set initial condition
                    for ni=1:refel.Np
                        for ej=1:mesh_data.nel
                            variables[vind].values[ni,ej] = prob.initial[vind].func(allnodes[ni,ej]);
                        end
                    end
                    variables[vind].ready = true;
                    log_entry("Built initial conditions for: "*string(variables[vind].symbol));
                end
            end
        end
    else
        #not ready
    end
end

function setup1D()
    nel = mesh_data.nel;
    global refel = build_refel(config.dimension, config.basis_order_min, 2, LOBATTO);
      
    global allnodes = ones(refel.Np,1)*mesh_data.nodes[mesh_data.elements[:,1]'] .+
                  0.5 .*(refel.r.+1).*(mesh_data.nodes[mesh_data.elements[:,2]'] .-
                  mesh_data.nodes[mesh_data.elements[:,1]']);
    global J = refel.Dr*allnodes;
    global rx = 1 ./J;
    global Fmask = [1;refel.Np];
    global normals = ones(refel.Nfaces, nel);
    normals[1,:] = -ones(nel);
    global Fscale = 1 ./J[Fmask,:];
    # The rest of this assumes simple left to right ordering of elements.
    # Should do it differently.
    global e2e = zeros(Int, nel, refel.Nfaces);
    e2e[2:end,1] = 1:(nel-1);
    e2e[1,1] = 1;
    e2e[:,2] = 1:nel;
    global e2f = ones(Int, nel, refel.Nfaces);
    e2f[2:end,1] = 2 .*e2f[2:end,1];
    e2f[end,2] = 2;
    global vmapm = zeros(Int, nel*refel.Nfaces);
    global vmapp = zeros(Int, nel*refel.Nfaces);
    for i=1:nel
        vmapm[i*2-1] = (i-1)*refel.Np + 1;
        vmapm[i*2] = i*refel.Np;
        vmapp[i*2-1] = vmapm[i*2-1] - 1;
        vmapp[i*2] = vmapm[i*2] + 1;
    end
    vmapp[1] = vmapm[1];
    vmapp[end] = vmapm[end];
    global mapb = [1; nel*refel.Nfaces];
    global vmapb = vmapm[mapb];
    global vmapI = 1;
    global vmapO = nel*refel.Np;
    global mapI = 1;
    global mapO = nel*2;
end

function solve()
    if prob.time_dependent
        #### TEMP this is just for the heat equation, must change
        rhs = heatRHS;
        ###########
        
        stepper = init_stepper(allnodes, config.stepper);
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        T = @elapsed(stepper.step(stepper.Nsteps, stepper.dt, rhs));
        log_entry("Completed time steps in "*string(T)*" seconds.");
    else
        
    end
end

##### TEMP ##### Just for testing
function heatRHS(vars, time)
    # Define field differences at faces
    du = (vars[1].values[vmapm]-vars[1].values[vmapp])./2;
    du = reshape(du,(2,mesh_data.nel));
    
    # impose boundary condition -- Dirichlet BC's
    uin  = -vars[1].values[vmapI];
    uout = -vars[1].values[vmapO];
    du[mapI] = (vars[1].values[vmapI] - uin)./2;
    du[mapO] = (vars[1].values[vmapO] - uout)./2;
    
    # Compute q and form differences at faces
    q = rx.*(refel.Dr*vars[1].values) .- refel.lift*(Fscale.*(normals.*du));
    dq = (q[vmapm]-q[vmapp])./2;
    dq = reshape(dq,(2,mesh_data.nel));
    
    # impose boundary condition -- Neumann BC's
    qin  = q[vmapI];
    qout = q[vmapO];
    dq[mapI] = (q[vmapI] - qin)./2;
    dq[mapO] = (q[vmapO] - qout)./2;
    
    # compute right hand side
    return [rx.*(refel.Dr*q).-refel.lift*(Fscale.*(normals.*dq)) q];
end

end #module