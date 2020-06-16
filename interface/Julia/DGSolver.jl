#=
# DG solver
=#
module DGSolver

export init_dgsolver, solve

import ..Femshop: JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET, NEUMANN, ROBIN,
        MSH_V2, MSH_V4,
        SCALAR, VECTOR, TENSOR
import ..Femshop: log_entry, printerr
import ..Femshop: config, prob, variables, mesh_data

using LinearAlgebra

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
# vmapI = [];
# vmapO = [];
# mapI = [];
# mapO = [];

include("refel.jl");
include("time_steppers.jl");
include("dg_operators.jl");

function init_dgsolver()
    if config.dimension == 1
        setup1D();
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
    global Fmask = [1;refel.Np]; # elemental nodes on the faces
    global normals = mesh_data.normals[:,1:2,1]';
    global Fscale = 1 ./J[Fmask,:];
    global e2e = mesh_data.neighbors[:,1:2];
    global e2f = ones(Int, nel, refel.Nfaces); # assumes left face is face 1, right is 2 (OK)
    e2f[2:end,1] = 2 .*e2f[2:end,1];
    e2f[end,2] = 2;
    global vmapm = zeros(Int, nel*refel.Nfaces);
    global vmapp = zeros(Int, nel*refel.Nfaces);
    for i=1:nel
        vmapm[i*2-1] = (i-1)*refel.Np + 1;
        vmapm[i*2] = i*refel.Np;
        
        j = mesh_data.neighbors[i,1];
        vmapp[i*2-1] = j*refel.Np;
        j = mesh_data.neighbors[i,2];
        vmapp[i*2] = (j-1)*refel.Np + 1;
    end
    # boundaries have their own index
    left = 1;
    right = nel;
    for i=1:nel
        if mesh_data.neighbors[i,1] == i
            left = i;
        end
        if mesh_data.neighbors[i,2] == i
            right = i;
        end
    end
    vl = (left-1)*refel.Nfaces + 1;
    vr = right*refel.Nfaces;
    vmapp[vl] = vmapm[vl];
    vmapp[vr] = vmapm[vr];
    global mapb = [vl; vr];
    global vmapb = vmapm[mapb];
    # global vmapI = vmapm[vl];
    # global vmapO = vmapm[vr];
    # global mapI = (left-1)*2 + 1;
    # global mapO = right*2;
end

function solve(lhs, lhspars, rhs, rhspars)
    if prob.time_dependent
        stepper = init_stepper(allnodes, config.stepper);
        log_entry("Beginning "*string(stepper.Nsteps)*" time steps.");
        T = @elapsed(stepper.step(stepper.Nsteps, stepper.dt, rhspars, rhs));
        log_entry("Completed time steps in "*string(T)*" seconds.");
    else
        
    end
end

function rhs_dg_1d(pars, v, vars, time)
    flux = zeros(2, mesh_data.nel, length(vars));
    rhsv = zeros(size(v));
    for vi=1:length(vars)
        if vars[vi].ready
            rhsv[:,:,vi] = v[:,:,vi];
        end
    end
    
    for i=1:length(pars.solve_order)
        next = pars.solve_order[i];
        for vi = 1:length(vars)
            if pars.dependence[next, vi] == 1
                subflux = (v[:,:,vi][vmapm]-v[:,:,vi][vmapp]).*0.5;
                subflux = reshape(subflux,(2,mesh_data.nel));
                # Impose boundary condition on flux
                apply_bc1D!(pars.bdry_type[vi], pars.bdry_func[vi].func, time, v[:,:,vi], subflux);
                flux[:,:,vi] = subflux;
            end
        end
        # evaluate the rhs expression
        rhsv[:,:,next] = pars.rhs_eq[next].func(v, flux, time);
        # If the LHS was like DT(u), change nothing. If it was like (u), put rhsv into v
        if !prob.lhs_time_deriv[next]
            v[:,:,next] = rhsv[:,:,next];
        end
    end
    
    return rhsv;
end

function apply_bc1D!(type, func, time, v, flux)
    x = allnodes[vmapb];
    if type == DIRICHLET
        flux[mapb] = v[vmapb] .+ func(x,time);
    elseif type == NEUMANN
        flux[mapb] = func(x,time) .* ones(size(flux[mapb]));
    elseif type == ROBIN
        # TODO
    end
end

end #module