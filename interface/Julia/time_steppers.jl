#=
# A set of time steppers
=#

mutable struct Stepper
    type;       # The constant for the stepper type
    step;       # a function handle to the appropriate stepper
    Nsteps;     # number of steps
    dt;         # step size
    cfl;        # CFL number
    
    Stepper(t, c) = new(t, euler_implicit, 0, 0, c);
end

function init_stepper(x, stepper)
    if stepper.type == EULER_EXPLICIT
        dxmin = abs(x[2,1]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.step = euler_explicit;
        
        return stepper;
        
    elseif stepper.type == EULER_IMPLICIT
        dxmin = abs(x[2,1]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.step = euler_implicit;
        
        return stepper;
        
    elseif stepper.type == CRANK_NICHOLSON
        dxmin = abs(x[2,1]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.step = crank_nicholson;
        
        return stepper;
        
    elseif stepper.type == LSRK4
        dxmin = abs(x[2,1]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.step = lsrk4;
        
        return stepper;
    end
end

################## TODO - possibly not needed - maybe built into the generated code
function euler_explicit()
    
end

function euler_implicit()
    
end

function crank_nicholson()
    
end

function rk4()
    
end

function lsrk4(Nsteps, dt, rhsparams, rhs)
    # rk4a = [0.0; -0.41789047449985195; -1.192151694642677; -1.6977846924715279; -1.5141834442571558];
    # rk4b = [0.14965902199922912; 0.37921031299962726;  0.8229550293869817; 0.6994504559491221; 0.15305724796815198];
    # rk4c = [0.0; 0.14965902199922912; 0.37040095736420475; 0.6222557631344432; 0.9582821306746903];
    
    # resu = zeros(size(variables[1].values)[1], size(variables[1].values)[2], length(variables)); #TODO only enough for scalars
    # rhsv = zeros(size(variables[1].values)[1], size(variables[1].values)[2], length(variables)); #TODO only enough for scalars
    # for vi=1:length(variables)
    #     rhsv[:,:,vi] = variables[vi].values;
    # end
    # time = 0;
    # rktime = 0;
    # for ti=1:Nsteps
    #     for rki=1:5
    #         rktime = time + rk4c[rki]*dt;
    #         rhsv = rhs(rhsparams, rhsv, variables, rktime);
    #         resu = rk4a[rki].*resu + dt.*rhsv;
    #         for vi=1:length(variables)
    #             variables[vi].values = variables[vi].values .+ rk4b[rki].*resu[:,:,vi];
    #             rhsv[:,:,vi] = variables[vi].values;
    #         end
    #     end
    #     time += dt;
    # end
end
    