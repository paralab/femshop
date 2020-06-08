#=
# A set of time steppers
=#

mutable struct Stepper
    step;       # a function handle to the appropriate stepper
    Nsteps;     # number of steps
    dt;         # step size
end

function init_stepper(x, type)
    if type == LSRK4
        dxmin = abs(x[2,1]-x[1,1]); # TODO this only works for 1D similar elements
        CFL = 0.25;
        dt = CFL*dxmin*dxmin;
        Nsteps = ceil(prob.end_time/dt);
        dt = prob.end_time/Nsteps;
        
        return Stepper(lsrk4, Nsteps, dt);
    end
end

function euler_explicit()
    
end

function euler_implicit()
    
end

function rk4()
    
end

function lsrk4(Nsteps, dt, rhsparams, rhs)
    rk4a = [0.0; -0.41789047449985195; -1.192151694642677; -1.6977846924715279; -1.5141834442571558];
    rk4b = [0.14965902199922912; 0.37921031299962726;  0.8229550293869817; 0.6994504559491221; 0.15305724796815198];
    rk4c = [0.0; 0.14965902199922912; 0.37040095736420475; 0.6222557631344432; 0.9582821306746903];
    
    resu = zeros(size(variables[1].values)[1], size(variables[1].values)[2], length(variables)); #TODO only enough for scalars
    rhsv = zeros(size(variables[1].values)[1], size(variables[1].values)[2], length(variables)); #TODO only enough for scalars
    time = 0;
    rktime = 0;
    for ti=1:Nsteps
        for rki=1:5
            rktime = time + rk4c[rki]*dt;
            rhsv = rhs(rhsparams, rhsv, variables, rktime);
            resu = rk4a[rki].*resu + dt.*rhsv;
            for vi=1:length(variables)
                variables[vi].values = variables[vi].values .+ rk4b[rki].*resu[:,:,vi];
            end
        end
        time += dt;
    end
end
    