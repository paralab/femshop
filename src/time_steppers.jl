#=
# A set of time steppers
=#

mutable struct Stepper
    type;       # The constant for the stepper type
    Nsteps;     # number of steps
    dt;         # step size
    cfl;        # CFL number
    stages;     # how many stages
    a;          # for RK steppers
    b;          # for RK steppers
    c;          # for RK steppers
    
    Stepper(t, c) = new(t, 0, 0, c, 0, [], [], []);
end

function init_stepper(x, stepper)
    if stepper.type == EULER_EXPLICIT
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == EULER_IMPLICIT
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = (prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == CRANK_NICHOLSON
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 1;
        
        return stepper;
        
    elseif stepper.type == LSRK4
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 5;
        stepper.a = [0; -0.41789047449985195; -1.192151694642677; -1.6977846924715279; -1.5141834442571558];
        stepper.b = [0.14965902199922912; 0.37921031299962726;  0.8229550293869817; 0.6994504559491221; 0.15305724796815198];
        stepper.c = [0.0; 0.14965902199922912; 0.37040095736420475; 0.6222557631344432; 0.9582821306746903];
        
        return stepper;
        
    elseif stepper.type == RK4
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        stepper.stages = 4;
        stepper.a = [   0   0   0   0;
                        0.5 0   0   0;
                        0   0.5 0   0;
                        0   0   1   0];
        stepper.b = [1/6; 1/3; 1/3; 1/6];
        stepper.c = [0.0; 0.5; 0.5; 1];
        
        return stepper;
    end
end

function reformat_for_stepper(lhs, rhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if length(rhs)>1 #multi dof
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], stepper);
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
            end
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = lhs[2][1][i]*dt; # dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            
            newlhs = copy(lhs[1][1]);
            append!(newlhs, lhs[2][1]); # lhs1 + dt*lhs2
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[1][1]);# dt*rhs + lhs1
            
        elseif stepper == EULER_EXPLICIT # lhs1 = dt*rhs - dt*lhs2 + lhs1
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]*dt; # -dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# dt*rhs - dt*lhs2
            append!(newrhs, lhs[1][1]);# dt*rhs - dt*lhs2 + lhs1
            
        elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
            lhs2l = copy(lhs[2][1]);
            lhs2r = copy(lhs[2][1]);
            for i=1:length(lhs[2][1])
                lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
                lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            
            newlhs = copy(lhs[1][1]);
            append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[1][1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            append!(newrhs, lhs2r);
            
        elseif stepper == LSRK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]; # -lhs2
            end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]; # rhs
            # end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# rhs - lhs2
            
        elseif stepper == RK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]; # -lhs2
            end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]; # rhs
            # end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# rhs - lhs2
            
        end
    end
    
    return (newlhs, newrhs);
end

function reformat_for_stepper(lhs, rhs, face_lhs, face_rhs,stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    newfacelhs = [];
    newfacerhs = [];
    
    if length(rhs)>1 #multi dof
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        newfacelhs = copy(face_rhs);
        newfacerhs = copy(face_rhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi], newfacelhs[vi], newfacerhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], face_lhs[vi], face_rhs[vi], stepper);
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
                newfacelhs[vi] = face_lhs[vi];
                newfacerhs[vi] = face_rhs[vi];
            end
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = lhs[2][1][i]*dt; # dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs[1])
                face_rhs[1][i] = face_rhs[1][i]*dt; # dt*facelhs
            end
            for i=1:length(face_lhs[1])
                face_lhs[1][i] = face_lhs[1][i]*dt; # dt*facerhs
            end
            
            newlhs = copy(lhs[1][1]);
            append!(newlhs, lhs[2][1]); # lhs1 + dt*lhs2
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[1][1]);# dt*rhs + lhs1
            newfacerhs = copy(face_rhs[1]);# dt*facerhs
            newfacelhs = copy(face_lhs[1]);# dt*facelhs
            
        elseif stepper == EULER_EXPLICIT # lhs1 = dt*rhs - dt*lhs2
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]*dt; # -dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs[1])
                face_rhs[1][i] = face_rhs[1][i]*dt; # dt*facerhs
            end
            for i=1:length(face_lhs[1])
                face_lhs[1][i] = -face_lhs[1][i]*dt; # -dt*facelhs
            end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# dt*rhs - dt*lhs2 + lhs1
            #append!(newrhs, lhs[1][1]);
            newfacerhs = copy(face_rhs[1]);
            append!(newfacerhs, face_lhs[1]);# dt*facerhs - dt*facelhs
            newfacelhs = [];
            
        elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
            lhs2l = copy(lhs[2][1]);
            lhs2r = copy(lhs[2][1]);
            facelhs2l = copy(face_lhs[1]);
            facelhs2r = copy(face_lhs[1]);
            for i=1:length(lhs[2][1])
                lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
                lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(face_lhs[1])
                facelhs2l[i] = Basic(0.5)*facelhs2l[i]*dt; # 0.5*dt*lhs2
                facelhs2r[i] = Basic(-0.5)*facelhs2r[i]*dt; # -0.5*dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs[1])
                face_rhs[1][i] = face_rhs[1][i]*dt; # dt*facelhs
            end
            
            newlhs = copy(lhs[1][1]);
            append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[1][1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            append!(newrhs, lhs2r);
            newfacelhs = facelhs2l;
            newfacerhs = face_rhs[1];
            append!(newfacerhs, facelhs2r);
            
        elseif stepper == LSRK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]; # -lhs2
            end
            for i=1:length(face_lhs[1])
                face_lhs[1][i] = -face_lhs[1][i]; # -facelhs
            end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# rhs - lhs2
            newfacerhs = copy(face_rhs[1]);
            append!(newfacerhs, face_lhs[1]);# facerhs - facelhs
            newfacelhs = [];
            
        elseif stepper == RK4 # (lhs1) : rhs - lhs2
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]; # -lhs2
            end
            for i=1:length(face_lhs[1])
                face_lhs[1][i] = -face_lhs[1][i]; # -facelhs
            end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# rhs - lhs2
            newfacerhs = copy(face_rhs[1]);
            append!(newfacerhs, face_lhs[1]);# facerhs - facelhs
            newfacelhs = [];
            
        end
    end
    
    return (newlhs, newrhs, newfacelhs, newfacerhs);
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv(dtvar, flhs, frhs, slhs, srhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newflhs = []; # flux
    newfrhs = [];
    newslhs = []; # source
    newsrhs = [];
    
    if length(dtvar)>1 #multi dof
        newflhs = copy(frhs);
        newfrhs = copy(frhs);
        newslhs = copy(srhs);
        newsrhs = copy(srhs);
        for vi=1:length(dtvar);
            (newflhs[vi], newfrhs[vi], newslhs[vi], newsrhs[vi]) = reformat_for_stepper_fv(dtvar[vi], flhs[vi], frhs[vi], slhs[vi], srhs[vi], stepper);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_IMPLICIT # lhs1 + dt*lhs2 = dt*rhs + lhs1
            # for i=1:length(lhs[2][1])
            #     lhs[2][1][i] = lhs[2][1][i]*dt; # dt*lhs2
            # end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            # end
            
            # newlhs = copy(lhs[1][1]);
            # append!(newlhs, lhs[2][1]); # lhs1 + dt*lhs2
            # newrhs = copy(rhs[1]);
            # append!(newrhs, lhs[1][1]);# dt*rhs + lhs1
            
        elseif stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4
            # (lhs1 - lhs1)/dt = rhs - lhs2
            for i=1:length(flhs[1])
                flhs[1][i] = -flhs[1][i]; # lhs2
            end
            # for i=1:length(slhs[1])
            #     slhs[1][i] = slhs[1][i]; # -lhs2
            # end
            
            newflhs = [];# lhs1
            newfrhs = copy(frhs[1]);
            append!(newfrhs, flhs[1]);# rhs - lhs2
            
            newslhs = [];# lhs1
            newsrhs = copy(srhs[1]);
            append!(newsrhs, slhs[1]);# rhs - lhs2
            
        elseif stepper == CRANK_NICHOLSON # lhs1 + 0.5*dt*lhs2 = dt*rhs - 0.5*dt*lhs2 + lhs1
            # lhs2l = copy(lhs[2][1]);
            # lhs2r = copy(lhs[2][1]);
            # for i=1:length(lhs[2][1])
            #     lhs2l[i] = Basic(0.5)*lhs2l[i]*dt; # 0.5*dt*lhs2
            #     lhs2r[i] = Basic(-0.5)*lhs2r[i]*dt; # -0.5*dt*lhs2
            # end
            # for i=1:length(rhs[1])
            #     rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            # end
            
            # newlhs = copy(lhs[1][1]);
            # append!(newlhs, lhs2l); # lhs1 + 0.5*dt*lhs2
            # newrhs = copy(rhs[1]);
            # append!(newrhs, lhs[1][1]);# dt*rhs - 0.5*dt*lhs2 + lhs1
            # append!(newrhs, lhs2r);
            
        end
    end
    
    return (newflhs, newfrhs, newslhs, newsrhs);
end

# Special version for FV. Flux term only
function reformat_for_stepper_fv_flux(flhs, frhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newflhs = []; # flux
    newfrhs = [];
    
    if length(frhs)>1 #multi dof
        newflhs = copy(frhs);
        newfrhs = copy(frhs);
        for vi=1:length(frhs);
            (newflhs[vi], newfrhs[vi]) = reformat_for_stepper_fv_flux(flhs[vi], frhs[vi], stepper);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4
            # du/dt = -(flhs + frhs)
            # note flhs and frhs do not mean lhs and rhs at this point
            for i=1:length(flhs[1])
                flhs[1][i] = -flhs[1][i];
            end
            # for i=1:length(frhs)
            #     frhs[i] = -frhs[i]; # don't change sign for rhs, because already done by parser
            # end
            
            newflhs = []; # put everything in rhs for explicit steppers
            newfrhs = copy(frhs[1]);
            append!(newfrhs, flhs[1]);
            
        elseif stepper == EULER_IMPLICIT
            
        elseif stepper == CRANK_NICHOLSON 
            
        end
    end
    
    return (newflhs, newfrhs);
end

# Special version for FV. Assumes a Dt(u) term that is not explicitly included.
function reformat_for_stepper_fv_source(slhs, srhs, stepper)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newslhs = []; # source
    newsrhs = [];
    
    if length(srhs)>1 #multi dof
        newslhs = copy(srhs);
        newsrhs = copy(srhs);
        for vi=1:length(srhs);
            (newslhs[vi], newsrhs[vi]) = reformat_for_stepper_fv_source(slhs[vi], srhs[vi], stepper);
        end
    else
        # reformat depending on stepper type
        if stepper == EULER_EXPLICIT || stepper == RK4 || stepper == LSRK4
            # du/dt = slhs + srhs
            # note slhs and srhs do not mean lhs and rhs at this point
            for i=1:length(srhs)
                srhs[i] = -srhs[i]; # change sign for rhs, because done by parser
            end
            
            newslhs = [];# lhs is empty for explicit
            newsrhs = copy(srhs[1]);
            append!(newsrhs, slhs[1]); # rhs
            
        elseif stepper == EULER_IMPLICIT
            
        elseif stepper == CRANK_NICHOLSON
            
        end
    end
    
    return (newslhs, newsrhs);
end