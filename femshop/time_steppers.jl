#=
# A set of time steppers
=#

mutable struct Stepper
    type;       # The constant for the stepper type
    Nsteps;     # number of steps
    dt;         # step size
    cfl;        # CFL number
    
    Stepper(t, c) = new(t, 0, 0, c);
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
        
        return stepper;
        
    elseif stepper.type == EULER_IMPLICIT
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 1;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = (prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        
        return stepper;
        
    elseif stepper.type == CRANK_NICHOLSON
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        
        return stepper;
        
    elseif stepper.type == LSRK4
        dxmin = abs(x[1,2]-x[1,1]); # TODO this only works for similar, square elements
        if stepper.cfl == 0
            stepper.cfl = 0.25;
        end
        stepper.dt = stepper.cfl*dxmin*dxmin;
        stepper.Nsteps = ceil(prob.end_time/stepper.dt);
        stepper.dt = prob.end_time/stepper.Nsteps;
        
        return stepper;
    end
end

function reformat_for_stepper(lhs, rhs, stepper, wrap=true)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    
    if length(rhs)>1 #multi dof
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], stepper, false);
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
            append!(newrhs, lhs[2][1]);# dt*rhs - dt*lhs2 + lhs1
            append!(newrhs, lhs[1][1]);
            
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
            
        end
    end
    
    # wrap in an array if needed (to avoid doing it on recursive steps)
    # if wrap
    #     newlhs = [newlhs];
    #     newrhs = [newrhs];
    # end
    
    return (newlhs, newrhs);
end

function reformat_for_stepper(lhs, rhs, face_lhs, face_rhs,stepper, wrap=true)
    # rebuild the expressions depending on type of time stepper
    dt = symbols("dt");
    newlhs = [];
    newrhs = [];
    newfacelhs = [];
    newfacerhs = [];
    
    if length(rhs)>1 #multi dof
        newlhs = copy(rhs);
        newrhs = copy(rhs);
        newfacelhs = copy(facerhs);
        newfacerhs = copy(facerhs);
        for vi=1:length(rhs)
            if length(lhs[1][vi]) > 0 # this dof has a time derivative term
                (newlhs[vi], newrhs[vi]) = reformat_for_stepper((lhs[1][vi], lhs[2][vi]), rhs[vi], face_lhs[vi], face_rhs[vi], stepper, false);
            else # no time derivative for this dof
                newlhs[vi] = lhs[2][vi];
                newrhs[vi] = rhs[vi];
                newfacelhs[vi] = facelhs[vi];
                newfacerhs[vi] = facerhs[vi];
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
            
        elseif stepper == EULER_EXPLICIT # lhs1 = dt*rhs - dt*lhs2 + lhs1
            for i=1:length(lhs[2][1])
                lhs[2][1][i] = -lhs[2][1][i]*dt; # -dt*lhs2
            end
            for i=1:length(rhs[1])
                rhs[1][i] = rhs[1][i]*dt; # dt*rhs
            end
            for i=1:length(face_rhs[1])
                face_rhs[1][i] = -face_rhs[1][i]*dt; # -dt*facelhs
            end
            for i=1:length(face_lhs[1])
                face_lhs[1][i] = -face_lhs[1][i]*dt; # -dt*facerhs
            end
            
            newlhs = copy(lhs[1][1]);# lhs1
            newrhs = copy(rhs[1]);
            append!(newrhs, lhs[2][1]);# dt*rhs - dt*lhs2 + lhs1
            append!(newrhs, lhs[1][1]);
            newfacerhs = copy(face_rhs[1]);
            append!(newfacerhs, face_lhs[1][1]);# dt*facerhs + dt*facelhs
            newfacelhs = [Basic(0)];
            
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
            
        end
    end
    
    # wrap in an array if needed (to avoid doing it on recursive steps)
    # if wrap
    #     newlhs = [newlhs];
    #     newrhs = [newrhs];
    # end
    
    return (newlhs, newrhs, newfacelhs, newfacerhs);
end
