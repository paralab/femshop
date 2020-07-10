#=
# Operators used in generated expressions
=#
export mass_operator, stiffness_operator

function zero_operator(ex, args)
    solver.zero_operator(ex,args);
end

function one_operator(ex,args)
    solver.one_operator(ex,args);
end

function mass_operator(ex, args)
    solver.mass_operator(ex,args);
end

function stiffness_operator(ex, args)
    solver.stiffness_operator(ex,args);
end

