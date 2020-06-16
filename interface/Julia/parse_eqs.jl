#=
# Parse LHS and RHS expressions to:
# 1. Extract variable dependencies
# 2. Construct equations in terms of symbolic operators
=#

function get_lhs_rhs(ex)
    
end

# Finds variable symbols to determine dependencies
function find_vars(ex)
    vars = [];
    if typeof(ex) == Symbol
        for v in variables
            if ex === v.symbol
                vars = [vars; v];
            end
        end
    elseif ex.head === :call
        for i=2:length(ex.args)
            if typeof(ex.args[i]) == Expr
                vars = [vars; find_vars(ex.args[i])];
            elseif typeof(ex.args[i]) == Symbol
                for v in variables
                    if ex.args[i] === v.symbol
                        vars = [vars; v];
                    end
                end
            end
        end
    end
    return vars;
end
