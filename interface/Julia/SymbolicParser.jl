#=
# A set of tools for parsing the variational forms into expressions 
# with symbolic operators implemented by Femshop solvers.
=#
module SymbolicParser
export sp_parse, sp_term

# These lists need to be expanded
test_ops = [:grad, :div, :curl]; # operators that may be applied to the test function
var_ops = [:Dt, :grad, :div, :curl]; # operators that may be applied to variables or coefficients

# A single term from the variational form.
# Has an arity and an array of multiplied factors
# example: where v is a test function symbol and this is from an equation for u with variable symbols [u,r]
# 1*v               ->  arity=1, subexpressions = [1, nothing, v]
# 3*u*v             ->  arity=2, subexpressions = [3, u, v]
# grad(r)*grad(v)   ->  arity=1, subexpressions = [grad(r), nothing, grad(v)] (arity=1 because r is considered a coefficient in an equation for u)
# surface(u*v)      ->  arity=2, subexpressions = [nothing, surface(u), v] (surface is a special symbol for surface integrals)
struct sp_term
    arity::Int               # usually 1 or 2 for linear and bilinear terms
    factors::Array{Any,1}    # array [coefficient part, variable part,  test function part]
end

# include some utility functions
include("parser_utils.jl");

# Parses a variational form into symbolic expressions in terms of operators 
# implemented by Femshop solvers.
# Takes an Expr, variable symbol, and test function symbol
# Returns a tuple of expressions (lhs, rhs) for the equation lhs = rhs
# lhs contains terms with arity > 1 (bilinear+)
# rhs contains terms with arity =0,1 (linear or functionals)
function sp_parse(ex, var, test)
    terms = get_all_terms(expand(ex));
    spterms = Array{sp_term,1}(undef,length(terms));
    dtterms = []; # Will hold index of terms that are wrapped in Dt()
    for i=1:length(terms)
        # split the test, var and coef factors
        # Note: The number and symbol cases shouldn't happen, but are included for completeness
        if typeof(terms[i]) <: Number
            spterms[i] = sp_term(0, [terms[i], nothing, nothing]);
        elseif typeof(terms[i]) == Symbol
            if terms[i] === test
                spterms[i] = sp_term(1, [nothing, nothing, terms[i]]);
            elseif terms[i] === var
                spterms[i] = sp_term(2, [nothing, terms[i], nothing]);
            else
                spterms[i] = sp_term(0, [terms[i], nothing, nothing]);
            end
        elseif typeof(terms[i]) == Expr
            facs = get_all_factors(terms[i]);
            tf = []; # test function factors
            cf = []; # coefficient factors
            vf = []; # variable factors
            
            # handle Dt factors Dt(a*u*v) -> a*Dt(u*v)
            newfacs = copy(facs);
            refactor = false;
            for j=1:length(facs)
                # get past possible negative sign
                if typeof(facs[j]) == Expr && facs[j].args[1] === :-
                    if typeof(facs[j].args[2]) == Expr && facs[j].args[2].args[1] === :Dt
                        newfacs[j].args[2] = facs[j].args[2].args[2]; # wrap the whole term in Dt
                        refactor = true;
                    end
                else
                    if typeof(facs[j]) == Expr && facs[j].args[1] === :Dt
                        newfacs[j] = facs[j]. args[2]; # wrap the whole term in Dt
                        refactor = true;
                    end
                end
            end
            if refactor
                # Dt() is stripped off, but noted
                facs = get_all_factors(assemble_term(newfacs));
                push!(dtterms, i);
            end
            
            for j=1:length(facs)
                if is_test_op(facs[j], test)
                    tf = [tf ; facs[j]];
                elseif is_var_op(facs[j], var)
                    vf = [vf ; facs[j]];
                else
                    cf = [cf ; facs[j]];
                end
            end
            a = 0; # arity
            mulex = :(a*b);
            if length(tf) > 0
                a = 1;
                if length(tf) > 1
                    mulex.args = [:* ; tf];
                    tf = copy(mulex);
                else
                    tf = tf[1];
                end
            else
                tf = nothing;
            end
            if length(vf) > 0
                a = 2;
                if length(vf) > 1
                    mulex.args = [:* ; vf];
                    vf = copy(mulex);
                else
                    vf = vf[1];
                end
            else
                vf = nothing;
            end
            if length(cf) > 1
                mulex.args = [:* ; cf];
                cf = copy(mulex);
            elseif length(cf) == 1
                cf = cf[1];
            else
                cf = nothing;
            end
            
            spterms[i] = sp_term(a, [cf, vf, tf]);
        end
    end
    
    # Translate each term into expressions like p1*operator(p2)
    lhsterms = [];
    lhsdtterms = [];
    rhsterms = [];
    for i=1:length(spterms)
        adddt = false
        for j=1:length(dtterms)
            if dtterms[j] == i
                adddt = true;
            end
        end
        if spterms[i].arity > 1
            if adddt
                lhsdtterms = [lhsdtterms ; translate_term(spterms[i], var, test)];
            else
                lhsterms = [lhsterms ; translate_term(spterms[i], var, test)];
            end
        else
            rhsterms = [rhsterms ; distribute_negative(translate_term(spterms[i], var, test))];
        end
    end
    
    if length(lhsterms) > 1
        lhs = :(a+b);
        lhs.args = [:+ ; lhsterms];
    elseif length(lhsterms) == 1
        lhs = lhsterms[1];
    else
        lhs = 0; # this should never happen
    end
    
    if length(lhsdtterms) > 1
        lhsdt = :(a+b);
        lhsdt.args = [:+ ; lhsdtterms];
    elseif length(lhsdtterms) == 1
        lhsdt = lhsdtterms[1];
    else
        lhsdt = 0;
    end
    
    if length(rhsterms) > 1
        rhs = :(a+b);
        rhs.args = [:+ ; rhsterms];
    elseif length(rhsterms) == 1
        rhs = rhsterms[1];
    else
        rhs = 0;
    end
    
    if length(lhsdtterms) > 0
        return ((lhsdt, lhs), rhs);
    else
        return (lhs, rhs);
    end
end

# Translate a term into an expression using symbolic operators
# example with test function v and variable u:
# u*v  ->  mass(u)
# grad(u)*v  ->  advection(u)
# u*grad(v)  ->  advection_transpose(u)
# grad(u)*grad(v)  ->  stiffness(u)
function translate_term(t, var, test)
    facs = [];
    term = nothing;
    borl = :RHS;
    if t.arity == 0
        # If arity=0, there's nothing to translate.
        return t.factors[1];
    elseif t.arity == 2
        borl = :LHS;
    end
    
    # if there's a negative, strip it off and put it around everything -a*b -> -(a*b)
    negative = false;
    vpart = t.factors[2];
    if typeof(vpart) == Expr && vpart.args[1] === :-
        negative = !negative;
        vpart = t.factors[2].args[2];
    end
    tpart = t.factors[3];
    if typeof(tpart) == Expr && tpart.args[1] === :-
        negative = !negative;
        tpart = t.factors[3].args[2];
    end
    cpart = t.factors[1];
    if typeof(cpart) == Expr && cpart.args[1] === :-
        negative = !negative;
        cpart = cpart.args[2];
    end
    
    # What parts are present? Split them into p1*op(p2)
    p1 = nothing;
    p2 = nothing;
    op = nothing;
    if !(tpart === nothing) # There is a test function part
        if typeof(cpart) <: Number
            p1 = cpart;
            p2 = vpart;
        else
            p2 = assemble_term([cpart ; vpart]);
            if p2 == 0
                p2 = nothing;
            else
                (p1, p2) = extract_constant(p2);
            end
        end
        if typeof(tpart) == Symbol
            # (?,t)
            if typeof(p2) == Expr && p2.args[1] === :grad
                op = :advection_operator;
                p2 = p2.args[2];
            elseif !(p2 === nothing)
                op = :mass_operator;
            else
                op = :test_int
                p2 = tpart;
            end
        elseif typeof(tpart) == Expr && tpart.args[1] === :grad
            # (?, grad(t))
            if typeof(p2) == Expr && p2.args[1] === :grad
                op = :stiffness_operator;
                p2 = p2.args[2];
            elseif !(p2 === nothing)
                op = :advection_transpose_operator;
            else
                op = :test_int
                p2 = tpart;
            end
        else
            op = :unknown_operator;
        end
    else # There is no test function part?? This shouldn't happen
        p1 = assemble_term([cpart, vpart]);
    end
    
    if p1 === nothing
        if op === nothing || p2 === nothing
            return nothing;
        end
        term = :(a(b));
        #term.args[1] = op;
        #term.args[2] = p2;
        #term.args = [op ; p2 ; :x ; :refel ; borl];
        term.args = [op ; p2 ; :args];
    else
        if op === nothing || p2 === nothing
            term = p1;
        else
            term = :(a*b(c));
            term.args[2] = p1;
            #term.args[3].args[1] = op;
            #term.args[3].args[2] = p2;
            #term.args[3].args = [op ; p2 ; :x ; :refel ; borl];
            term.args[3].args = [op ; p2 ; :args];
        end
    end
    
    # negative sign
    if negative
        negex = :(-a);
        negex.args[2] = term;
        term = negex;
    end
    
    return term;
end

##### not ready ########

# # Turns something like surface(B*u*v) into something like (operator_surface_int(B*flux[:,:,u.index]))
# function replace_surface_terms(ex)
#     if !(typeof(ex) == Expr)
#         return ex;
#     elseif ex.head === :call
#         if ex.args[1] === :surface
#             # The expr in args[2] should be something like (B*u*v) for variable u and trial function v
#             st = surface_term(ex.args[2]);
#             ex = :(operator_surface_int(a));
#             ex.args[2] = st;
#         else
#             for i=2:length(ex.args)
#                 ex.args[i] = replace_surface_terms(ex.args[i]);
#             end
#         end
#     end
    
#     return ex;
# end

# # This gets the ex in surface(ex)
# # remove the *v part
# # turn var into flux[:,:,var.index]
# function surface_term(ex)
#     if typeof(ex) == Symbol
#         # is it a var or a test function?
#         for vi = 1:var_count
#             if ex === variables[vi].symbol
#                 return :(flux[:,:,$ex.index]);
#             end
#         end
            
#     elseif typeof(ex) == Expr && ex.head === :call
#         # work recursively on all sub expressions
#         for i=2:length(ex.args)
#             ex.args[i] = surface_term(ex.args[i]);
#         end
#         # if one of the args is the test function, remove it
#         # (a*b*v) -> (a*b) , (a*v) -> (a)
#         for i=2:length(ex.args)
#             ex.args[i] = surface_term(ex.args[i]);
#             if ex.args[i] === test_function_symbol
#                 # if (a*b*v)
#                 if length(ex.args) > 3
#                     newargs = copy(ex.args);
#                     newind = 1;
#                     for j=1:length(ex.args)
#                         if j!=i
#                             newargs[newind] = ex.args[j];
#                             newind += 1;
#                         end
#                     end
#                     ex.args = newargs;
#                     return ex;
#                 elseif i == 2
#                     return ex.args[3];
#                 else
#                     return ex.args[2];
#                 end
#             end 
#         end
#     end
    
#     return ex;
# end

end # module