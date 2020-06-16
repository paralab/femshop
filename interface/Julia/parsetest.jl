# test the parser
if !@isdefined(SymbolicParser)
    include("SymbolicParser.jl")
    using .SymbolicParser
end

ex = [
    :(0),
    :(5),
    :(a),
    :(u),
    :(v),
    :(2*a),
    :(2*u),
    :(2*v),
    :(a*u),
    :(a*v),
    :(a*b),
    :(u*v),
    :(u*u),
    :(0 + 1),
    :(0 + a),
    :(1 + a),
    :(1 + u),
    :(a + u),
    :(5 + a + u),
    :(5 + 7 + a + b + u + v),
    :(5*a + b),
    :(5*a + 7*u),
    :(5*a*u + 7*b*v),
    :(grad(2)),
    :(grad(2) + 0 + 0*b + a*grad(2) + a*grad(2)*grad(u)),
    :(grad(a)),
    :(grad(u)),
    :(grad(v)),
    :(grad(-u)),
    :(grad(-(a+b))),
    :(-grad(-2*(-u))),
    :(grad(2*u)),
    :(grad(a*u)),
    :(grad(1+2)),
    :(grad(2+a)),
    :(grad(a+u)),
    :(grad(a*u + 3*v)),
    :(grad(a*b - c*d)),
    :((1-r*u)*v + (a+b)*(c-d)*v - 5*grad(v) + 6 + u - v + 7*grad(3*u - a*v)),
    :(u*v),
    :(u*grad(v)),
    :(grad(u)*v),
    :(grad(u)*grad(v)),
    :(Dt(u*v) - grad(u)*v + u*grad(v) - grad(u)*grad(v)),
    :(Dt(u*v) + q*grad(v) + surface(flux(q)*v)),
    :(q*v + u*grad(v) + surface(flux(u)*v)),
    :(Dt(u*v) + a*q*grad(v) - grad(b*q)*v + 7*grad(u)*grad(v) + surface(flux(q)*v))
];

for i=1:length(ex)
    println("input:    "*string(ex[i]));
    println("expanded: "*string(SymbolicParser.expand(ex[i])));
    (lhs, rhs) = sp_parse(ex[i], :u, :v);
    println("lhs: "*string(lhs));
    println("rhs: "*string(rhs));
    println("");
end
