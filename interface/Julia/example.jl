# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end

@domain(1, IRREGULAR, UNSTRUCTURED);

@solver(DG)

@functionSpace(LEGENDRE)
@order(4)

c = Femshop.config
c.geometry == IRREGULAR

@language(CPP, "tryit", "Some header text");

@finalize
