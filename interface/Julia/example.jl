# Just a scratch file for testing things

if !@isdefined(Femshop)
    include("Femshop.jl");
    using .Femshop
end
# Try making an optional log
@useLog("logfile")

# Set up the configuration (order doesn't matter)
@domain(2, IRREGULAR, UNSTRUCTURED) # dimension, geometry, decomposition
@mesh("circle_2d.msh")              # .msh file or generate our own
@solver(DG)                         # DG, CG, etc.
@functionSpace(LEGENDRE, 4)         # function, order (or use testFunction and trialFunction)
@nodes(LOBATTO)                     # elemental node arrangement

# check
c = Femshop.config
log_dump_config(c);
@assert(c.geometry == IRREGULAR)
@assert(Femshop.mesh_data.nx > 0)

# Try making some files
@language(CPP, "tryitcpp", "Just an empty file for testing");
@finalize