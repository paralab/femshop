#=
Define some constants
Use string values to make printing/interpretation easier
=#
export JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING, EULER_EXPLICIT, EULER_IMPLICIT, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT, DIRICHLET

# Languages for generated code
const JULIA = "Julia";
const CPP = "C++";
const MATLAB = "Matlab";

# Domain types
const SQUARE = "square";
const IRREGULAR = "irregular";

# Domain decomposition
const TREE = "tree";
const UNSTRUCTURED = "unstructured";

# Solver type
const CG = "CG";
const DG = "DG";
const HDG = "HDG";

const NODAL = "nodal";
const MODAL = "modal";

# Function space
const LEGENDRE = "Legendre";

# Element node positions
const UNIFORM = "uniform";
const GAUSS = "Gauss";
const LOBATTO = "Lobatto";

# Nonlinear solver methods
const NONLINEAR_NEWTON = "Newton";
const NONLINEAR_SOMETHING = "something";

# Time steppers
const EULER_EXPLICIT = "Euler-explicit";
const EULER_IMPLICIT = "Euler-implicit";
const RK4 = "RK4";
const LSRK4 = "LSRK4";
const ABM4 = "ABM4";

# Linear system solvers/structures
const OURS = "ours";
const PETSC = "PETSC";

# Output format
const VTK = "vtk";
const RAW_OUTPUT = "raw";
const CUSTOM_OUTPUT = "custom";

#BC
const DIRICHLET = "Dirichlet";
