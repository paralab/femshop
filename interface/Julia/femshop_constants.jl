#=
Define some constants
=#
export JULIA, CPP, MATLAB, SQUARE, IRREGULAR, TREE, UNSTRUCTURED, CG, DG, HDG,
        NODAL, MODAL, LEGENDRE, UNIFORM, GAUSS, LOBATTO, NONLINEAR_NEWTON,
        NONLINEAR_SOMETHING_ELSE, EULER_EXPLICIT, EULER_IMPLICIT, RK4, LSRK4,
        ABM4, OURS, PETSC, VTK, RAW_OUTPUT, CUSTOM_OUTPUT

# Languages for generated code
const JULIA = 1;
const CPP = 2;
const MATLAB = 3;

# Domain types
const SQUARE = 11;
const IRREGULAR = 12;

# Domain decomposition
const TREE = 21;
const UNSTRUCTURED = 22;

# Solver type
const CG = 31;
const DG = 32;
const HDG = 33;

const NODAL = 34;
const MODAL = 35;

# Function space
const LEGENDRE = 41;

# Element node positions
const UNIFORM = 51;
const GAUSS = 52;
const LOBATTO = 53;

# Nonlinear solver methods
const NONLINEAR_NEWTON = 61;
const NONLINEAR_SOMETHING_ELSE = 62;

# Time steppers
const EULER_EXPLICIT = 71;
const EULER_IMPLICIT = 72;
const RK4 = 73;
const LSRK4 = 74;
const ABM4 = 75;

# Linear system solvers/structures
const OURS = 81;
const PETSC = 82;

# Output format
const VTK = 91;
const RAW_OUTPUT = 92;
const CUSTOM_OUTPUT = 93;
