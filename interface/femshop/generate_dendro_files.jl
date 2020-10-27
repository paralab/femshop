#=
# Dendro file generation functions
=#

# builds a c++ functional for a constant
function dendro_number_to_function(name, val)
    indent = "";
    args = ["x"; "y"; "z"; "var"];
    argtypes = ["double"; "double"; "double"; "double*"];
    ret = "";
    rettype = "void";
    captures = "gridX_to_X,gridY_to_Y,gridZ_to_Z";
    content = ["var[0] = "*string(val)*";"];
    return cpp_functional(indent, name, args, argtypes, ret, rettype, captures, content);
end
# changes symbol "a" to symbol "b" in expression ex
function swap_symbol(a, b, ex)
    if typeof(ex) == Symbol
        if ex === a
            return b;
        else
            return ex;
        end
    elseif typeof(ex) <: Number
        return ex;
    elseif typeof(ex) == Expr && length(ex.args) > 1
        swapped = copy(ex);
        for i=1:length(ex.args)
            swapped.args[i] = swap_symbol(a,b,ex.args[i]);
        end
        return swapped;
    else
        return ex;
    end
end

function dendro_genfunction_to_string(genfun)
    newex = swap_symbol(:x, :(gridX_to_X(x)), genfun.expr); # swap x for gridX_to_X(x)
    newex = swap_symbol(:y, :(gridY_to_Y(y)), newex); # swap y for gridY_to_Y(y)
    newex = swap_symbol(:z, :(gridZ_to_Z(z)), newex); # swap z for gridZ_to_Z(z)
    newex = swap_symbol(:pi, :M_PI, newex); # swap pi for M_PI
    s = string(newex);
    ns = replace(s, r"([\d)])([(A-Za-z])" => s"\1*\2"); # explicitly multiply with "*"
    return ns;
end

function dendro_config_file(dparams)
    file = genfiles.config;
    # dparams has (maxdepth, wavelet_tol, partition_tol, solve_tol, solve_max_iters)
    println(file, "m_uiMaxDepth = "*string(dparams[1])*";           // mesh refinement depth");
    println(file, "double wavelet_tol = "*string(dparams[2])*";     // tolerance for approximating functions(f) determines mesh fineness");
    println(file, "double partition_tol = "*string(dparams[3])*"; // load balancing parameter");
    println(file, "double solve_tol = "*string(dparams[4])*";            // tol used by cgsolve for stopping iterations");
    println(file, "unsigned int solve_max_iters = "*string(dparams[5])*"; // used by cgsolve");
    # From the config struct
    println(file, "unsigned int eOrder  = "*cpp_gen_string(config.basis_order_min)*";");
    
    println(file, "int config_dimension = "*cpp_gen_string(config.dimension)*";"); 
    
    println(file, "const char* config_solver = "*cpp_gen_string(config.solver_type)*";");
    println(file, "const char* config_trial_function = "*cpp_gen_string(config.trial_function)*";");
    println(file, "const char* config_test_function = "*cpp_gen_string(config.test_function)*";");
    println(file, "const char* config_elemental_nodes = "*cpp_gen_string(config.elemental_nodes)*";");
    println(file, "const char* config_quadrature= "*cpp_gen_string(config.quadrature)*";");
end

function dendro_genfunction_file()
    file = genfiles.genfunction;
    indent = "";
    args = ["x"; "y"; "z"; "var"];
    argtypes = ["double"; "double"; "double"; "double*"];
    rettype = "void";
    
    for i = 1:length(genfunctions)
        if isfunction(genfunctions[i])
            str = dendro_genfunction_to_string(genfunctions[i]);
            content = "var[0] = "*str*";";
            #lines = cpp_function_def(indent, genfunctions[i].name, args, argtypes, rettype, content);
            lines = ["std::function<void(double,double,double,double*)> "*genfunctions[i].name*" = [gridX_to_X, gridY_to_Y, gridZ_to_Z](const double x,const double y,const double z,double *var){", content, "};"];
            for j=1:length(lines)
                println(file, lines[j]);
            end
        end
    end
end

function dendro_prob_file()
    file = genfiles.problem;
    
    # variable and coefficient DOFs
    dofs = 0;
    varpart = "";
    varnames = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar.vals)
            varpart = varpart*"M_UI"*string(variables[i].symvar.vals[j]);
            varnames = varnames*"\"m_ui"*string(variables[i].symvar.vals[j])*"\", ";
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
            dofs = dofs+1;
        end
    end
    coefpart = "";
    coefnames = "";
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar.vals)
            coefpart = coefpart*"M_UI"*string(coefficients[i].symvar.vals[j])*", ";
            coefnames = coefnames*"\"m_ui"*string(coefficients[i].symvar.vals[j])*"\", ";
            dofs = dofs+1;
        end
    end
    dofs += 1; # for linear part
    println(file, "enum VAR{"*varpart*coefpart*"M_UI_RHS}; // variables, coefficients, linear part");
    println(file, "const char * VAR_NAMES[]={"*varnames*coefnames*"\"m_uiRHS\"};");
    println(file, "const unsigned int DOF="*string(dofs)*";"); 
    
    # From the prob struct
    println(file, "const char* prob_bc_type[] = "*cpp_gen_string(prob.bc_type)*";");
    println(file, "int prob_bid[] = "*cpp_gen_string(prob.bid)*";");
    
    # for time dependent problems
    if prob.time_dependent
        println(file, "double tBegin = 0;");
        println(file, "double tEnd = "*string(prob.end_time)*";");
        println(file, "double dt = "*string(time_stepper.dt)*";");
        println(file, "int tN = "*string(time_stepper.Nsteps)*";");
    end
end

function dendro_output_file()
    file = genfiles.output;
    if config.output_format == VTK
        println(file, "octDA->vecTopvtu(uSolVecPtr,\""*genFileName*"\",(char**)VAR_NAMES,false,false,DOF);");
    end
end

function dendro_bilinear_file(bl)
    file = genfiles.bilinear;
    
    print(file, bl);
end

function dendro_linear_file(l)
    file = genfiles.linear;
    
    print(file, l);
end

function dendro_stepper_file()
    file = genfiles.stepper;
end

function dendro_main_file()
    file = genfiles.main;
    # also write the two skeleton files linear_skel.cpp and bilinear_skel.cpp and their headers
    linfile = open(genDir*"/src/linear_skel.cpp", "w");
    bilinfile = open(genDir*"/src/bilinear_skel.cpp", "w");
    linheaderfile = open(genDir*"/include/linear_skel.h", "w");
    bilinheaderfile = open(genDir*"/include/bilinear_skel.h", "w");
    dendro_linear_skeleton(linfile, linheaderfile);
    dendro_bilinear_skeleton(bilinfile, bilinheaderfile);
    close(linfile);
    close(bilinfile);
    close(linheaderfile);
    close(bilinheaderfile);
    # And write a cmakelists file and README with simple instructions for compiling
    cmakefile = open(genDir*"/CMakeLists.txt", "w");
    readmefile = open(genDir*"/README.txt", "w");
    dendro_cmake_file(cmakefile);
    dendro_readme_file(readmefile);
    close(cmakefile);
    close(readmefile);
    
    varname = "variable_"* string(variables[1].symbol);
    valvec = string(variables[1].symvar.vals[1]);
    
    alloc_dof = "";
    set_dof = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar.vals)
            nam = string(variables[i].symvar.vals[j]);
            # if initial condition, use it, otherwise zero
            if length(prob.initial) >= i && !(prob.initial[i] === nothing)
                if typeof(prob.initial[i]) <: Array
                    fun = cpp_gen_string(prob.initial[i][j]);
                else
                    fun = cpp_gen_string(prob.initial[i]);
                end
            else
                fun = "zero_init"
            end
            alloc_dof = alloc_dof*"double * "*nam*"=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI"*nam*", false,false);\n";
            set_dof = set_dof*"octDA->setVectorByFunction("*nam*","*fun*",false,false,1);\n";
        end
    end
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar.vals)
            nam = string(coefficients[i].symvar.vals[j]);
            fun = cpp_gen_string(coefficients[i].value[j]);
            alloc_dof = alloc_dof*"double * "*nam*"=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI"*nam*", false,false);\n";
            set_dof = set_dof*"octDA->setVectorByFunction("*nam*","*fun*",false,false,1);\n";
        end
    end
    
    solvepart = "";
    if prob.time_dependent
        solvepart = "double currentT = 0.0;
        for(int ti=0; ti<tN; ti++){
            rhsVec.computeVec("*valvec*",rhs,1.0);
            lhsMat.cgSolve("*valvec*",rhs,solve_max_iters,solve_tol,0);
            currentT += dt;
        }";
    else
        solvepart = "// This uses the generated RHS code to compute the RHS vector.
        rhsVec.computeVec("*valvec*",rhs,1.0);
        
        // Solve the linear system. 
        lhsMat.cgSolve("*valvec*",rhs,solve_max_iters,solve_tol,0);";
    end
    
    # Just write the whole skeleton
    content = """
    #include "TreeNode.h"
    #include "mpi.h"
    #include "genPts_par.h"
    #include "sfcSort.h"
    #include "mesh.h"
    #include "dendro.h"
    #include "dendroIO.h"
    #include "octUtils.h"
    #include "functional"
    #include "fdCoefficient.h"
    #include "stencil.h"
    #include "rkTransport.h"
    #include "refel.h"
    #include "operators.h"
    #include "cg.h"
    
    #include "linear_skel.h"
    #include "bilinear_skel.h"

    int main (int argc, char** argv)
    {

        MPI_Init(&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;

        int rank, npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &npes);
        
        m_uiMaxDepth = 4; // a default value, but should be set in config.cpp
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Config.cpp"
        ////////////////////////////////////////////////////////////////////////////
        
        Point domain_min(0,0,0);
        Point domain_max(1,1,1);
        
        Point grid_min(0, 0, 0);
        Point grid_max((1u << m_uiMaxDepth), (1u << m_uiMaxDepth), (1u << m_uiMaxDepth));
        
        double Rg_x=(grid_max.x()-grid_min.x());
        double Rg_y=(grid_max.y()-grid_min.y());
        double Rg_z=(grid_max.z()-grid_min.z());

        double Rd_x=(domain_max.x()-domain_min.x());
        double Rd_y=(domain_max.y()-domain_min.y());
        double Rd_z=(domain_max.z()-domain_min.z());

        const Point d_min=domain_min;
        const Point d_max=domain_max;

        const Point g_min=grid_min;
        const Point g_max=grid_max;
        
        std::function<double(double)> gridX_to_X = [d_min,g_min,Rd_x,Rg_x](const double x){
            return d_min.x() + (x-g_min.x())*Rd_x/Rg_x;
        };
        
        std::function<double(double)> gridY_to_Y = [d_min,g_min,Rd_y,Rg_y](const double y){
            return d_min.y() + (y-g_min.y())*Rd_y/Rg_y;
        };
        
        std::function<double(double)> gridZ_to_Z = [d_min,g_min,Rd_z,Rg_z](const double z){
            return d_min.z() + (z-g_min.z())*Rd_z/Rg_z;
        };
        
        std::function<void(double,double,double,double*)> zero_init = [](const double x,const double y,const double z,double *var){
            var[0]=0;
        };
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Genfunction.cpp"
        #include "Problem.cpp"
        /////////////////////////////////////////////////////////////////////////////
        
        // Uncomment to display various parameters
        if (!rank) {
        //     std::cout << YLW << "maxDepth: " << m_uiMaxDepth << NRM << std::endl;
        //     std::cout << YLW << "wavelet_tol: " << wavelet_tol << NRM << std::endl;
        //     std::cout << YLW << "partition_tol: " << partition_tol << NRM << std::endl;
        //     std::cout << YLW << "eleOrder: " << eOrder << NRM << std::endl;
        }

        _InitializeHcurve(m_uiDim);
        RefElement refEl(m_uiDim,eOrder);
        
        // This is the tricky part. Octree generation could be based on a function or other variable. This will need to be generated.
        // But for now just use this
        ot::DA* octDA=new ot::DA("""* cpp_gen_string(coefficients[1].value[1]) *""",1,comm,eOrder,wavelet_tol,100,partition_tol,ot::FEM_CG);
        
        // Variable info will also be generated, but for now assume a single scalar variable
        std::vector<double> uSolVec;
        octDA->createVector(uSolVec,false,false,DOF);
        double *uSolVecPtr=&(*(uSolVec.begin()));

        FemshopDendroSkeleton::LHSMat lhsMat(octDA,1);
        lhsMat.setProblemDimensions(domain_min,domain_max);
        lhsMat.setGlobalDofVec(uSolVecPtr);
        

        FemshopDendroSkeleton::RHSVec rhsVec(octDA,1);
        rhsVec.setProblemDimensions(domain_min,domain_max);
        rhsVec.setGlobalDofVec(uSolVecPtr);
        
        // This assumes some things
        lhsMat.setBdryFunction("""* cpp_gen_string(prob.bc_func[1,1]) *""");
        rhsVec.setBdryFunction("""* cpp_gen_string(prob.bc_func[1,1]) *""");
        
        // Allocate dofs
        """*alloc_dof*"""
        double * rhs=octDA->getVecPointerToDof(uSolVecPtr,VAR::M_UI_RHS, false,false); // linear part
        
        // Init dofs
        """*set_dof*"""
        octDA->setVectorByFunction(rhs,zero_init,false,false,1); // zeros
        
        """*solvepart*"""
        
        // Output
        //////////////will be generated/////////////////////////////////////////////
        #include "Output.cpp"
        ////////////////////////////////////////////////////////////////////////////
        
        octDA->destroyVector(uSolVec);

        if(!rank)
            std::cout<<" End of computation. "<<std::endl;

        delete octDA;

        MPI_Finalize();
        return 0;
    }
    """
    print(file, content);
end

function dendro_linear_skeleton(file, headerfile)
    varpart = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar.vals)
            varpart = varpart*"M_UI"*string(variables[i].symvar.vals[j]);
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
        end
    end
    coefpart = "";
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar.vals)
            coefpart = coefpart*"M_UI"*string(coefficients[i].symvar.vals[j])*", ";
        end
    end
    
    content = """
    #include "linear_skel.h"

    FemshopDendroSkeleton::RHSVec::RHSVec(ot::DA* da,unsigned int dof) : feVector(da,dof)
    {
        const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
        imV1=new double[nPe];
        imV2=new double[nPe];
    }

    FemshopDendroSkeleton::RHSVec::~RHSVec()
    {
        delete [] imV1;
        delete [] imV2;

        imV1=NULL;
        imV2=NULL;
    }

    void FemshopDendroSkeleton::RHSVec::elementalComputVec(const VECType* in,VECType* out, double*coords,double scale)
    {
        const RefElement* refEl=m_uiOctDA->getReferenceElement();
        const double * Q1d=refEl->getQ1d();
        const double * QT1d=refEl->getQT1d();
        const double * Dg=refEl->getDg1d();
        const double * W1d=refEl->getWgq();

        const unsigned int eleOrder=refEl->getOrder();
        const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
        const unsigned int nrp=eleOrder+1;

        Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
        Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

        const double refElSz=refEl->getElementSz();
        
        const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
        const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
        const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());
        
        const double Jx = 1.0/(refElSz/(double (szX)));
        const double Jy = 1.0/(refElSz/(double (szY)));
        const double Jz = 1.0/(refElSz/(double (szZ)));
        
        double* Qx = {0};
        double* Qy = {0};
        double* Qz = {0};
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Linear.cpp"
        ////////////////////////////////////////////////////////////////////////////
    }
    
    void FemshopDendroSkeleton::RHSVec::setBdryFunction(std::function<void(double,double,double,double*)> bdry){
        bdry_function = bdry;
    }
    
    void FemshopDendroSkeleton::RHSVec::setGlobalDofVec(double* gdv){
        grandDofVecPtr = gdv;
    }

    bool FemshopDendroSkeleton::RHSVec::preComputeVec(const VECType* in,VECType* out, double scale)
    {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);
        
        double x,y,z,val;
        for(unsigned int i=0;i<bdyIndex.size();i++){
            x = bdyCoords[i*3+0];
            y = bdyCoords[i*3+1];
            z = bdyCoords[i*3+2];
            bdry_function(x,y,z,&val);
            
            out[bdyIndex[i]] = val;
        }

        return true;
    }

    bool FemshopDendroSkeleton::RHSVec::postComputeVec(const VECType* in,VECType* out, double scale) {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        double x,y,z,val;
        for(unsigned int i=0;i<bdyIndex.size();i++){
            x = bdyCoords[i*3+0];
            y = bdyCoords[i*3+1];
            z = bdyCoords[i*3+2];
            bdry_function(x,y,z,&val);
            
            out[bdyIndex[i]] = val;
        }

        return true;
    }


    double FemshopDendroSkeleton::RHSVec::gridX_to_X(double x)
    {
        double Rg_x=((1u<<m_uiMaxDepth)-0);
        return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
    }

    double FemshopDendroSkeleton::RHSVec::gridY_to_Y(double y)
    {
        double Rg_y=((1u<<m_uiMaxDepth)-0);
        return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
    }


    double FemshopDendroSkeleton::RHSVec::gridZ_to_Z(double z)
    {
        double Rg_z=((1u<<m_uiMaxDepth)-0);
        return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
    }
    """
    print(file, content);
    content = """
    #ifndef DENDRO_5_0_LINEAR_SKEL_H
    #define DENDRO_5_0_LINEAR_SKEL_H

    #include "oda.h"
    #include "feVector.h"

    namespace FemshopDendroSkeleton
    {
        class RHSVec : public feVector<RHSVec>{

        private:

            double * imV1;
            double * imV2;
            // function for boundary
            std::function<void(double,double,double,double*)> bdry_function;
            // coefficient vectors
            double* grandDofVecPtr;
            """*"enum VAR{"*varpart*coefpart*"M_UI_RHS};"*"""

        public:
            RHSVec(ot::DA* da,unsigned int dof=1);
            ~RHSVec();
            
            /**@biref elemental compute vec for rhs*/
            virtual void elementalComputVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);
            
            /**@brief set boundary function*/	
            void setBdryFunction(std::function<void(double,double,double,double*)> bdry);
            
            /**@brief set pointer to global dof vector*/	
                void setGlobalDofVec(double* gdv);
            
            bool preComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            bool postComputeVec(const VECType* in,VECType* out, double scale=1.0);
            
            /**@brief octree grid x to domin x*/
            double gridX_to_X(double x);
            /**@brief octree grid y to domin y*/
            double gridY_to_Y(double y);
            /**@brief octree grid z to domin z*/
            double gridZ_to_Z(double z);
        };
    }

    #endif //DENDRO_5_0_LINEAR_SKEL_H
    """
    print(headerfile, content);
end

function dendro_bilinear_skeleton(file, headerfile)
    varpart = "";
    for i=1:length(variables)
        for j=1:length(variables[i].symvar.vals)
            varpart = varpart*"M_UI"*string(variables[i].symvar.vals[j]);
            if i==1 && j==1
                varpart = varpart*"=0";
            end
            varpart = varpart*", ";
        end
    end
    coefpart = "";
    for i=1:length(coefficients)
        for j=1:length(coefficients[i].symvar.vals)
            coefpart = coefpart*"M_UI"*string(coefficients[i].symvar.vals[j])*", ";
        end
    end
    
    content = """
    #include "bilinear_skel.h"

    FemshopDendroSkeleton::LHSMat::LHSMat(ot::DA* da,unsigned int dof) : feMatrix(da,dof)
    {
        const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
        imV1=new double[nPe];
        imV2=new double[nPe];

        Qx=new double[nPe];
        Qy=new double[nPe];
        Qz=new double[nPe];
    }

    FemshopDendroSkeleton::LHSMat::~LHSMat()
    {
        delete [] imV1;
        delete [] imV2;

        delete [] Qx;
        delete [] Qy;
        delete [] Qz;

        imV1=NULL;
        imV2=NULL;

        Qx=NULL;
        Qy=NULL;
        Qz=NULL;
    }

    void FemshopDendroSkeleton::LHSMat::elementalMatVec(const VECType* in,VECType* out, double*coords,double scale)
    {
        const RefElement* refEl=m_uiOctDA->getReferenceElement();

        const double * Q1d=refEl->getQ1d();
        const double * QT1d=refEl->getQT1d();
        const double * Dg=refEl->getDg1d();
        const double * DgT=refEl->getDgT1d();
        const double * W1d=refEl->getWgq();

        const unsigned int eleOrder=refEl->getOrder();
        const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
        const unsigned int nrp=eleOrder+1;

        Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
        Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

        const double refElSz=refEl->getElementSz();
        
        const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
        const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
        const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


        const double Jx = 1.0/(refElSz/(double (szX)));
        const double Jy = 1.0/(refElSz/(double (szY)));
        const double Jz = 1.0/(refElSz/(double (szZ)));
        
        //////////////will be generated/////////////////////////////////////////////
        #include "Bilinear.cpp"
        ////////////////////////////////////////////////////////////////////////////
    }
    
    void FemshopDendroSkeleton::LHSMat::setBdryFunction(std::function<void(double,double,double,double*)> bdry){
        bdry_function = bdry;
    }
    
    void FemshopDendroSkeleton::LHSMat::setGlobalDofVec(double* gdv){
        grandDofVecPtr = gdv;
    }

    bool FemshopDendroSkeleton::LHSMat::preMatVec(const VECType* in,VECType* out,double scale)
    {
        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        for(unsigned int i=0;i<bdyIndex.size();i++){
            
            out[bdyIndex[i]] = in[bdyIndex[i]]; // Dirichlet BC
        }

        return true;
    }

    bool FemshopDendroSkeleton::LHSMat::postMatVec(const VECType* in,VECType* out,double scale) {

        // apply boundary conditions.
        std::vector<unsigned int> bdyIndex;
        std::vector<double> bdyCoords;

        m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

        for(unsigned int i=0;i<bdyIndex.size();i++){
            
            out[bdyIndex[i]] = in[bdyIndex[i]]; // Dirichlet BC
        }

        return true;
    }


    double FemshopDendroSkeleton::LHSMat::gridX_to_X(double x)
    {
        double Rg_x=((1u<<m_uiMaxDepth)-0);
        return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
    }

    double FemshopDendroSkeleton::LHSMat::gridY_to_Y(double y)
    {
        double Rg_y=((1u<<m_uiMaxDepth)-0);
        return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
    }


    double FemshopDendroSkeleton::LHSMat::gridZ_to_Z(double z)
    {
        double Rg_z=((1u<<m_uiMaxDepth)-0);
        return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
    }

    int FemshopDendroSkeleton::LHSMat::cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var)
    {
        double resid,alpha,beta,rho,rho_1;
        int status=1; // 0 indicates it has solved the system within the specified max_iter, 1 otherwise.

        const unsigned int local_dof=m_uiOctDA->getLocalNodalSz();

        MPI_Comm globalComm=m_uiOctDA->getGlobalComm();

        if(m_uiOctDA->isActive())
        {

            int activeRank=m_uiOctDA->getRankActive();
            int activeNpes=m_uiOctDA->getNpesActive();

            MPI_Comm activeComm=m_uiOctDA->getCommActive();

            double* p;
            double* z;
            double* q;
            double* Ax;
            double* Ap;
            double* r0;
            double* r1;

            m_uiOctDA->createVector(p);
            m_uiOctDA->createVector(z);
            m_uiOctDA->createVector(q);

            m_uiOctDA->createVector(Ax);
            m_uiOctDA->createVector(Ap);
            m_uiOctDA->createVector(r0);
            m_uiOctDA->createVector(r1);

            double normb = normLInfty(b,local_dof,activeComm);
            par::Mpi_Bcast(&normb,1,0,activeComm);

            if(!activeRank)
                std::cout<<"normb = "<<normb<<std::endl;

            matVec(x,Ax);

            /*char fPrefix[256];
            sprintf(fPrefix,"%s_%d","cg",0);
            const char * varNames[]={"U"};
            const double * var[]={Ax};
            io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);
            */
            for(unsigned int i=0;i<local_dof;i++)
            {
                r0[i]=b[i]-Ax[i];
                p[i]=r0[i];
            }

            if (normb == 0.0)
                normb = 1;

            double normr=normLInfty(r0,local_dof,activeComm);
            par::Mpi_Bcast(&normr,1,0,activeComm);
            if(!activeRank) std::cout<<"initial residual : "<<(normr/normb)<<std::endl;

            if ((resid = normr / normb) <= tol) {
                tol = resid;
                max_iter = 0;

                m_uiOctDA->destroyVector(p);
                m_uiOctDA->destroyVector(z);
                m_uiOctDA->destroyVector(q);

                m_uiOctDA->destroyVector(Ax);
                m_uiOctDA->destroyVector(Ap);
                m_uiOctDA->destroyVector(r0);
                m_uiOctDA->destroyVector(r1);

                status=0;
            }

            if(status!=0)
            {
                for(unsigned int i=1;i<=max_iter;i++)
                {
                    matVec(p,Ap);

                    alpha=(dot(r0,r0,local_dof,activeComm)/dot(p,Ap,local_dof,activeComm));
                    par::Mpi_Bcast(&alpha,1,0,activeComm);

                    //if(!activeRank) std::cout<<"rank: " <<activeRank<<" alpha: "<<alpha<<std::endl;
                    for(unsigned int e=0;e<local_dof;e++)
                    {
                        x[e]+=alpha*p[e];
                        r1[e]=r0[e]-alpha*Ap[e];
                    }

                    normr=normLInfty(r1,local_dof,activeComm);
                    par::Mpi_Bcast(&normr,1,0,activeComm);

                    if((!activeRank) && (i%10==0)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;

                    if ((resid = normr / normb) <= tol) {

                        if((!activeRank)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;
                        tol = resid;
                        m_uiOctDA->destroyVector(p);
                        m_uiOctDA->destroyVector(z);
                        m_uiOctDA->destroyVector(q);

                        m_uiOctDA->destroyVector(Ax);
                        m_uiOctDA->destroyVector(Ap);
                        m_uiOctDA->destroyVector(r0);
                        m_uiOctDA->destroyVector(r1);

                        status=0;
                        break;
                    }

                    beta=(dot(r1,r1,local_dof,activeComm)/dot(r0,r0,local_dof,activeComm));
                    par::Mpi_Bcast(&beta,1,0,activeComm);

                    //if(!activeRank) std::cout<<"<r_1,r_1> : "<<dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)<<" <r_0,r_0>: "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" beta "<<beta<<std::endl;

                    for(unsigned int e=0;e<local_dof;e++)
                    {
                        p[e]=r1[e]+beta*p[e];
                        r0[e]=r1[e];
                    }
                }

                if(status!=0)
                {
                    tol = resid;
                    m_uiOctDA->destroyVector(p);
                    m_uiOctDA->destroyVector(z);
                    m_uiOctDA->destroyVector(q);

                    m_uiOctDA->destroyVector(Ax);
                    m_uiOctDA->destroyVector(Ap);
                    m_uiOctDA->destroyVector(r0);
                    m_uiOctDA->destroyVector(r1);
                    status=1;
                }
            }
        }

        // bcast act as a barrier for active and inactive meshes.
        par::Mpi_Bcast(&tol,1,0,globalComm);
        return status;
    }
    """
    print(file, content);
    content = """
    #ifndef DENDRO_5_0_BILINEAR_SKEL_H
    #define DENDRO_5_0_BILINEAR_SKEL_H

    #include "oda.h"
    #include "feMatrix.h"

    namespace FemshopDendroSkeleton
    {
        class LHSMat : public feMatrix<LHSMat>{

        private:
            // some additional work space variables to perform elemental MatVec
            double* imV1;
            double* imV2;
            double* Qx;
            double* Qy;
            double* Qz;
            // function for boundary
            std::function<void(double,double,double,double*)> bdry_function;
            // coefficient vectors
            double* grandDofVecPtr;
            """*"enum VAR{"*varpart*coefpart*"M_UI_RHS};"*"""

        public:
            /**@brief: constructor*/
            LHSMat(ot::DA* da,unsigned int dof=1);

            /**@brief default destructor*/
            ~LHSMat();

            /**@biref elemental matvec*/
            virtual void elementalMatVec(const VECType* in,VECType* out, double*coords=NULL,double scale=1.0);
            
            /**@brief set boundary function*/	
            void setBdryFunction(std::function<void(double,double,double,double*)> bdry);
            
            /**@brief set pointer to global dof vector*/	
                void setGlobalDofVec(double* gdv);
            
            /**@brief things need to be performed before matvec (i.e. coords transform)*/
            bool preMatVec(const VECType* in,VECType* out,double scale=1.0);

            /**@brief things need to be performed after matvec (i.e. coords transform)*/
            bool postMatVec(const VECType* in,VECType* out,double scale=1.0);

            /**@brief octree grid x to domin x*/
            double gridX_to_X(double x);
            /**@brief octree grid y to domin y*/
            double gridY_to_Y(double y);
            /**@brief octree grid z to domin z*/
            double gridZ_to_Z(double z);

            int cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var=0);
        };
    }

    #endif //DENDRO_5_0_BILINEAR_SKEL_H
    """
    print(headerfile, content);
end

function dendro_cmake_file(file)
    content = """
cmake_minimum_required(VERSION 2.8)
project("""*genFileName*""")

set("""*genFileName*"""_INC include/linear_skel.h
            include/bilinear_skel.h
        )

set("""*genFileName*"""_SRC src/linear_skel.cpp
            src/bilinear_skel.cpp
        )

set(SOURCE_FILES src/"""*genFileName*""".cpp \${"""*genFileName*"""_INC} \${"""*genFileName*"""_SRC})
add_executable("""*genFileName*""" \${SOURCE_FILES})
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/include/test)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories("""*genFileName*""" PRIVATE \${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories("""*genFileName*""" PRIVATE \${MPI_INCLUDE_PATH})
target_include_directories("""*genFileName*""" PRIVATE \${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories("""*genFileName*""" PRIVATE \${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries("""*genFileName*""" dendro5 \${LAPACK_LIBRARIES} \${MPI_LIBRARIES} m)

""";
    print(file, content);
end

function dendro_readme_file(file)
    content = """
Basic instructions for compiling this generated code with dendro.
    1. Place this generated directory in the Dendro directory.
    
    2. Append the following line to the Dendro CMakeLists.txt file:
        add_subdirectory("""*uppercasefirst(genFileName)*""")
    
    3. Remake Dendro as usual. For example, to build to the directory "build",
       enter the Dendro directory and do the following.
        a. If it doesn't exist, make the build directory: "mkdir build"
        b. "cd build"
        c. "ccmake .."
        d. "make"
        e. "cd """*uppercasefirst(genFileName)*""" "
    
    4. Run the executable with "./"""*genFileName*""" "
    
    *Note: The directory will have an upper case first letter. 
           The executable will match the name passed to @language().
"""
    print(file, content);
end