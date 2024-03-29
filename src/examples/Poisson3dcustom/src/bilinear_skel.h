//This file was generated by Femshop.

/*

*/
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
        enum VAR{M_UI_u_1=0, M_UI_f_1, M_UI_RHS};
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
