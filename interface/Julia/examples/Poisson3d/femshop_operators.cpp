
#include "femshop_operators.h"

VECType* femshopDendroOp_negative(VECType* out, const RefElement* refEl){
    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    for(unsigned int i=0;i<nPe;i++){
        out[i] = -out[i];
    }
    return out;
}

VECType* femshopDendroOp_mass_operator(const VECType* in,VECType* out, const RefElement* refEl, const double Jx, 
                                                const double Jy, const double Jz, double* imV1, double* imV2, 
                                                double* Qx, double* Qy, double* Qz){
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * W1d=refEl->getWgq();
    
    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    const unsigned int nrp=eleOrder+1;
    
    // interpolate to quadrature points.
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,out);
    
    for(unsigned int k=0;k<(eleOrder+1);k++){
        for(unsigned int j=0;j<(eleOrder+1);j++){
            for(unsigned int i=0;i<(eleOrder+1);i++){
                out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=(Jx*Jy*Jz*W1d[i]*W1d[j]*W1d[k]);
            }
        }
    }
    
    // apply transpose operator
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,out);
    
    return out;
}

VECType* femshopDendroOp_stiffness_operator(const VECType* in,VECType* out, const RefElement* refEl, const double Jx, 
                                                const double Jy, const double Jz, double* imV1, double* imV2, 
                                                double* Qx, double* Qy, double* Qz){
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * DgT=refEl->getDgT1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    const unsigned int nrp=eleOrder+1;
    
    //x derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qx);

    //y derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qy);

    //z derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,Qz);

    for(unsigned int k=0;k<(eleOrder+1);k++){
        for(unsigned int j=0;j<(eleOrder+1);j++){
            for(unsigned int i=0;i<(eleOrder+1);i++){
                Qx[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jy*Jz)/Jx)*W1d[i]*W1d[j]*W1d[k]);
                Qy[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jz)/Jy)*W1d[i]*W1d[j]*W1d[k]);
                Qz[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jy)/Jz)*W1d[i]*W1d[j]*W1d[k]);
            }
        }
    }
    
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,Qx,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qx);
    
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qy,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qy);
    
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qz,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,Qz);
    
    for(unsigned int i=0;i<nPe;i++){
        out[i]=Qx[i]+Qy[i]+Qz[i];
    }
    
    return out;
}