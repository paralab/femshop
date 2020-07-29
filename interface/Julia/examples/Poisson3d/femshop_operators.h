#ifndef DENDRO_5_0_FEMSHOP_OPERATORS_H
#define DENDRO_5_0_FEMSHOP_OPERATORS_H

#include "oda.h"
#include "refel.h"
#include "tensor.h"

VECType* negative(VECType* out, const RefElement* refEl);

VECType* mass_operator(const VECType* in,VECType* out, const RefElement* refEl, const double Jx, 
                        const double Jy, const double Jz, double* imV1, double* imV2, 
                        double* Qx, double* Qy, double* Qz);
                                            
VECType* stiffness_operator(const VECType* in,VECType* out, const RefElement* refEl, const double Jx, 
                            const double Jy, const double Jz, double* imV1, double* imV2, 
                            double* Qx, double* Qy, double* Qz);

#endif //DENDRO_5_0_FEMSHOP_OPERATORS_H