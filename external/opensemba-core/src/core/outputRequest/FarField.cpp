#include "FarField.h"

namespace semba {
namespace outputRequest {

FarField::FarField(const Domain& domain,
                   const std::string& name,
                   const Target& elem,
                   const math::Real iTh, const math::Real fTh, const math::Real sTh,
                   const math::Real iPhi, const math::Real fPhi, const math::Real sPhi)
:   OutputRequest(Type::electricFarField, name, domain, elem) 
{
    initialTheta = iTh;
    finalTheta = fTh;
    stepTheta = sTh;
    initialPhi = iPhi;
    finalPhi = fPhi;
    stepPhi = sPhi;
}

} 
} 
