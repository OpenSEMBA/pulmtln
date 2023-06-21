#pragma once

#include "OutputRequest.h"

namespace semba {
namespace outputRequest {

class FarField : public virtual OutputRequest {
public:
    FarField(const Domain& domain,
             const std::string& name,
             const Target& box,
             const math::Real iTh, const math::Real fTh, const math::Real sTh,
             const math::Real iPhi, const math::Real fPhi,
             const math::Real sPhi);
    virtual ~FarField() = default;

    std::unique_ptr<OutputRequest> clone() const override {
        return std::make_unique<FarField>(*this);
    }

    math::Real initialTheta, finalTheta, stepTheta;
    math::Real initialPhi, finalPhi, stepPhi;
};

} 
} 

