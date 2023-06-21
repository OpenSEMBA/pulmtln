#pragma once

#include "TEM.h"

namespace semba {
namespace source {
namespace port {

class TEMCoaxial : public TEM {
public:
    TEMCoaxial(
            const std::unique_ptr<Magnitude::Magnitude>& magnitude,
            const Target& elem,
            const ExcitationMode excMode,
            const math::CVecR3& origin,
            const math::Real innerRadius,
            const math::Real outerRadius);
    TEMCoaxial(const TEMCoaxial& rhs);
    virtual ~TEMCoaxial() = default;

    virtual std::unique_ptr<Source> clone() const override {
        return std::make_unique<TEMCoaxial>(*this);
    }

    //void set(const Target&);

    std::string getName() const override { return "Coaxial_TEM_port"; }

    //math::CVecR3 getOrigin() const override;
    //math::CVecR3 getWeight(const math::CVecR3& pos) const override;

private:
    math::CVecR3 origin_;
    math::Real innerRadius_, outerRadius_;
};

}
}
} 

