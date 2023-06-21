#pragma once

#include "Source.h"

namespace semba {
namespace source {

class Dipole : public Source {
public:
    Dipole(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
           const Target& elem,
           math::Real   length,
           math::CVecR3 orientation,
           math::CVecR3 position);
    
    virtual std::unique_ptr<Source> clone() const override {
        return std::make_unique<Dipole>(*this);
    }

    std::string getName() const override { return "Dipole"; }
private:
    math::Real length_;
    math::CVecR3 orientation_;
    math::CVecR3 position_;
};

}
} 

