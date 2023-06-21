#pragma once

#include "Waveguide.h"
#include "core/physicalModel/Bound.h"

namespace semba {
namespace source {
namespace port {

class WaveguideRectangular: public Waveguide {
public:
    WaveguideRectangular(
            const std::unique_ptr<Magnitude::Magnitude>& magnitude,
            const Target& elem,
            const ExcitationMode excMode,
            const std::pair<size_t,size_t> mode);
    WaveguideRectangular(const WaveguideRectangular&) = default;
    virtual ~WaveguideRectangular() = default;

    virtual std::unique_ptr<Source> clone() const override {
        return std::make_unique<WaveguideRectangular>(*this);
    }

    //void set(const Target&);

    std::string getName() const override;
    //math::Real getWidth() const;
    //math::Real getHeight() const;

    //math::CVecR3 getOrigin() const override;
    //math::CVecR3 getWeight(const math::CVecR3& pos) const override;
};

}
}
} 

