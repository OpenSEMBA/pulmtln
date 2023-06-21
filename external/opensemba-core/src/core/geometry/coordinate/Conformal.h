#pragma once

#include "Coordinate.h"

namespace semba {
namespace geometry {
namespace coordinate {

class Conformal : public virtual Coordinate<math::Int,3> {
public:
    Conformal();
    Conformal(const Id id_,
              const math::CVecI3& pos,
              const math::Constants::CartesianAxis dir,
              const math::Real length);
    Conformal(const math::Constants::CartesianAxis dir,
              const math::Real length);
    Conformal(const Conformal& rhs);
    virtual ~Conformal();

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Conformal>(*this);
    }

    Conformal& operator=(const Conformal& rhs);

    bool operator==(const Base& rhs) const override;

    math::Constants::CartesianAxis getDir   () const { return dir_;    }
    math::Real                     getLength() const { return length_; }

    CoordR3* toUnstructured(const Grid3&) const override;

private:
    math::Constants::CartesianAxis dir_;
    math::Real                     length_;
};

} 

typedef coordinate::Conformal CoordConf;

} 
} 

