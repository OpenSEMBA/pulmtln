#pragma once

#include "Coordinate.h"

namespace semba {
namespace geometry {
namespace coordinate {

class Relative : public virtual Coordinate<math::Int,3> {
public:
    Relative() = default;
    Relative(const Id, const math::CVecR3&);
    Relative(const Id, const math::CVecI3&, const math::CVecR3&);
    Relative(const math::CVecR3&);
    Relative(const Relative&);
    virtual ~Relative() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Relative>(*this);
    }

    Relative& operator=(const Relative& rhs);

    bool operator==(const Base& rhs) const override;

    math::CVecR3&       rel()       { return rel_; }
    const math::CVecR3& rel() const { return rel_; }

    CoordR3* toUnstructured(const Grid3&) const override;

private:
    math::CVecR3 rel_;
};

} 

typedef coordinate::Relative CoordRel;

} 
} 

