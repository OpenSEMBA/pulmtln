

#pragma once

#include "Triangle.h"

namespace semba {
namespace geometry {
namespace element {

class Triangle3 : public Triangle {
public:
    Triangle3(const Id id,
              const CoordR3* v[3],
              const Layer* lay = nullptr,
              const Model* mat = nullptr);
    Triangle3(const Triangle3& rhs);
    virtual ~Triangle3() = default;

    static const std::size_t sizeOfCoordinates = 3;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Triangle3>(*this);
    }

    std::size_t numberOfCoordinates() const override { return sizeOfCoordinates; }

    std::size_t numberOfSideCoordinates(const std::size_t f = 0) const override {
        return 2;
    }

    const CoordR3* getV     (const std::size_t i) const override { return v_[i]; }
    const CoordR3* getVertex(const std::size_t i) const override;

    const CoordR3* getSideV     (const std::size_t face,
                                 const std::size_t i) const override;
    const CoordR3* getSideVertex(const std::size_t face,
                                 const std::size_t i) const override;

    math::Real getArea() const;

    void setV(const std::size_t i, const CoordR3*) override;

    void check() const;

    virtual std::unique_ptr<Element<math::Int >> toStructured(
        const CoordI3Group&,
        const Grid3&,
        const math::Real = Grid3::tolerance) const override;

    virtual std::unique_ptr<Element<math::Real>> toUnstructured(
        const CoordR3Group&,
        const Grid3&) const override;

protected:
    static const math::simplex::Triangle<1> geo;

    std::array<const CoordR3*,3> v_;
};

} 

typedef element::Triangle3 Tri3;

} 
} 

