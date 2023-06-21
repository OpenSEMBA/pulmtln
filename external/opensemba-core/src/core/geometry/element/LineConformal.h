#pragma once

#include "core/geometry/coordinate/Conformal.h"

#include "Line2.h"

namespace semba {
namespace geometry {
namespace element {

class LineConformal : public virtual Line2<math::Int> {
public:
    LineConformal(const Id id,
                  std::array<const coordinate::Coordinate<math::Int,3>*,2>,
                  const math::CVecR3& norm,
                  const Layer* lay = nullptr,
                  const Model* mat = nullptr);
    LineConformal(const Id id,
                  const coordinate::Coordinate<math::Int, 3>* v[2],
                  const math::CVecR3& norm,
                  const Layer* lay = nullptr,
                  const Model* mat = nullptr);
    LineConformal(std::array<const coordinate::Coordinate<math::Int, 3>*, 2>,
                  const math::CVecR3& norm,
                  const Layer* lay = nullptr,
                  const Model* mat = nullptr);
    LineConformal(const LineConformal& rhs);
    virtual ~LineConformal() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<LineConformal>(*this);
    }

    math::CVecR3 getNorm () const { return norm_;  }

    const CoordConf* getV(const std::size_t i) const override;

    void setV(const std::size_t i, const CoordI3* coord) override;

    std::unique_ptr<ElemI> toStructured(const CoordI3Group&,
        const Grid3&,
        const math::Real = Grid3::tolerance) const override;
    std::unique_ptr<ElemR> toUnstructured(const CoordR3Group&,
        const Grid3&) const override;

private:
    void checkCoordinates();
    math::CVecR3 norm_;
};

} 

typedef element::LineConformal LinConf;

} 
} 

