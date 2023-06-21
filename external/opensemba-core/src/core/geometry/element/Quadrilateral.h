#pragma once

#include "Surface.h"

namespace semba {
namespace geometry {
namespace element {

class QuadrilateralBase : public virtual SurfaceBase {
public:
    virtual ~QuadrilateralBase() = default;
    std::size_t numberOfFaces   () const { return 4; }
    std::size_t numberOfVertices() const { return 4; }

    std::size_t numberOfSideVertices(const std::size_t f = 0) const { 
        return 2; 
    }
};

template<class T>
class Quadrilateral : public virtual Surface<T>,
                      public virtual QuadrilateralBase {
public:
    virtual ~Quadrilateral() = default;
};

} 

typedef element::QuadrilateralBase         Qua;
typedef element::Quadrilateral<math::Real> QuaR;
typedef element::Quadrilateral<math::Int > QuaI;

} 
} 

