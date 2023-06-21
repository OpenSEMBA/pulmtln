#include "Geometric.h"

namespace semba {
namespace geometry {
namespace mesh {

Geometric::Geometric(const Grid3& grid)
:   grid_(grid) {

}

Geometric::Geometric(const Grid3& grid,
                     const CoordR3Group& cG,
                     const ElemRGroup& elem,
                     const LayerGroup& layers)
:   Unstructured(cG, elem, layers),
    grid_(grid) {

}

Geometric::Geometric(const Geometric& rhs)
:   Unstructured(rhs),
    grid_(rhs.grid_) {

}

Geometric::Geometric(const Grid3& grid, const Unstructured& unstructured) :
	Unstructured(unstructured),
	grid_(grid) {}

Geometric& Geometric::operator=(const Geometric& rhs) {
    if(this == &rhs) {
        return *this;
    }
    Unstructured::operator=(rhs);
    grid_ = rhs.grid_;

    return *this;
}

Structured* Geometric::getMeshStructured(const math::Real tol) const {
    return Unstructured::getMeshStructured(grid_, tol);
}

void Geometric::applyScalingFactor(const math::Real factor) {
    Unstructured::applyScalingFactor(factor);
    grid_.applyScalingFactor(factor);
}

} 
} 
} 
