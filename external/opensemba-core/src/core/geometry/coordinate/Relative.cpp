#include "Relative.h"

namespace semba {
namespace geometry {
namespace coordinate {

Relative::Relative(const Id id,
                   const math::CVecR3& rel)
:   Identifiable<Id>(id) {
    for (std::size_t d = 0; d < 3; d++) {
        this->pos()(d) = (int)std::floor(rel(d));
        this->rel_ (d) = rel(d) - this->pos()(d);
    }
}

Relative::Relative(const Id id,
                   const math::CVecI3& pos,
                   const math::CVecR3& rel)
:   Identifiable<Id>(id),
    math::CVecI3(pos) {
    rel_ = rel;
}

Relative::Relative(const math::CVecR3& rel) {
    for (std::size_t d = 0; d < 3; d++) {
        this->pos()(d) = (int)std::floor(rel(d));
        this->rel_(d) = rel(d) - this->pos()(d);
    }
}

Relative::Relative(const Relative& rhs)
:   Identifiable<Id>(rhs),
    math::CVecI3(rhs) {
    rel_ = rhs.rel_;
}

Relative& Relative::operator=(const Relative& rhs) {
    if (this == &rhs)
        return *this;

    CoordI3::operator=(rhs);
    rel_ = rhs.rel_;
    return *this;
}

bool Relative::operator==(const Base& rhs) const {
    if (!Coordinate<math::Int,3>::operator==(rhs)) {
        return false;
    }
    const Relative* rhsPtr = rhs.castTo<Relative>();
    bool res = true;
    res &= (this->rel_ == rhsPtr->rel_);
    return res;
}

CoordR3* Relative::toUnstructured(const Grid3& grid) const {
    math::CVecR3 pos = grid.getPos(*this);
    for (std::size_t d = 0; d < 3; d++) {
        math::Real length = rel_(d);
        math::Int cellDir;
        if (this->pos()(d) == grid.getNumCells()(d)) {
            cellDir = this->pos()(d);
        }
        else {
            cellDir = this->pos()(d) + 1;
        }
        math::Real posAux = grid.getPos(d, cellDir);
        math::Real step = posAux - pos(d);
        pos(d) += step * length;
    }
    return new CoordR3(this->getId(), pos);
}

} 
} 
} 
