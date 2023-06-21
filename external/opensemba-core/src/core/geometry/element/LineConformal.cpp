#include "LineConformal.h"

#include <iterator>
#include <algorithm>

namespace semba {
namespace geometry {
namespace element {

LineConformal::LineConformal(const Id id,
                             std::array<const coordinate::Coordinate<math::Int, 3>*, 2> v,
                             const math::CVecR3& norm,
                             const Layer* lay,
                             const Model* mat)
:   Identifiable<Id>(id),
    Elem(lay, mat),
    LinI2(v) 
{
    checkCoordinates();
    norm_  = norm;
}

LineConformal::LineConformal(
    const Id id,
    const coordinate::Coordinate<math::Int, 3>* v[2],
    const math::CVecR3& norm,
    const Layer* lay,
    const Model* mat)
{
    std::array<const coordinate::Coordinate<math::Int, 3>*, 2> vArr;
    std::copy(v, v + 2, vArr.begin());
    *this = LineConformal{id, vArr, norm, lay, mat};
}

LineConformal::LineConformal(
    std::array<const coordinate::Coordinate<math::Int, 3>*, 2> v,
    const math::CVecR3& norm,
    const Layer* lay,
    const Model* mat):   
    LinI2(ElemId(0), v, lay, mat) 
{
    checkCoordinates();
    norm_  = norm;
}

LineConformal::LineConformal(const LineConformal& rhs)
:   Identifiable<Id>(rhs),
    Elem(rhs),
    LinI2(rhs) 
{
    norm_  = rhs.norm_;
}

const CoordConf* LineConformal::getV(const std::size_t i) const {
	const coordinate::Coordinate<math::Int,3>* coord;
	coord = Line2<math::Int>::getV(i);
    return coord->castTo<CoordConf>();
}

void LineConformal::setV(const std::size_t i, const CoordI3* coord) {
    LinI2::setV(i, coord);
    checkCoordinates();
}

void LineConformal::checkCoordinates() {
    for(std::size_t i = 0; i < this->numberOfCoordinates(); i++) {
        if (!this->getV(i)->is<CoordConf>()) {
            throw std::logic_error("Coord not conformal");
        }
    }
}

std::unique_ptr<ElemI> LineConformal::toStructured(
    const CoordI3Group& cG,
    const Grid3& grid, const math::Real tol) const 
{
    return std::make_unique<LineConformal>(
        this->getId(),
        this->vertexToStructured(cG, grid, tol).data(),
        norm_,
        this->getLayer(),
        this->getModel());
}

std::unique_ptr<ElemR> LineConformal::toUnstructured(
    const CoordR3Group& cG,
    const Grid3& grid) const {
    throw std::logic_error("LineConformal::toUnstructured operation not permitted");
}

} 
} 
} 
