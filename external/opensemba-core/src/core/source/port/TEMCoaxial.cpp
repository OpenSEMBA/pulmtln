#include "TEMCoaxial.h"
#include "core/geometry/Bound.h"

namespace semba {
namespace source {
namespace port {

TEMCoaxial::TEMCoaxial(
        const std::unique_ptr<Magnitude::Magnitude>& magnitude,
        const Target& elem,
        const ExcitationMode excMode,
        const math::CVecR3& origin,
        const math::Real innerRadius,
        const math::Real outerRadius) :
        TEM(magnitude, elem, excMode) 
{
    origin_ = origin;
    innerRadius_ = innerRadius;
    outerRadius_ = outerRadius;
}

TEMCoaxial::TEMCoaxial(const TEMCoaxial& rhs) : 
    TEM(rhs) 
{
    origin_ = rhs.origin_;
    innerRadius_ = rhs.innerRadius_;
    outerRadius_ = rhs.outerRadius_;
}
//
//math::CVecR3 TEMCoaxial::getOrigin() const {
//    return origin_;
//}
//
//math::CVecR3 TEMCoaxial::getWeight(
//        const math::CVecR3& pos) const {
//    // Return normalized weights for electric field components.
//    const math::Real rho = (pos - getOrigin()).norm();
//    switch (getExcitationMode()) {
//    case ExcitationMode::voltage:
//    {
//        const math::CVecR3 rhoHat = (pos - getOrigin()).normalize();
//        return rhoHat / (rho * log(outerRadius_/innerRadius_));
//    }
//    case ExcitationMode::current:
//    {
//        const math::CVecR3 phiHat = (math::CVecR3(0,0,1) ^ pos).normalize();
//        return phiHat / (2.0 * math::Constants::pi * rho);
//    }
//    default:
//        throw std::logic_error("Unsupported excitation mode.");
//    }
//}
//
//void TEMCoaxial::set(const Target& elemGroup) {
//    // Reescales internal dimensions.
//    geometry::BoxR3 box = geometry::getBound(elemGroup.begin(), elemGroup.end());
//    const math::CVecR3 diagonal = box.getMax()-box.getMin();
//    if (!diagonal.isContainedInPlane(math::Constants::CartesianPlane::xy)) {
//        throw std::logic_error("Port is not contained in a XY plane");
//    }
//    const math::Real width  = box.getMax()(math::Constants::x) -
//                              box.getMin()(math::Constants::x);
//    const math::Real height = box.getMax()(math::Constants::y) -
//                              box.getMin()(math::Constants::y);
//    const math::Real averageNewRadius = (width + height)/4;
//    innerRadius_ *= averageNewRadius;
//    outerRadius_ *= averageNewRadius;
//    const math::CVecR3 averageNewOrigin = (box.getMax() + box.getMin()) / 2;
//    origin_ = averageNewOrigin;
//    //
//    Source::setTarget(elemGroup);
//}

}
}
} 
