#include "WaveguideRectangular.h"
#include "core/geometry/Bound.h"

namespace semba {
namespace source {
namespace port {

WaveguideRectangular::WaveguideRectangular(
        const std::unique_ptr<Magnitude::Magnitude>& magn,
        const Target& elem,
        const ExcitationMode excMode,
        const std::pair<size_t,size_t> mode): 
    Waveguide(magn, elem, excMode, mode) 
{
    if (mode.first == 0 && mode.second == 0) {
        throw std::logic_error("At least one mode must be non-zero.");
    }
}

std::string WaveguideRectangular::getName() const {
    return "Rectangular_waveguide_port";
}

//math::CVecR3 WaveguideRectangular::getWeight(
//        const math::CVecR3& pos) const {
//    // Return normalized weights for electric field components.
//    static const math::Real pi = acos(-1.0);
//    math::CVecR3 res;
//    math::CVecR3 rPos = pos - getOrigin();
//    const math::Real m = pi * getMode().first / getWidth();
//    const math::Real n = pi * getMode().second / getHeight();
//    math::Real normFactor = m;
//    if (n > m) {
//        normFactor = n;
//    }
//    //const math::Real betaC = sqrt(pow(m,2) + pow(n,2));
//    if (getExcitationMode() == Waveguide::ExcitationMode::TE) {
//        res(math::Constants::x) =   n * cos(m * rPos(math::Constants::x)) *
//                                        sin(n * rPos(math::Constants::y)) /
//                                        normFactor;
//        res(math::Constants::y) =   m * sin(m * rPos(math::Constants::x)) *
//                                        cos(n * rPos(math::Constants::y)) /
//                                        normFactor;
//        res(math::Constants::z) = (math::Real) 0.0;
//    } else {
//        res(math::Constants::x) = - m * cos(m * rPos(math::Constants::x)) *
//                                        sin(n * rPos(math::Constants::y)) /
//                                        normFactor;
//        res(math::Constants::y) = - m * sin(m * rPos(math::Constants::x)) *
//                                        cos(n * rPos(math::Constants::y)) /
//                                        normFactor;
//        res(math::Constants::z) = (math::Real) 0.0;
//    }
//    return res;
//}

//math::Real WaveguideRectangular::getWidth() const {
//    math::CVecR3 origin = getOrigin();
//    math::CVecR3 max = box_.getMax();
//    return max(math::Constants::x) - origin(math::Constants::x);
//}
//
//math::Real WaveguideRectangular::getHeight() const {
//    math::CVecR3 origin = getOrigin();
//    math::CVecR3 max = box_.getMax();
//    return max(math::Constants::y) - origin(math::Constants::y);
//}
//
//void WaveguideRectangular::set(
//    const Target& constGroupElems) {
//    Waveguide::setTarget(constGroupElems);
//    box_ = geometry::getBound(constGroupElems.begin(), constGroupElems.end());
//}
//
//math::CVecR3 WaveguideRectangular::getOrigin() const {
//    return box_.getMin();
//}

}
}
} 
