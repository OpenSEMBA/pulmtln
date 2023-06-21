#include "PlaneWave.h"

namespace semba {
namespace source {

PlaneWave::PlaneWave(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
                     const Target& elem,
                     const math::CVecR3& direction,
                     const math::CVecR3& polarization) :   
    Source(magnitude, elem)
{
    init_(direction, polarization);
}

PlaneWave::PlaneWave(
        const std::unique_ptr<Magnitude::Magnitude>& magnitude,
        const Target& elem,
        std::pair<math::Real, math::Real> directionAngles,
        std::pair<math::Real, math::Real> polarizationAngles) :   
    Source(magnitude, elem) 
{

    math::Real theta = directionAngles.first;
    math::Real phi   = directionAngles.second;
    math::Real alpha_ = polarizationAngles.first;
    math::Real beta_  = polarizationAngles.second;

    math::CVecR3 dirVec = polarToCartesian(theta, phi);
    math::CVecR3 polVec = polarToCartesian(alpha_, beta_);
    init_(dirVec, polVec);
}

PlaneWave::PlaneWave(
        const std::unique_ptr<Magnitude::Magnitude>& magnitude,
        const Target& elem,
        math::Int numberOfRandomPlanewaves,
        math::Real relativeVariationOfRandomDelay) :   
    Source(magnitude, elem) 
{
    randomic_ = true;
    numberOfRandomPlanewaves_ = numberOfRandomPlanewaves;
    relativeVariationOfRandomDelay_ = relativeVariationOfRandomDelay;
}

const math::CVecR3& PlaneWave::getPolarization() const {
    return polarization_;
}

const math::CVecR3& PlaneWave::getDirection() const {
    return direction_;
}

math::Real PlaneWave::getTheta() const {
    return cartesianToPolar (direction_).first;
}

math::Real PlaneWave::getPhi() const {
    return cartesianToPolar (direction_).second;
}

math::Real PlaneWave::getAlpha() const {
    return cartesianToPolar (polarization_).first;
}

math::Real PlaneWave::getBeta() const {
    return cartesianToPolar (polarization_).second;
}

bool PlaneWave::isRandomic() const {
    return randomic_;
}

math::Int PlaneWave::getNumberOfRandomPlanewaves() const {
    return numberOfRandomPlanewaves_;
}

math::Real PlaneWave::getRelativeVariationOfRandomDelay() const {
    return relativeVariationOfRandomDelay_;
}

math::CVecR3 PlaneWave::getElectricField(const math::Real time) const {
    math::CVecR3 res = polarization_ * getMagnitude()->evaluate(time);
    return res;
}

std::pair<math::CVecR3, math::CVecR3>
PlaneWave::getElectromagneticField(const math::Real time) const {
    math::CVecR3 electric = getElectricField(time);
    math::CVecR3 magnetic = (direction_ ^ electric) *
                            math::Constants::VACUUM_ADMITANCE;
    return std::pair<math::CVecR3,math::CVecR3>(electric, magnetic);
}

std::pair<math::Real,math::Real> PlaneWave::cartesianToPolar(
        const math::CVecR3& v) {
    if (v.norm() == 0.0) {
        return std::pair<math::Real,math::Real>(0.0, 0.0);
    }
    math::Real theta, phi;
    theta = std::acos(v(math::Constants::z) / v.norm());
    if (v(math::Constants::x) == 0.0) {
        if (v(math::Constants::y) == 0.0) {
            phi = 0.0;
        } else if (v(math::Constants::y) > 0.0) {
            phi = math::Constants::pi_2;
        } else {
            phi = 3.0 * math::Constants::pi_2;
        }
    } else {
        phi = std::atan2(v(math::Constants::y), v(math::Constants::x));
    }
    return std::pair<math::Real,math::Real>(theta,phi);
}

void PlaneWave::init_(math::CVecR3 direction, math::CVecR3 polarization) {
    direction_ = direction;
    polarization_ = polarization;

    randomic_ = false;
    numberOfRandomPlanewaves_ = 0;
    relativeVariationOfRandomDelay_ = 0.0;

    if (polarization_.norm() == 0) {
        throw Error::PlaneWave::ZeroPolarization();
    }
    if (direction_.norm() == 0) {
        throw Error::PlaneWave::ZeroMagnitude();
    }
    math::Real dotProd = direction.dot(polarization);
    if (math::notEqual(dotProd, 0.0)) {
        throw Error::PlaneWave::NotPerpendicular();
    }
}

math::CVecR3 PlaneWave::polarToCartesian(math::Real theta, math::Real phi) {
    return math::CVecR3(
            std::sin(theta)*std::cos(phi),
            std::sin(theta)*std::sin(phi),
            std::cos(theta));
}

math::Real PlaneWave::reduceRadians(const math::Real radianIn) {
    math::Real nVueltas, nVueltasComp, radianOut, Val2Pi;
    Val2Pi = (math::Real) 2.0 * (math::Real) acos((math::Real) 0.0);
    nVueltas = radianIn/(Val2Pi);
    nVueltasComp = (math::Real) floor(nVueltas);
    radianOut = radianIn - nVueltasComp*Val2Pi;
    return  radianOut;
}

}
} 
