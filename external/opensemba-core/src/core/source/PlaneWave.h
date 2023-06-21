#pragma once

#include <exception>

#include "Source.h"

namespace semba {
namespace source {

class PlaneWave : public Source {
public:
    PlaneWave(const std::unique_ptr<Magnitude::Magnitude>&,
              const Target& elem,
              const math::CVecR3& directionVector,
              const math::CVecR3& polarizationVector);
    PlaneWave(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
              const Target& elem,
              std::pair<math::Real, math::Real> directionAngles,
              std::pair<math::Real, math::Real> polarizationAngles);
    // TODO: Probably would have sense as another source type
    PlaneWave(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
              const Target& elem,
              math::Int numberOfRandomPlanewaves,
              math::Real relativeVariationOfRandomDelay);
    
    virtual ~PlaneWave() = default;

    virtual std::unique_ptr<Source> clone() const override {
        return std::make_unique<PlaneWave>(*this);
    }

    std::string getName() const override { return "PlaneWave"; };
    const math::CVecR3& getPolarization() const;
    const math::CVecR3& getDirection() const;
    
    math::Real getTheta() const;
    math::Real getPhi() const;
    math::Real getAlpha() const;
    math::Real getBeta() const;
    
    bool isRandomic() const;
    math::Int getNumberOfRandomPlanewaves() const;
    math::Real getRelativeVariationOfRandomDelay() const;

    math::CVecR3 getElectricField(const math::Real time) const;
    std::pair<math::CVecR3,math::CVecR3> getElectromagneticField(const math::Real time) const;

private:
    math::CVecR3 direction_;
    math::CVecR3 polarization_;

    bool randomic_ = false;
    math::Int numberOfRandomPlanewaves_ = 0;
    math::Real relativeVariationOfRandomDelay_ = 0.0;

    void init_(math::CVecR3 direction, math::CVecR3 polarization);
    static std::pair<math::Real,math::Real> cartesianToPolar(const math::CVecR3& vec);
    static math::CVecR3 polarToCartesian(math::Real theta, math::Real phi);
    static math::Real reduceRadians(const math::Real radianIn);
};

// TODO: Borrar namespace de error y refactor usos
namespace Error {
namespace PlaneWave {

class Error : public std::exception {
public:
    Error() {}
    virtual ~Error() throw() {}
};

class ZeroPolarization : public Error {
public:
    ZeroPolarization() {}
    virtual ~ZeroPolarization() throw() {}

    const char* what() const throw() {
        return "PlaneWave: Polarization can't be zero.";
    }
};

class ZeroMagnitude : public Error {
public:
    ZeroMagnitude() {}
    virtual ~ZeroMagnitude() throw() {}

    const char* what() const throw() {
        return "PlaneWave: W. direction can't be zero.";
    }
};

class NotPerpendicular : public Error {
public:
    NotPerpendicular() {}
    virtual ~NotPerpendicular() throw() {}

    const char* what() const throw() {
        return ("PlaneWave: W. direction is not "
                "perpendicular to polarization.");
    }
};

}
} 
}
} 

