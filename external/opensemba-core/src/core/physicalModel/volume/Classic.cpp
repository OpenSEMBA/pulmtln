#include <physicalModel/volume/Classic.h>
#include <cmath>

namespace semba {
namespace physicalModel {
namespace volume {

Classic::Classic(const Id matId,
                 const std::string& name,
                 const math::Real relativePermittivity,
                 const math::Real relativePermeability,
                 const math::Real electricConductivity,
                 const math::Real magneticConductivity)
:   Identifiable<Id>(matId), 
    PhysicalModel(name) {
    rEps_ = relativePermittivity;
    rMu_ = relativePermeability;
    electricConductivity_ = electricConductivity;
    magneticConudctivity_ = magneticConductivity;
    if (relativePermittivity< 1.0 ||
        relativePermeability < 1.0 ||
        electricConductivity < 0.0 ||
        magneticConductivity < 0.0
     ) {
        std::cerr << std::endl << "WARNING: "
            << "Material " << matId << ": " << name << " has wrong "
            << "permittivity, permeability, or conductivity." << std::endl;
    }
}

Classic::Classic(const Classic& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs) {
    rEps_ = rhs.rEps_;
    rMu_ = rhs.rMu_;
    electricConductivity_ = rhs.electricConductivity_;
    magneticConudctivity_ = rhs.magneticConudctivity_;
}

Classic::~Classic() {

}

math::Real Classic::getImpedance() const {
    if (rEps_ <= 0.0) {
        return std::numeric_limits<math::Real>::infinity();
    }
    return std::sqrt((rMu_ * math::Constants::mu0) /
                     (rEps_ * math::Constants::eps0));
}

math::Real Classic::getAdmitance() const {
    if (rMu_ <= 0.0) {
        return std::numeric_limits<math::Real>::infinity();
    }
    return (1.0 / getImpedance());
}

math::Real Classic::getRelativePermittivity() const {
    return rEps_;
}

math::Real Classic::getRelativePermeability() const {
    return rMu_;
}

math::Real Classic::getPermittivity() const {
    return (rEps_ * math::Constants::eps0);
}

math::Real Classic::getPermeability() const {
    return (rMu_ * math::Constants::mu0);
}

math::Real Classic::getElectricConductivity() const {
    return electricConductivity_;
}

math::Real Classic::getMagneticConductivity() const {
    return magneticConudctivity_;
}

bool Classic::isVacuum() const {
    return (rEps_ == 1.0
            && rMu_ == 1.0
            && electricConductivity_ == 0.0
            && magneticConudctivity_ == 0.0);
}

}
} 
} 
