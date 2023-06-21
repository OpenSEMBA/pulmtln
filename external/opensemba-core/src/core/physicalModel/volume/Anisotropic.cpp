#include "Anisotropic.h"

namespace semba {
namespace physicalModel {
namespace volume {

Anisotropic::Anisotropic(const math::LocalAxis& localAxe) {
    localAxe_ = localAxe;
}

Anisotropic::Anisotropic(const Anisotropic& rhs) {
    localAxe_ = rhs.localAxe_;
}

math::LocalAxis Anisotropic::getLocalAxe() const {
    return localAxe_;
}

AnisotropicCrystal::AnisotropicCrystal(
        const Id matId,
        const std::string& name,
        const math::LocalAxis& local,
        const math::CVecR3& principalAxesRelativePermittivity,
        const math::Real relativePermeability)
:   Identifiable<Id>(matId),
    PhysicalModel(name),
    Anisotropic(local) {
    principalAxesRelativePermittivity_ = principalAxesRelativePermittivity;
    relativePermeability_ = relativePermeability;
}

AnisotropicCrystal::AnisotropicCrystal(
    const AnisotropicCrystal& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs),
    Anisotropic(rhs) {
    principalAxesRelativePermittivity_ =
        rhs.principalAxesRelativePermittivity_;
    relativePermeability_ = rhs.relativePermeability_;
}

const math::CVecR3
    AnisotropicCrystal::getPrincipalAxesRelativePermittivity() const {
    return principalAxesRelativePermittivity_;
}

math::Real AnisotropicCrystal::getRelativePermeability() const {
    return relativePermeability_;
}


math::MatR33 AnisotropicCrystal::getRelPermittivityMatR() const {
    math::MatR33 local;
    local.setInDiagonal(principalAxesRelativePermittivity_);
    return getLocalAxe().convertToGlobal(local);
}

math::MatR33 AnisotropicCrystal::getRelPermeabilityMatR() const {
    math::MatR33 local;
    local.setInDiagonal(math::CVecR3(relativePermeability_));
    return getLocalAxe().convertToGlobal(local);
}

math::MatR33 AnisotropicCrystal::getElectricConductivityMat() const {
    return math::MatR33();
}

math::MatR33 AnisotropicCrystal::getMagneticConductivityMat() const {
    return math::MatR33();
}


AnisotropicFerrite::AnisotropicFerrite(
        const Id matId,
        const std::string& name,
        const math::LocalAxis& local,
        const math::Real kappa,
        const math::Real relativePermeability,
        const math::Real relativePermittivity)
:   Identifiable<Id>(matId),
    PhysicalModel(name),
    Anisotropic(local) {
    kappa_ = kappa;
    relativePermeability_ = relativePermeability;
    relativePermittivity_ = relativePermittivity;
}

AnisotropicFerrite::AnisotropicFerrite(
    const AnisotropicFerrite& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs),
    Anisotropic(rhs) {
    kappa_ = rhs.kappa_;
    relativePermeability_ = rhs.relativePermeability_;
    relativePermittivity_ = rhs.relativePermittivity_;
}

math::MatR33 AnisotropicFerrite::getRelPermittivityMatR() const {
    return math::MatR33().setInDiagonal(math::CVecR3(relativePermittivity_));
}

math::MatR33 AnisotropicFerrite::getRelPermeabilityMatR() const {
    math::MatR33 local;
    math::CVecR3 principalAxis(relativePermeability_,
                               relativePermeability_, 1.0);
    local.setInDiagonal(principalAxis);
    return getLocalAxe().convertToGlobal(local);
}

math::MatR33 AnisotropicFerrite::getRelPermeabilityMatI() const {
    math::MatR33 local;
    local(0,1) =   kappa_;
    local(1,0) = - kappa_;
    return getLocalAxe().convertToGlobal(local);
}

math::MatR33 AnisotropicFerrite::getElectricConductivityMat() const {
    return math::MatR33();
}

math::MatR33 AnisotropicFerrite::getMagneticConductivityMat() const {
    return math::MatR33();
}


}
} 
} 
