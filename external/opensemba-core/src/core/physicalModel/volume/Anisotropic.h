#pragma once

#include "Volume.h"
#include "core/math/LocalAxis.h"

namespace semba {
namespace physicalModel {
namespace volume {

class Anisotropic : public virtual Volume {
public:
    enum class Model {
        crystal,
        ferrite
    };

    Anisotropic(const math::LocalAxis& local);
    Anisotropic(const Anisotropic& rhs);
    virtual ~Anisotropic() = default;

    math::LocalAxis getLocalAxe() const;
    virtual math::MatR33 getRelPermittivityMatR() const = 0;
    virtual math::MatR33 getRelPermeabilityMatR() const = 0;
    virtual math::MatR33 getElectricConductivityMat() const = 0;
    virtual math::MatR33 getMagneticConductivityMat() const = 0;

private:
    math::LocalAxis localAxe_;
};


// Described in: https://courses.cit.cornell.edu/ece303/Lectures/lecture17.pdf
class AnisotropicCrystal: public Anisotropic {
public:
    AnisotropicCrystal(
            const Id matId,
            const std::string& name,
            const math::LocalAxis& local,
            const math::CVecR3& principalAxesRelativePermittivity,
            const math::Real relativePermeability);
    AnisotropicCrystal(const AnisotropicCrystal&);
    virtual ~AnisotropicCrystal() = default;

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<AnisotropicCrystal>(*this);
    }

    const math::CVecR3 getPrincipalAxesRelativePermittivity() const;
    math::Real getRelativePermeability() const;

    math::MatR33 getRelPermittivityMatR() const override;
    math::MatR33 getRelPermeabilityMatR() const override;
    math::MatR33 getElectricConductivityMat() const override;
    math::MatR33 getMagneticConductivityMat() const override;
private:
    math::CVecR3 principalAxesRelativePermittivity_;
    math::Real relativePermeability_;
};

class AnisotropicFerrite: public Anisotropic {
public:
    AnisotropicFerrite(const Id matId,
                             const std::string& name,
                             const math::LocalAxis& local,
                             const math::Real kappa,
                             const math::Real relativePermeability,
                             const math::Real relativePermittivity);
    AnisotropicFerrite(const AnisotropicFerrite&);
    virtual ~AnisotropicFerrite() = default; 

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<AnisotropicFerrite>(*this);
    }

    math::MatR33 getRelPermittivityMatR() const override;

    math::MatR33 getRelPermeabilityMatR() const override; // math::Real part.
    math::MatR33 getRelPermeabilityMatI() const; // Imaginary part.

    math::MatR33 getElectricConductivityMat() const override;
    math::MatR33 getMagneticConductivityMat() const override;
private:
    math::Real kappa_;
    math::Real relativePermeability_;
    math::Real relativePermittivity_;
};

}
} 
} 

