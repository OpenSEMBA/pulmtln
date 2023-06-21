#pragma once

#include <limits>

#include "Volume.h"

namespace semba {
namespace physicalModel {
namespace volume {

class Classic : public virtual Volume {
public:
    Classic(const Id matId,
                  const std::string& name,
                  const math::Real relativePermittivity = 1.0,
                  const math::Real relativePermeability = 1.0,
                  const math::Real electricConductivity = 0.0,
                  const math::Real magneticConductivity = 0.0);
    Classic(const Classic&);
    virtual ~Classic();

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Classic>(*this);
    }

    math::Real getRelativePermittivity() const;
    math::Real getPermittivity() const;
    math::Real getRelativePermeability() const;
    math::Real getPermeability() const;
    math::Real getImpedance() const;
    math::Real getAdmitance() const;
    math::Real getElectricConductivity() const;
    math::Real getMagneticConductivity() const;
    bool isVacuum() const;

private:
    math::Real rEps_; // Rel. permittivity @ infte. freq.
    math::Real rMu_; // Rel. permeability @ infte. freq.
    math::Real electricConductivity_;
    math::Real magneticConudctivity_;
};

}
} 
} 

