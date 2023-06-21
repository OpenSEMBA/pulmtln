#pragma once

#include "core/geometry/Box.h"
#include "core/math/LocalAxis.h"

#include "Volume.h"

namespace semba {
namespace physicalModel {
namespace volume {

class PML : public virtual Volume {
public:
    PML(const Id id, const std::string& name, const math::LocalAxis orientation);
    PML(const PML& rhs);
    virtual ~PML() = default;
  
    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<PML>(*this);
    }

    const math::LocalAxis getOrientation() const;
    const math::CVecR3 getGlobalZAxis() const;

private:
    const math::LocalAxis orientation_;
};

}
} 
} 

