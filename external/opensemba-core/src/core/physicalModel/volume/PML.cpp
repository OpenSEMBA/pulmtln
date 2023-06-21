

#include <physicalModel/volume/PML.h>

namespace semba {
namespace physicalModel {
namespace volume {

PML::PML(const Id id,
         const std::string& name,
         const math::LocalAxis orientation)
:   Identifiable<Id>(id),
    PhysicalModel(name),
    orientation_(orientation) {
}

PML::PML(const PML& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs),
    orientation_(rhs.orientation_) {
}

const math::LocalAxis PML::getOrientation() const {
    return orientation_;
}

const math::CVecR3 PML::getGlobalZAxis() const {
    math::CVecR3 localZ(0.0,0.0,1.0);
    math::CVecR3 res = getOrientation().convertToGlobal(localZ)
            - getOrientation().getOrigin();
    return res;
}

}
} 
} 

