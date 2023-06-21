

#include <physicalModel/surface/SIBC.h>

namespace semba {
namespace physicalModel {
namespace surface {

SIBC::SIBC(const Id id,
        const std::string& name,
        const math::MatC22& Zinfinite,
        const math::MatC22& ZLinear,
        const std::vector<PoleResidue>& poleImpedance)
:   Identifiable<Id>(id),
    PhysicalModel(name),
    ZInfinity_(Zinfinite),
    ZLinear_(ZLinear),
    poleZ_(poleImpedance) {
}

SIBC::~SIBC() {

}

std::size_t SIBC::getNumberOfPoles() const {
    return poleZ_.size();
}

} 
} 
} 
