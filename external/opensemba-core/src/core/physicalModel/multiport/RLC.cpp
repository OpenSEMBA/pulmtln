

#include <physicalModel/multiport/RLC.h>

namespace semba {
namespace physicalModel {
namespace multiport {

RLC::RLC(const Id id,
                           const std::string name,
                           const Multiport::Type type,
                           const math::Real resistance,
                           const math::Real inductance,
                           const math::Real capacitance)
:   Identifiable<Id>(id),
    PhysicalModel(name) {
    type_ = type;
    R_ = resistance;
    L_ = inductance;
    C_ = capacitance;
}

RLC::RLC(const RLC& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs) {
    type_ = rhs.type_;
    R_ = rhs.R_;
    L_ = rhs.L_;
    C_ = rhs.C_;
}

RLC::~RLC() {

}

math::Real RLC::getR() const {
    return R_;
}

math::Real RLC::getL() const {
    return L_;
}

math::Real RLC::getC() const {
    return C_;
}


}
} 
} 
