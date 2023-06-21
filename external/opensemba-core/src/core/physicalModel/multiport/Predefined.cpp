#include "core/physicalModel/multiport/Predefined.h"

namespace semba {
namespace physicalModel {
namespace multiport {

Predefined::Predefined(const Id id,
        const std::string name,
        const Multiport::Type type)
:   Identifiable<Id>(id),
    PhysicalModel(name) {
    type_ = type;
}

Predefined::Predefined(const Predefined& rhs)
: Identifiable<Id>(rhs),
  PhysicalModel(rhs) {
    type_ = rhs.type_;
}



}
} 
} 
