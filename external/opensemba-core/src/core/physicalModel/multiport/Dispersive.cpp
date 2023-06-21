

#include "Dispersive.h"

namespace semba {
namespace physicalModel {
namespace multiport {

Dispersive::Dispersive(const Id id,
                       const std::string name,
                       const std::string filename)
:   Multiport(id, name, Multiport::Type::dispersive) {
    filename_ = filename;
}

Dispersive::Dispersive(const Dispersive& rhs)
:   Identifiable<Id>(rhs),
    PhysicalModel(rhs) {
    filename_ = rhs.filename_;
    type_ = rhs.type_;
}


std::string Dispersive::getFilename() const {
    return filename_;
}


} 
} 
} 
