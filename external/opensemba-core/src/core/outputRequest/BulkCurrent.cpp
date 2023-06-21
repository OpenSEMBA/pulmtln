

#include "BulkCurrent.h"

namespace semba {
namespace outputRequest {

BulkCurrent::BulkCurrent(
        const Domain& domain,
        const std::string& name,
        const Target& elem,
        const math::Constants::CartesianAxis& dir,
        const math::UInt& skip) :
    OutputRequest(Type::bulkCurrentElectric, name, domain, elem),
    dir_(dir),
    skip_(skip)
{}

math::Constants::CartesianAxis BulkCurrent::getDir() const {
    return dir_;
}

math::UInt BulkCurrent::getSkip() const {
    return skip_;
}

} 
} 
