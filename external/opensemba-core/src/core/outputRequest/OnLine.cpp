#include "OnLine.h"

namespace semba {
namespace outputRequest {
    OnLine::OnLine(
        const Type& outputType,
        const Domain& domain,
        const std::string& name,
        const Target& elem
    )
        : OutputRequest(outputType, name, domain, elem)
    {}

} 
} 
