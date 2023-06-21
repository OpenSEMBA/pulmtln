#include "OnLayer.h"

namespace semba {
namespace outputRequest {

    OnLayer::OnLayer(
        const Type& outputType,
        const Domain& domain,
        const std::string& name,
        const Target& elem
    )
        : OutputRequest(outputType, name, domain, elem)
    {}

} 
} 
