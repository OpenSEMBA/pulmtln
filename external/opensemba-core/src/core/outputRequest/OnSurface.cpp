#include "OnSurface.h"

namespace semba {
namespace outputRequest {

OnSurface::OnSurface(
    const Type& outputType,
    const Domain& domain,
    const std::string& name,
    const Target& elem
) : OutputRequest(outputType, name, domain, elem)
{}

} 
} 
