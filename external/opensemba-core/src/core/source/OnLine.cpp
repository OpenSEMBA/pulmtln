#include "OnLine.h"

namespace semba {
namespace source {

OnLine::OnLine(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
               const Target& elem,
               const Type& sourceType,
               const Hardness& sourceHardness)
:   Source(magnitude, elem)
{
    type_ = sourceType;
    hardness_ = sourceHardness;
}

}
} 
