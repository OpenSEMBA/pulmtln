#include "Dipole.h"

namespace semba {
namespace source {

Dipole::Dipole(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
               const Target& elem,
               math::Real length,
               math::CVecR3 orientation,
               math::CVecR3 position)
:   Source(magnitude, elem)
{
    length_ = length;
    orientation_ = orientation;
    position_ = position;
}


}
} 
