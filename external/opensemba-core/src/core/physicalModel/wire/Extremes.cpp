#include "core/physicalModel/wire/Extremes.h"

namespace semba {
namespace physicalModel {
namespace wire {

Extremes::Extremes(const std::string& name,
                   const Wire& wire,
                   const multiport::Multiport* extremeL,
                   const multiport::Multiport* extremeR)
:   Wire(wire)
{
    setName(name);
    setId(wire.getId());
    if (extremeL != nullptr) {
        extreme_[0].reset(extremeL->clone().release()->castTo<multiport::Multiport>());
    }
    if (extremeR != nullptr) {
        extreme_[1].reset(extremeR->clone().release()->castTo<multiport::Multiport>());
    }
}

Extremes::Extremes(const Extremes& rhs) : 
    Identifiable<Id>(rhs),
    PhysicalModel(rhs),
    Wire(rhs)
{
    if (rhs.extreme_[0] != nullptr) {
        extreme_[0].reset(rhs.extreme_[0]->clone().release()->castTo<multiport::Multiport>());
    }
    if (rhs.extreme_[1] != nullptr) {
        extreme_[1].reset(rhs.extreme_[1]->clone().release()->castTo<multiport::Multiport>());
    }
}

void Extremes::setExtreme(const std::size_t i,
                          const multiport::Multiport* extreme) {
    if (extreme == nullptr) {
        extreme_[i].reset();
    }
    else {
        extreme_[i].reset(extreme->clone().release()->castTo<multiport::Multiport>());
    }
}

void Extremes::swapExtremes() {
    std::swap(extreme_[0], extreme_[1]);
}

}
} 
} 
