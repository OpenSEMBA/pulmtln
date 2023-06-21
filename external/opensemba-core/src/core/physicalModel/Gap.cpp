#include "core/physicalModel/Gap.h"

namespace semba {
namespace physicalModel {

Gap::Gap(const Id id, const std::string name, const math::Real width) :   
    Identifiable<Id>(id),
    PhysicalModel(name) 
{
    width_ = width;
}

Gap::Gap(const Gap& rhs) :   
    Identifiable<Id>(rhs),
    PhysicalModel(rhs) 
{
    width_ = rhs.width_;
}

math::Real Gap::getWidth() const 
{
    return width_;
}

} 
} 
