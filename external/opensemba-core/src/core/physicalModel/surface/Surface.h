#pragma once

#include "core/physicalModel/PhysicalModel.h"

namespace semba {
namespace physicalModel {
namespace surface {

class Surface : public virtual PhysicalModel {
public:
    Surface() = default;
    virtual ~Surface() = default;
};

} 
} 
} 

