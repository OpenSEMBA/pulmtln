#pragma once

#include "core/source/Source.h"
#include "core/physicalModel/Bound.h"

namespace semba {
namespace source {
namespace port {

typedef std::array<std::array<const physicalModel::Bound*,2>,3> Bound3;

class Port : public Source {
public:
    Port(const std::unique_ptr<Magnitude::Magnitude>& magnitude, 
         const Target& elem);
    virtual ~Port() = default;
    
    std::string getName() const { return "Port"; };
    
    //virtual math::CVecR3 getOrigin() const = 0;
    //virtual math::CVecR3 getWeight(const math::CVecR3& pos) const = 0;
};

}
}
} 

