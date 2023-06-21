#pragma once

#include "Wire.h"
#include "core/physicalModel/multiport/Multiport.h"

namespace semba {
namespace physicalModel {
namespace wire {

class Extremes : public virtual Wire  {
public:
    Extremes(const std::string&,
             const Wire& wire,
             const multiport::Multiport* extremeL,
             const multiport::Multiport* extremeR);
    Extremes(const Extremes& rhs);
    virtual ~Extremes() = default;

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Extremes>(*this);
    }

    const multiport::Multiport *getExtreme(const std::size_t i) const {
        return extreme_[i].get();
    }

    void setExtreme(const std::size_t i, const multiport::Multiport* extreme);
    void swapExtremes();

private:
    // TODO: Remove plain array, maybe GroupIdentifiableUnique / Physical Model Group?
    std::unique_ptr<const multiport::Multiport> extreme_[2];
};

}
} 
} 

