

#pragma once

#include "Multiport.h"

namespace semba {
namespace physicalModel {
namespace multiport {

class RLC : public virtual Multiport {
public:
    RLC(const Id idIn,
                 const std::string nameIn,
                 const Multiport::Type typeIn,
                 const math::Real resistance,
                 const math::Real inductance,
                 const math::Real capacitance);
    RLC(const RLC&);
    virtual ~RLC();

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<RLC>(*this);
    }

    virtual math::Real getR() const;
    virtual math::Real getL() const;
    virtual math::Real getC() const;
private:
    math::Real R_, L_, C_;
};

}
} 
} 

