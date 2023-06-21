#pragma once

#include "core/math/matrix/Static.h"
#include "core/util/ProjectFile.h"

#include "Surface.h"

namespace semba {
namespace physicalModel {
namespace surface {

class SIBC : public virtual Surface {
public:
    typedef std::pair<std::complex<math::Real>, math::MatC22> PoleResidue;

    SIBC(const Id id,
            const std::string& name,
            const math::MatC22& Zinfinite,
            const math::MatC22& Zstatic,
            const std::vector<PoleResidue>& poleImpedance);
    virtual ~SIBC();

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<SIBC>(*this);
    }

    virtual std::size_t getNumberOfPoles() const;

    std::vector<PoleResidue> getPoleZ() const {return poleZ_;}
    math::MatC22 getZInfinity() const {return ZInfinity_;}
    math::MatC22 getZLinear() const {return ZLinear_;}

private:
    math::MatC22 ZInfinity_, ZLinear_;
    std::vector<PoleResidue> poleZ_;
};

} 
} 
} 

