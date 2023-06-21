#pragma once

#include "PhysicalModel.h"

namespace semba {
namespace physicalModel {

class Gap : public virtual PhysicalModel {
public:
    Gap(const Id id, const std::string name, const math::Real width);
    Gap(const Gap&);
    virtual ~Gap() = default;

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Gap>(*this);
    }

    math::Real getWidth() const;

private:
    math::Real width_;
};

} 
} 

