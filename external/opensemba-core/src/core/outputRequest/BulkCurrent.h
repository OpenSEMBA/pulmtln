#pragma once

#include "OutputRequest.h"

namespace semba {
namespace outputRequest {

class BulkCurrent : public virtual OutputRequest {
public:
    BulkCurrent(const Domain& domain,
                const std::string& name,
                const Target& box,
                const math::Constants::CartesianAxis& dir,
                const math::UInt& skip);
    virtual ~BulkCurrent() = default;

    std::unique_ptr<OutputRequest> clone() const override {
        return std::make_unique<BulkCurrent>(*this);
    }

    math::Constants::CartesianAxis getDir() const;
    math::UInt                     getSkip() const;

private:
    math::Constants::CartesianAxis dir_;
    math::UInt                     skip_;
};

} 
} 

