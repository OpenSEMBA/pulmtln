#pragma once

#include "OutputRequest.h"

namespace semba {
namespace outputRequest {
    class OnSurface : public virtual OutputRequest {
    public:
        OnSurface(
            const Type& outputType,
            const Domain& domain,
            const std::string& name,
            const Target& box
        );
        virtual ~OnSurface() = default;

        std::unique_ptr<OutputRequest> clone() const override {
            return std::make_unique<OnSurface>(*this);
        }
    };

} 
} 
