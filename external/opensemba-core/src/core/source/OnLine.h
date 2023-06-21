#pragma once

#include "Source.h"

namespace semba {
namespace source {

class OnLine : public Source {
public:
    enum class Type {
        electric = 1,
        magnetic = 2
    };
    enum class Hardness {
        hard = 1,
        soft = 2
    };
    
    OnLine(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
           const Target& elem,
           const Type& sourceType,
           const Hardness& sourceHardness);
    virtual ~OnLine() = default;
        
    virtual std::unique_ptr<Source> clone() const override {
        return std::make_unique<OnLine>(*this);
    }

    std::string getName() const override { return "OnLine"; };
    Type getType() const { return type_; };
    Hardness getHardness() const { return hardness_; };

private:
    Type type_;
    Hardness hardness_;
};

}
} 

