#pragma once

#include "Port.h"

namespace semba {
namespace source {
namespace port {

class Waveguide : public Port {
public:
    enum class ExcitationMode {
        TE,
        TM
    };

    Waveguide(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
              const Target& elem,
              const ExcitationMode excMode,
              const std::pair<size_t,size_t> mode);
    Waveguide(const Waveguide& rhs);
    virtual ~Waveguide() = default;

    //virtual std::unique_ptr<Source> clone() const override {
    //    return std::make_unique<Waveguide>(*this);
    //}

    std::string getName() const { return "Waveguide"; };

    ExcitationMode getExcitationMode() const;
    std::pair<size_t, size_t> getMode() const;

private:
    ExcitationMode excitationMode_;
    std::pair<size_t,size_t> mode_;
};

}
}
} 

