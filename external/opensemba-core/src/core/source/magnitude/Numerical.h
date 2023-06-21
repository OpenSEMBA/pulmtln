#pragma once

#include "Magnitude.h"
#include "core/util/ProjectFile.h"

namespace semba {
namespace source {
namespace Magnitude {

class Numerical : public virtual Magnitude {
public:
    Numerical() = default;
    Numerical(const Numerical&) = default;
    Numerical(const util::ProjectFile& filename);
    Numerical(const util::ProjectFile& filename,
              const Magnitude& mag,
              const math::Real timeStep,
              const math::Real finalTime);
    
    std::unique_ptr<Magnitude> clone() const override {
        return std::make_unique<Numerical>(*this);
    }

    bool operator==(const Numerical&) const;
    math::Real evaluate(const math::Real time) const;

    util::ProjectFile getFile() const { return file; };

private:
    static const std::size_t defaultNumberOfSteps = 1000;
    util::ProjectFile file;

};

}
}
} 

