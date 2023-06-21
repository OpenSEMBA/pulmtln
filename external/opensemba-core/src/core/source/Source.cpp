#include "Source.h"

#include <iostream>

namespace semba {
namespace source {

Source::Source(
    const std::unique_ptr<Magnitude::Magnitude>& magnitude, 
    const Target& target) :
    magnitude_(magnitude->clone()),
    target_(target)
{
}

Source::Source(const Source& rhs) 
{
    target_ = rhs.target_;
    if (rhs.magnitude_ != nullptr) {
        magnitude_ = rhs.magnitude_->clone();
    }
}

Source& Source::operator=(const Source& rhs)
{
    target_ = rhs.target_;
    if (rhs.magnitude_ != nullptr) {
        magnitude_ = rhs.magnitude_->clone();
    }

    return *this;
}


void Source::convertToNumerical(
    const util::ProjectFile& file,
    const math::Real step,
    const math::Real finalTime) 
{
    if(magnitude_->is<Magnitude::Numerical>()) {
        return;
    }
    auto newMagnitude = std::make_unique<Magnitude::Numerical>(file, *magnitude_, step, finalTime);
    magnitude_ = std::move(newMagnitude);
}

Magnitude::Numerical Source::exportToFile(const util::ProjectFile& file,
                                         const math::Real step,
                                         const math::Real finalTime) const {
    if(magnitude_->is<Magnitude::Numerical>()) {
        return Magnitude::Numerical(*magnitude_->castTo<Magnitude::Numerical>());
    }
    return Magnitude::Numerical(file, *magnitude_, step, finalTime);
}

}
} 
