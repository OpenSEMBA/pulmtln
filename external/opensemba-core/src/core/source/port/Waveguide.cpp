#include "Waveguide.h"
#include "core/geometry/Bound.h"

namespace semba {
namespace source {
namespace port {

Waveguide::Waveguide(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
        const Target& elem,
        const ExcitationMode excMode,
        const std::pair<size_t,size_t> mode)
:   Port(magnitude, elem) {

    excitationMode_ = excMode;
    mode_ = mode;
    // Performs checks
    //if (!geometry::getBound(elem.begin(), elem.end()).isSurface()) {
    //    throw std::logic_error("Waveport elements must be contained "
    //                           "in a coplanar geometry::Surface");
    //}

    //math::CVecR3 diagonal = geometry::getBound(elem.begin(), elem.end()).getMax() -
    //    geometry::getBound(elem.begin(), elem.end()).getMin();
    //if (!diagonal.isContainedInPlane(math::Constants::xy)) {
    //    throw std::logic_error("Waveport must be contained in plane xy.");
    //}

    //if (elem.size() == 0) {
    //    throw std::logic_error("Waveport must contain some elements.");
    //}
}

Waveguide::Waveguide(const Waveguide& rhs) :
                Port(rhs) {
    excitationMode_ = rhs.excitationMode_;
    mode_ = rhs.mode_;
}

Waveguide::ExcitationMode Waveguide::getExcitationMode() const {
    return excitationMode_;
}

std::pair<size_t,size_t> Waveguide::getMode() const {
    return mode_;
}

}
}
} 
