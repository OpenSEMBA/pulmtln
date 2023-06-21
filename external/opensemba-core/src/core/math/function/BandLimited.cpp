#include "BandLimited.h"

#include <cmath>
#include <iostream>

namespace semba {
namespace math {
namespace function {

BandLimited::BandLimited(
    const Real mininumFrequency, 
    const Real maximumFrequency) : 
    minimumFrequency_(mininumFrequency), 
    maximumFrequency_(maximumFrequency) {
    if (mininumFrequency <= 0.0 || maximumFrequency <= 0.0) {
        throw std::logic_error(
            "Band limited signals must have positive frequencies"
        );
    }
    if (mininumFrequency >= maximumFrequency) {
        throw std::logic_error(
            "In band limited signals, min. frequency must be smaller max. frequency"
        );
    }
}

Real BandLimited::operator ()(const Real& t) const {
    static const Real pi = Constants::pi;
                
    const Real carrier =
        (maximumFrequency_ + minimumFrequency_) / 2.0;
    const Real bTuned = 
        (maximumFrequency_ - minimumFrequency_) / 2.0;
    const Real delay = 20.0 / minimumFrequency_;
    const Real spread = 10.0 / maximumFrequency_;
    Real tD = t - delay;
    if (tD == 0.0) {
        tD += std::numeric_limits<Real>::epsilon();
    }
                
    Real res = std::sin(2.0 * pi* bTuned * tD) / (2 * pi*bTuned*tD) *
        exp(- std::pow(tD,2) / std::pow(spread,2) / 2.0) *
            cos(2.0 * pi * carrier * tD);
    return res;
}

bool BandLimited::operator==(const Base& rhs) const {
    if (typeid(*this) != typeid(rhs)) {
        return false;
    }
    const BandLimited* rhsPtr = dynamic_cast<const BandLimited*>(&rhs);
    return ((this->maximumFrequency_ == rhsPtr->maximumFrequency_) &&
        (this->minimumFrequency_ == rhsPtr->minimumFrequency_)
    );
}


} 
} 
} 
