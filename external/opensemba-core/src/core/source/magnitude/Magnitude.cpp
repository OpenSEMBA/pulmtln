#include "Magnitude.h"

namespace semba {
namespace source {
namespace Magnitude {

Magnitude::Magnitude(math::FunctionRR* mathFunction) {
    mathFunction_ = mathFunction;
}

Magnitude::Magnitude(const Magnitude& rhs) {
    mathFunction_ = dynamic_cast<math::FunctionRR*>(rhs.mathFunction_->clone());
}

Magnitude& Magnitude::operator=(const Magnitude& rhs) {
    mathFunction_ = dynamic_cast<math::FunctionRR*>(rhs.mathFunction_->clone());
    return *this;
}

bool Magnitude::operator ==(const Magnitude& rhs) const {
    return *mathFunction_ == *rhs.mathFunction_;
}

math::Real Magnitude::evaluate(const math::Real time) const {
    return mathFunction_->operator()(time);
}

}
}
} 
