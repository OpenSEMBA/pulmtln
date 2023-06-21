#include "Coordinate.h"

namespace semba {
namespace geometry {
namespace coordinate {

bool Base::operator ==(const Base& rhs) const {
    if (typeid(*this) == typeid(rhs)) {
        return true;
    }
    return false;
}

bool Base::operator !=(const Base& rhs) const {
    return !(*this == rhs);
}

} 
} 
} 
