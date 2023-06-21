#pragma once

#include "core/math/matrix/Static.h"
#include "CartesianVector.h"

namespace semba {
namespace math {

class LocalAxis {
public:
    LocalAxis() = default;
    LocalAxis(CVecR3 eulerAngles, CVecR3 origin = CVecR3());
    

    MatR33 getTransformationMatrix() const;
    const CVecR3 getEulerAngles() const;
    const CVecR3 getOrigin() const;

    MatR33 convertToGlobal(const MatR33& local) const;
    CVecR3 convertToGlobal(const CVecR3& local) const;

private:
    CVecR3 eulerAngles_; // Euler angles in radians.
    CVecR3 origin_;
};

} 
} 

