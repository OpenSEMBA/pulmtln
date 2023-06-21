#include <iomanip>

#include "core/math/function/LinearInterpolation.h"
#include "core/math/Real.h"

#include "Numerical.h"

namespace semba {
namespace source {
namespace Magnitude {

Numerical::Numerical(const util::ProjectFile& fileIn) :   
    Magnitude(new math::function::LinearInterpolation<math::Real, math::Real>(fileIn)),
    file(fileIn)
{}

Numerical::Numerical(const util::ProjectFile& fileIn,
                     const Magnitude& mag,
                     const math::Real timeStep,
                     const math::Real finalTime) 
{
    file = fileIn;
    if(mag.is<Numerical>()) {
        operator=(*mag.castTo<Numerical>());
        return;
    }
    std::size_t nSteps;
    if (timeStep != 0.0) {
        nSteps = (std::size_t)std::abs(finalTime / timeStep);
    } else {
        nSteps = defaultNumberOfSteps;
        std::cerr << "WARNING @ Numerical: "
                  << "Attempting to build a "
                  << "numerical magnitude with a 0.0 step."
                  << "Using default number of steps instead: " << nSteps
                  << std::endl;
    }
    std::ofstream out;
    out << std::scientific;
    out.precision(10);
    out.open(file.c_str());

    math::Real time = 0.0;
    for (std::size_t i = 0; i < nSteps; i++) {
        // Determines if neigh values are aligned with current.
        std::vector<std::pair<math::Real,math::Real>> preAndPost;
        const math::Real tPre = time - timeStep;
        const math::Real tPost = time + timeStep;
        preAndPost.push_back(std::pair<math::Real,math::Real>(tPre, mag.evaluate(tPre)));
        preAndPost.push_back(std::pair<math::Real,math::Real>(tPost, mag.evaluate(tPost)));
        const math::Real interpolated =
            math::function::LinearInterpolation<math::Real,math::Real>(
                preAndPost)(time);
        const math::Real current = mag.evaluate(time);
        bool isAligned = math::equal(current, interpolated,
                0.0, std::numeric_limits<math::Real>::epsilon());
        if (!isAligned) {
            out << std::setw(16) << std::setfill(' ') << time 
                << std::setw(18) << std::setfill(' ') << current 
                << std::endl;
        }
        //
        time += timeStep;
    }
    out.close();

    Magnitude::operator=(
        Magnitude(
            new math::function::LinearInterpolation<math::Real,math::Real>(
                file)));
}

bool Numerical::operator==(const Numerical& rhs) const {
    bool areEqual = true;
    areEqual &= Magnitude::operator==(rhs);
    areEqual &= (bool) file.compare(rhs.file);
    return areEqual;
}

math::Real Numerical::evaluate(const math::Real time) const {
    throw std::logic_error("Numerical::evaluate not implemented");
}

}
}
} 
