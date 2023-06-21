#pragma once

#include "core/geometry/element/Group.h"
#include "core/util/ProjectFile.h"
#include "core/util/Class.h"
#include "core/util/Identifiable.h"
#include "core/util/Identification.h"
#include "magnitude/Magnitude.h"
#include "magnitude/Numerical.h"

namespace semba {
namespace source {

class Source;
    typedef util::Identification<Source> Id;

class Source : public virtual util::Identifiable<Id>,
               public virtual util::Class {
public:
    using Target = std::vector<geometry::ElemId>;

    Source(const std::unique_ptr<Magnitude::Magnitude>&, const Target&);
    Source(const Source&);
    virtual ~Source() = default;

    Source& operator=(const Source&);
    
    virtual std::unique_ptr<Source> clone() const = 0;

    virtual std::string getName() const = 0;

    void convertToNumerical(const util::ProjectFile& file,
                            const math::Real step,
                            const math::Real finalTime);
    
    Magnitude::Numerical exportToFile(const util::ProjectFile& file,
                                       const math::Real step,
                                       const math::Real finalTime) const;
    const Magnitude::Magnitude* getMagnitude() const { return magnitude_.get(); }

    Target getTarget() const { return target_; }
    void setTarget(const Target& target) { target_ = target; }

protected:
    std::unique_ptr<Magnitude::Magnitude> magnitude_;
    Target target_;
};

}
} 


