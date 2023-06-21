#pragma once

#include "core/math/function/Function.h"
#include "core/util/Class.h"

namespace semba {
namespace source {
namespace Magnitude {

class Magnitude : public virtual util::Class {
public:
    Magnitude() = default;
    Magnitude(math::FunctionRR* mathFunction);
    Magnitude(const Magnitude& rhs);
    virtual ~Magnitude() = default;

    Magnitude& operator=(const Magnitude& rhs);

    virtual std::unique_ptr<Magnitude> clone() const {
        return std::make_unique<Magnitude>(*this);
    }

    virtual bool operator==(const Magnitude&) const;

    math::Real evaluate(const math::Real time) const;

private:
    math::FunctionRR* mathFunction_;
};

}
}
} 

