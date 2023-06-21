#pragma once

#include "core/math/Types.h"
#include "core/geometry/Box.h"
#include "core/util/Class.h"
#include "core/util/Identifiable.h"
#include "core/util/Identification.h"
#include "core/geometry/element/Element.h"
#include "core/geometry/element/Group.h"
#include "core/physicalModel/Group.h"

namespace semba {
namespace geometry {
namespace mesh {

class Mesh;
    typedef util::Identification<Mesh> Id;

class Mesh : public virtual util::Identifiable<Id>,
             public virtual util::Class {
public:
    Mesh() = default;
    virtual ~Mesh() = default;

    virtual void applyScalingFactor(const math::Real factor) = 0;
    virtual BoxR3 getBoundingBox() const = 0;
    virtual void reassignPointers(const PMGroup& matGr = PMGroup()) = 0;

    virtual std::unique_ptr<Mesh> clone() const = 0;

    virtual ElemView reassign(const ElemView&) = 0;
};

} 
} 
} 

