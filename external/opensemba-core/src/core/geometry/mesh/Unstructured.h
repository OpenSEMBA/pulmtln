#pragma once

#include "Mesh.h"

#include "core/geometry/coordinate/Group.h"
#include "core/geometry/element/Group.h"
#include "core/geometry/Layer.h"

namespace semba {
namespace geometry {
namespace mesh {

class Structured;

class Unstructured : public Mesh {
public:
    Unstructured() = default;
    Unstructured(const CoordR3Group& cG,
                 const ElemRGroup& elem,
                 const LayerGroup& = LayerGroup());
    Unstructured(const Unstructured& param);
    Unstructured(Unstructured&&) = default;
    virtual ~Unstructured() = default;

    Unstructured& operator=(const Unstructured&);
    Unstructured& operator=(Unstructured&&) = default;

    virtual std::unique_ptr<Mesh> clone() const override {
        return std::make_unique<Unstructured>(*this);
    }

    CoordR3Group& coords() { return coords_; }
    ElemRGroup&      elems () { return elems_; }
    LayerGroup& layers() { return layers_; }

    const CoordR3Group& coords() const { return coords_; }
    const ElemRGroup&      elems () const { return elems_; }
    const LayerGroup& layers() const { return layers_; }

    Structured* getMeshStructured(
            const Grid3& grid,
            const math::Real tol = Grid3::tolerance) const;

    void applyScalingFactor(const math::Real factor) override;
    BoxR3 getBoundingBox() const override;

    void reassignPointers(const PMGroup& matGr = PMGroup()) override;

    virtual ElemView reassign(const ElemView&) override;

private:
	CoordR3Group coords_;
    ElemRGroup elems_;
	LayerGroup layers_;
};

} 

using UnstructuredMesh = mesh::Unstructured;

} 
} 

