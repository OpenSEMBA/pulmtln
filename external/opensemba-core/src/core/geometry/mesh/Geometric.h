#pragma once

#include "core/geometry/Grid.h"

#include "Unstructured.h"

namespace semba {
namespace geometry {
namespace mesh {

class Geometric : public Unstructured {
public:
    Geometric() = default;
    Geometric(const Grid3& grid);
    Geometric(const Grid3& grid, 
              const CoordR3Group& cG,
              const ElemRGroup& elem,
              const LayerGroup& = LayerGroup());
    Geometric(const Geometric&);
    Geometric(const Grid3& grid, const Unstructured& unstructured);
    Geometric(Geometric&&) = default;
    virtual ~Geometric() = default;

    Geometric& operator=(const Geometric& rhs);
    Geometric& operator=(Geometric&&) = default;

    virtual std::unique_ptr<Mesh> clone() const override {
        return std::make_unique<Geometric>(*this);
    }

    Grid3&       grid()       { return grid_; }
    const Grid3& grid() const { return grid_; }

    Structured* getMeshStructured(const math::Real = Grid3::tolerance) const;

    void applyScalingFactor(const math::Real factor) override;

private:
	Grid3 grid_;
};

} 
} 
} 

