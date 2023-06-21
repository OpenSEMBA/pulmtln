#pragma once

#include "Tetrahedron.h"

namespace semba {
namespace geometry {
namespace element {

class Tetrahedron4 : public Tetrahedron {
public:
    static const std::size_t sizeOfCoordinates = 4;

    Tetrahedron4(const Id id,
                 const CoordR3* v[4],
                 const Layer* lay = nullptr,
                 const Model* mat = nullptr);
    Tetrahedron4(const Tetrahedron4& rhs);
    virtual ~Tetrahedron4() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Tetrahedron4>(*this);
    }

    bool isInnerPoint(const math::CVecR3& pos) const override;
    bool isCurvedFace(const std::size_t face) const override;
    bool isFaceContainedInPlane(
            const std::size_t face,
            const math::Constants::CartesianPlane plane) const override;

    std::size_t numberOfCoordinates() const override { return sizeOfCoordinates; }

    std::size_t numberOfSideCoordinates(const std::size_t f = 0) const override {
        return 3;
    }

    const CoordR3* getV(const std::size_t i) const override { return v_[i]; }
    const CoordR3* getSideV(const std::size_t f, const std::size_t i) const override;

    const CoordR3* getVertex(const std::size_t i) const override {
        return v_[tet.vertex(i)];
    }
    const CoordR3* getSideVertex(const std::size_t f,
                                 const std::size_t i) const override;

    math::Real getVolume() const override;
    const math::simplex::Simplex& getTet() const override { return tet; }
    math::Real getAreaOfFace(const std::size_t face) const override;

    void setV(const std::size_t i, const CoordR3*) override;
    void check() const;

    virtual std::unique_ptr<Element<math::Int >> toStructured(
        const CoordI3Group&,
        const Grid3&,
        const math::Real = Grid3::tolerance) const override;

    virtual std::unique_ptr<Element<math::Real>> toUnstructured(
        const CoordR3Group&,
        const Grid3&) const override;

private:
    static const math::simplex::Tetrahedron<1> tet;

    // TODO: Remove plain array
    const CoordR3* v_[4];

    bool hasZeroVolume() const;
};

} 

typedef element::Tetrahedron4 Tet4;

} 
} 

