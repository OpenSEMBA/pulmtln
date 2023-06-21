#pragma once

#include "core/math/simplex/Tetrahedron.h"
#include "Triangle3.h"
#include "Volume.h"

namespace semba {
namespace geometry {
namespace element {

class Tetrahedron : public Volume<math::Real> {
public:
    Tetrahedron();
    virtual ~Tetrahedron();

    virtual bool isCurvedFace(const std::size_t face) const = 0;
    virtual bool isFaceContainedInPlane(const std::size_t face,
            const math::Constants::CartesianPlane plane) const = 0;

    inline std::size_t numberOfFaces   () const { return 4; }
    inline std::size_t numberOfVertices() const { return 4; }
    inline std::size_t numberOfSideVertices(const std::size_t f = 0) const {
        return 3;
    }
    virtual const math::simplex::Simplex& getTet() const = 0;
    virtual math::Real getVolume() const = 0;
    virtual math::Real getAreaOfFace(const std::size_t face) const = 0;
    virtual Triangle3* getTri3Face(const std::size_t f) const;

};

} 

typedef element::Tetrahedron Tet;

} 
} 

