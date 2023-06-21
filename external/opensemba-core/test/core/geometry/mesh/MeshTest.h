#pragma once

#include <gtest/gtest.h>

#include "core/geometry/element/Tetrahedron4.h"
#include "core/geometry/mesh/Mesh.h"
#include "core/geometry/element/Group.h"

using namespace semba;
using namespace geometry;
using namespace math;

class MeshTest {
public:
    void SetUp() {
        CoordId coordId(1);
        cG_.add(std::make_unique<CoordR3>(coordId++, CVecR3(0.0, 0.0, 0.0)));
        cG_.add(std::make_unique<CoordR3>(coordId++, CVecR3(0.0, 0.0, 1.0)));
        cG_.add(std::make_unique<CoordR3>(coordId++, CVecR3(0.0, 1.0, 0.0)));
        cG_.add(std::make_unique<CoordR3>(coordId++, CVecR3(1.0, 0.0, 0.0)));
        
        const CoordR3* vTet[4] = {
                cG_.getId(CoordId(1)),
                cG_.getId(CoordId(2)),
                cG_.getId(CoordId(3)),
                cG_.getId(CoordId(4))
        };
        const CoordR3* vTri[3] = {
                cG_.getId(CoordId(2)),
                cG_.getId(CoordId(1)),
                cG_.getId(CoordId(3))
        };
        eG_ = ElemRGroup();
        eG_.add(std::make_unique<Tet4>(ElemId(1), vTet));
        eG_.add(std::make_unique<Tri3>(ElemId(2), vTri));
        lG_ = LayerGroup();
    }

protected:
    CoordR3Group cG_;
    ElemRGroup eG_;
    LayerGroup lG_;

};

