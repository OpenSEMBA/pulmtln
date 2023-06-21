#include "Exporter.h"

#include "core/geometry/element/Line2.h"
#include "core/geometry/element/Quadrilateral4.h"
#include "core/geometry/element/Hexahedron8.h"

namespace semba {
namespace exporters {

Exporter::Exporter(const std::string& name)
:   ProjectFile(name) {}

void Exporter::deleteExistentOutputFiles() const {
    std::string file;
    file = getFilename() + ".post.msh";
    std::remove(file.c_str());
    file = getFilename() + ".post.res";
    std::remove(file.c_str());
}

std::string Exporter::getOutputfilename() const {
    std::string folder = getFolder();
    std::string output = getOutputName();
    return folder + output;
}

geometry::ElemR* Exporter::getBoundary(
        const math::Constants::CartesianAxis dir,
        const math::Constants::CartesianBound pos,
        geometry::CoordR3Group& cG,
        const geometry::Grid3* grid,
        const geometry::mesh::Mesh* mesh) const {
    geometry::BoxR3 box;
    if (grid != nullptr) {
        box = grid->getFullDomainBoundingBox();
    } else {
        box = mesh->getBoundingBox();
    }
    geometry::BoxR3 quadBox = box.getBoundAsBox(dir,pos);

    if (!quadBox.isSurface()) {
        throw std::runtime_error("Not surface");
    }

    auto posVec = quadBox.getPos();
    std::vector<const geometry::CoordR3*> coords;
    for (std::size_t i = 0; i < geometry::Qua4::sizeOfCoordinates; i++) {
        coords.push_back(
            cG.addAndAssignId(
                std::make_unique<geometry::CoordR3>(geometry::CoordId(), posVec[i])
            )->get()
        );
    }

    return new geometry::QuaR4(geometry::ElemId(0), coords.data());
}

std::string Exporter::getBoundaryName(
        const geometry::mesh::Structured* mesh,
        const std::size_t i,
        const std::size_t j) {
    auto boundType = mesh->bounds()(i, j);

    std::string boundName;
    if (boundType == nullptr) {
        boundName = "Undefined";
    } else {
        boundName = boundType->castTo<physicalModel::PhysicalModel>()->getName();
    }
    return boundName + "@Boundary";
}

ElemRGroup Exporter::getGridElems(
        geometry::CoordR3Group& cG,
        const geometry::Grid3* grid
) const {
    auto elem = ElemRGroup();

    if (grid == nullptr) {
        return elem;
    }

    geometry::BoxR3 box = grid->getFullDomainBoundingBox();

    for (std::size_t d = 0; d < 3; d++) {
        // Generates grid as lines.
        for (std::size_t i = 0; i < 2; i++) {
            std::vector<math::Real> pos = grid->getPos((d+i+1)%3);
            for (std::size_t j = 0; j < pos.size(); j++) {
                math::CVecR3 pMin, pMax;
                pMin(d) = grid->getPos(d,math::Constants::L);
                pMin((d+i+1)%3) = pos[j];
                pMax = pMin;
                pMin((d-i+2)%3) = grid->getPos((d-i+2)%3).front();
                pMax((d-i+2)%3) = grid->getPos((d-i+2)%3).back();

                auto newBox = geometry::BoxR3(pMin, pMax);
                if(!newBox.isLine()) {
                    throw std::runtime_error("Not line");
                }

                auto posVec = newBox.getPos();
                std::vector<const geometry::CoordR3*> coords;
                for (std::size_t i = 0; i < geometry::LinR2::sizeOfCoordinates; i++) {
                    coords.push_back(
                        cG.addAndAssignId(
                            std::make_unique<geometry::CoordR3>(geometry::CoordId(), posVec[i])
                        )->get()
                    );
                }

                elem.addAndAssignId(
                    std::make_unique<geometry::LinR2>(
                        geometry::ElemId(0),
                        coords.data()
                    )
                );
            }
        }
    }

    return elem;
}

}
} 
