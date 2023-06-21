#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/ProblemDescription.h"
#include "core/geometry/mesh/Structured.h"

namespace semba {
namespace exporters {

using ElemRGroup = geometry::element::Group<geometry::ElemR>;
using ElemRView = std::vector<const geometry::ElemR*>;

class Exporter : public util::ProjectFile {

public:
    Exporter(const std::string& name);
    virtual ~Exporter() = default;
    
protected:
    void deleteExistentOutputFiles() const;
    std::size_t determineStepsSaved(const math::Real savingPeriod, const math::Real dt) const;
    
    std::string getOutputfilename() const;
	geometry::ElemR* getBoundary(
            const math::Constants::CartesianAxis dir,
            const math::Constants::CartesianBound bound,
            geometry::CoordR3Group& cG,
            const geometry::Grid3* grid,
            const geometry::mesh::Mesh* mesh) const;
    ElemRGroup getGridElems(
            geometry::CoordR3Group& cG,
            const geometry::Grid3* grid) const;
    static std::string getBoundaryName(
            const geometry::mesh::Structured* mesh,
            const std::size_t i,
            const std::size_t j);
};

} 
}