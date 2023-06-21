#pragma once

#include <fstream>
#include <utility>
#include <algorithm>

#include "Exporter.h"
#include "core/ProblemDescription.h"

namespace semba {
namespace exporters {

class ExporterVTK : public exporters::Exporter {
public:
    ExporterVTK(const UnstructuredProblemDescription&, const std::string&);
    
private:
    void writeMesh_(const UnstructuredProblemDescription&);
    void writeFile_(const exporters::ElemRView& elems,
                    const std::string& name,
                    std::ofstream& outMain,
                    std::size_t& part);
    std::pair<std::vector<math::CVecR3>, 
              std::map<geometry::CoordId, std::size_t>> getPoints_(
              const exporters::ElemRView& elems);
    void writePoints_(std::ofstream& outFile,
                      const std::vector<math::CVecR3>& pos);
    void writeCells_(
            std::ofstream& outFile,
            const exporters::ElemRView& elems,
            std::map<geometry::CoordId, std::size_t>& mapCoords);

    static std::string makeValid_(const std::string&);
};

}
} 