#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Parameters.h"

namespace pulmtln {

class Driver {
public:
    Driver(Model&& model, const DriverOptions& opts);
    
    PULParameters getMTLPUL() const;
    PULParametersByDomain getMTLPULByDomains() const;
    mfem::DenseMatrix getFloatingPotentialsMatrix() const;

    static Driver loadFromFile(const std::string& filename);

    void setExportFolder(const std::string folder) { opts_.exportFolder = folder; }
private:
    Model model_;
    DriverOptions opts_;
};

}