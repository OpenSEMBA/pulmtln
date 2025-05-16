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
    InCellParameters getInCellParameters() const;

    FloatingPotentials getFloatingPotentials() const;

    static Driver loadFromFile(const std::string& filename);
    static SolverParameters buildSolverParametersFromModel(
        const Model& model,
        bool ignoreDielectrics);

    static DenseMatrix getCMatrix(
        const Model& model, const DriverOptions& opts,
        bool ignoreDielectrics = false, bool generalized = false);
    static DenseMatrix getLMatrix(
        const Model& model, const DriverOptions& opts);
    static DenseMatrix getFloatingPotentialsMatrix(
        const Model& model, const DriverOptions& opts,
        const bool ignoreDielectrics);

    void setExportFolder(const std::string folder) { opts_.exportFolder = folder; }
private:
    Model model_;
    DriverOptions opts_;
};

}