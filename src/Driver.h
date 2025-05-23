#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Results.h"
#include "ElectrostaticSolver.h"

namespace pulmtln {

class Driver {
public:
    Driver(Model&& model, const DriverOptions& opts);
    
    PULParameters getMTLPUL() const;
    PULParametersByDomain getMTLPULByDomains() const;
    InCellPotentials getInCellPotentials() const;
    FloatingPotentials getFloatingPotentials() const;

    void setExportFolder(const std::string folder) { opts_.exportFolder = folder; }

    static SolverInputs buildSolverInputsFromModel(
        const Model& model,
        bool ignoreDielectrics);

    static Driver loadFromFile(const std::string& filename);
    static DenseMatrix getCMatrix(
        const Model& model, const DriverOptions& opts,
        bool ignoreDielectrics = false, bool generalized = false);
    static DenseMatrix getLMatrix(
        const Model& model, const DriverOptions& opts);
    static DenseMatrix getFloatingPotentialsMatrix(
        const Model& model, const DriverOptions& opts,
        const bool ignoreDielectrics);

    static double getInCellCapacitanceUsingInnerRegion(
        const InCellPotentials& potential, int i, int j);
    static double getInCellInductanceUsingInnerRegion(
        const InCellPotentials& potential, int i, int j);

private:
    Model model_;
    DriverOptions opts_;
};

}