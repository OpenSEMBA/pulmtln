#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

namespace pulmtln {

using namespace mfem;

mfem::DenseMatrix MTLPULParameters::getCapacitiveCouplingCoefficients() const
{
    mfem::DenseMatrix r(C.NumRows(), C.NumCols());
    for (auto i{ 0 }; i < C.NumRows(); ++i) {
        auto selfC{ C(i,i) };
        for (auto j{ 0 }; j < C.NumCols(); ++j) {
            r(i, j) = C(i, j) / selfC;
        }
    }
    return r;
}

Driver Driver::loadFromFile(const std::string& fn)
{
    Parser p{ fn };
    return Driver{ p.readModel(), p.readSolverOptions() };
}

Driver::Driver(const Model& model, const SolverOptions& opts) :
    model_{model},
    opts_{opts}
{}

int getNumberContainedInName(const std::string& name) 
{
    std::stringstream ss{ name.substr(name.find("_") + 1) };
    int res;
    ss >> res;
    return res;
}

BdrConditionValues initializeBdrConditionValuesTo(
    const MatNameToAttribute& matToAtt, double value)
{
    BdrConditionValues bcs;
    for (const auto& [mat, att] : matToAtt) {
        bcs[att] = value;
    }
    return bcs;
}

mfem::DenseMatrix solveCMatrix(
    const Model& model, 
    const SolverOptions& opts,
    bool ignoreDielectrics = false)
{
    const auto& mats{ model.getMaterials() };

    auto pecToBdrMap{ mats.getMatNameToAttributeMap<PEC>() };
    if (pecToBdrMap.size() < 2) {
        throw std::runtime_error(
            "The number of conductors must be greater than 2."
        );
    }
    
    int CSize{ (int)pecToBdrMap.size() - 1 };
    mfem::DenseMatrix C(CSize);

    std::map<int, double> domainToEpsr;

    for (const auto& d : mats.dielectrics) {
        if (ignoreDielectrics) {
            domainToEpsr[d.tag] = 1.0;
        } else {
            domainToEpsr[d.tag] = d.relativePermittivity;
        }
    }

    for (const auto& [nameI, bdrAttI] : pecToBdrMap) {
        auto numI{ getNumberContainedInName(nameI) };
        if (numI == 0) {
            continue;
        }
        
        BdrConditionValues bcs{ 
            initializeBdrConditionValuesTo(pecToBdrMap, 0.0) 
        };
        bcs[bdrAttI] = 1.0;
            
        Mesh mesh{ *model.getMesh() };
        ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
        s.Solve();
        
        for (const auto& [nameJ, bdrAttJ] : pecToBdrMap) {
            auto numJ{ getNumberContainedInName(nameJ) };
            if (numJ == 0) {
                continue;
            }
            C(numI - 1, numJ - 1) = s.chargeInBoundary(pecToBdrMap.at(nameJ));
        }

        if (opts.exportParaViewSolution) {
            std::string outputName{"ParaView/DriverResult"};
            if (ignoreDielectrics) {
                outputName += "_no_dielectrics_";
            }
            outputName += std::to_string(numI);
            ParaViewDataCollection pd{ outputName, s.getMesh() };
            s.writeParaViewFields(pd);
        }
    }
    return C;
}

mfem::DenseMatrix solveLMatrix(
    const Model& model, const SolverOptions& opts)
{
    // Inductance matrix can be computed from the 
    // capacitance obtained ignoring dielectrics as
    //          L = mu0 * eps0 * C^{-1}
    auto res{ solveCMatrix(model, opts, true) };
    res.Invert();
    res *= MU0 * EPSILON0;
    return res;
}


MTLPULParameters Driver::getMTLPUL() const
{
    MTLPULParameters res;

    // Computes matrices in natural units.
    res.C = solveCMatrix(model_, opts_);
    res.L = solveLMatrix(model_, opts_);
    
    // Converts to SI units.
    res.C *= EPSILON0_SI;
    res.L *= MU0_SI;

    return res;
}




}