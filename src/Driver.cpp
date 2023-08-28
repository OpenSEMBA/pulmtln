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
    return Driver{ 
        p.readModel(), 
        p.readSolverOptions() 
    };
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

AttrToValueMap buildAttrToValueMap(
    const NameToAttrMap& matToAtt, double value)
{
    AttrToValueMap bcs;
    for (const auto& [mat, att] : matToAtt) {
        bcs[att] = value;
    }
    return bcs;
}

AttrToValueMap buildBdrVoltagesWithZero(const Materials& mats)
{
    auto pecVoltages{ 
        buildAttrToValueMap(mats.buildNameToAttrMap<PEC>(), 0.0)};
    auto openRegionVoltages{
        buildAttrToValueMap(mats.buildNameToAttrMap<OpenBoundary>(), 0.0)
    };
    auto bdrVoltages{ pecVoltages };
    bdrVoltages.insert(openRegionVoltages.begin(), openRegionVoltages.end());
    return bdrVoltages;
}

mfem::DenseMatrix solveCMatrix(
    const Model& model, 
    const SolverOptions& opts,
    bool ignoreDielectrics = false)
{
    const auto& mats{ model.getMaterials() };
    if (mats.pecs.size() < 2) {
        throw std::runtime_error(
            "The number of conductors must be greater than 2."
        );
    }
    
    auto domainToEpsr{
        buildAttrToValueMap(mats.buildNameToAttrMap<Dielectric>(), 1.0)
    };
    if (!ignoreDielectrics) {
        for (const auto& d : mats.dielectrics) {
            domainToEpsr.at(d.tag) = d.relativePermittivity;
        }
    }

    mfem::DenseMatrix C((int)mats.pecs.size() - 1);

    // Solves a electrostatic problem for each conductor besides the
    // reference conductor (Conductor_0).
    const auto pecToBdrMap{ mats.buildNameToAttrMap<PEC>() };  
    for (const auto& [nameI, bdrAttI] : pecToBdrMap) {
        auto numI{ getNumberContainedInName(nameI) };
        if (numI == 0) {
            continue;
        }
        
        Mesh mesh{ *model.getMesh() };
        auto bdrVoltages{ buildBdrVoltagesWithZero(mats) };
        bdrVoltages[bdrAttI] = 1.0;

        ElectrostaticSolver s(mesh, bdrVoltages, domainToEpsr, opts);
        s.Solve();
        
        // Fills row
        for (const auto& [nameJ, bdrAttJ] : pecToBdrMap) {
            auto numJ{ getNumberContainedInName(nameJ) };
            if (numJ == 0) {
                continue;
            }
            C(numI - 1, numJ - 1) = s.chargeInBoundary(pecToBdrMap.at(nameJ));
        }

        if (opts.exportParaViewSolution) {
            std::string outputName{ opts.exportFolder + "/" + "ParaView/Conductor_"};
            outputName += std::to_string(numI);
            if (ignoreDielectrics) {
                outputName += "_no_dielectrics";
            }
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