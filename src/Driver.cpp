#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

namespace pulmtln {

using namespace mfem;

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
    auto pecToBdrMap{ model.getMaterialsOfType(MaterialType::PEC) };
    if (pecToBdrMap.size() < 2) {
        throw std::runtime_error(
            "The number of conductors must be greater than 2."
        );
    }
    
    mfem::DenseMatrix C((int) pecToBdrMap.size() - 1);

    std::map<int, double> domainToEpsr{};
    if (!ignoreDielectrics) {
        // TODO
    }

    for (const auto& [name, bdrAtt] : pecToBdrMap) {
        auto num{ getNumberContainedInName(name) };
        if (num == 0) {
            continue;
        }
        
        BdrConditionValues bcs{ 
            initializeBdrConditionValuesTo(pecToBdrMap, 0.0) 
        };
        bcs[bdrAtt] = 1.0;
            
        Mesh mesh{ *model.getMesh() };
        ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
        s.Solve();
       
        C(num - 1, num - 1) = s.chargeInBoundary(pecToBdrMap.at(name));
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