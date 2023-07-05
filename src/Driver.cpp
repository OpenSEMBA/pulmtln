#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

namespace pulmtln {

using namespace mfem;

Driver loadFromFile(const std::string& fn)
{
    Parser p{ fn };
    return Driver{ p.readModel(), p.readSolverOptions() };
}

Driver::Driver(const Model& model, const SolverOptions& opts) :
    model_{model},
    opts_{opts}
{}

int getNameNumber(const std::string& name) 
{
    std::stringstream ss{ name.substr(name.find("_") + 1) };
    int res;
    ss >> res;
    return res;
}

BdrConditionValues initializeBdrsToValue(
    const MatNameToAttribute& matToAtt, double value)
{
    BdrConditionValues bcs;
    for (const auto& [mat, att] : matToAtt) {
        bcs[att] = value;
    }
    return bcs;
}

mfem::DenseMatrix solveCMatrix(const Model& model, const SolverOptions& opts)
{
    auto pecToBdrMap{ model.getMaterialsOfType(MaterialType::PEC) };
    if (pecToBdrMap.size() < 2) {
        throw std::runtime_error(
            "The number of conductors must be greater than 2."
        );
    }
    
    mfem::DenseMatrix C(pecToBdrMap.size() - 1);

    BdrConditionValues bcs{ initializeBdrsToValue(pecToBdrMap, 0.0) };

    for (const auto& [name, bdrAtt] : pecToBdrMap) {
        auto num{ getNameNumber(name) };
        if (num == 0) {
            continue;
        }
        
        bcs[bdrAtt] = 1.0;
    }
}

MTLPULParameters Driver::getMTLPUL() const
{
    MTLPULParameters res;

    res.C = solveCMatrix(model_, opts_);

    return res;
}




}