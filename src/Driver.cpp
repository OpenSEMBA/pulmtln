#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

namespace pulmtln {

using namespace mfem;



Driver Driver::loadFromFile(const std::string& fn)
{
    Parser p{ fn };
    return Driver{ 
        p.readModel(), 
        p.readDriverOptions() 
    };
}

Driver::Driver(const Model& model, const DriverOptions& opts) :
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

std::vector<int> getAttributesInMap(const NameToAttrMap& m)
{
    std::vector<int> res;
    for (const auto& [k,v] : m) {
        res.push_back(v);
    }
    return res;
}

mfem::DenseMatrix solveCMatrix(
    const Model& model, 
    const DriverOptions& opts,
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


    const bool isOpenProblem{model.isOpen()};
    // Solves a electrostatic problem for each conductor besides the
    // reference conductor (Conductor_0).
    const auto pecToBdrMap{ mats.buildNameToAttrMap<PEC>() };  
    for (const auto& [nameI, bdrAttI] : pecToBdrMap) {
        auto numI{ getNumberContainedInName(nameI) };
        if (numI == 0) {
            continue;
        }
        
        Mesh mesh{ *model.getMesh() };
        
        SolverParameters parameters;
        if (isOpenProblem) {
            parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMap<PEC>(), -1.0);
            parameters.dirichletBoundaries[bdrAttI] = 1.0;
        }
        else {
            parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMap<PEC>(), 0.0);
            parameters.dirichletBoundaries[bdrAttI] = 1.0;
        }
        parameters.domainPermittivities = domainToEpsr;
        parameters.openBoundaries = getAttributesInMap(mats.buildNameToAttrMap<OpenBoundary>());

        ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
        s.Solve();
        
        // Fills row
        for (const auto& [nameJ, bdrAttJ] : pecToBdrMap) {
            auto numJ{ getNumberContainedInName(nameJ) };
            if (numJ == 0) {
                continue;
            }
            auto charge{ s.chargeInBoundary(pecToBdrMap.at(nameJ)) };
            if (isOpenProblem) {
                C(numI - 1, numJ - 1) = charge / 2.0;
            }
            else {
                C(numI - 1, numJ - 1) = charge;
            }
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

        if (opts.exportVisItSolution) {
            std::string outputName{ opts.exportFolder + "/" + "VisIt/Conductor_" };
            outputName += std::to_string(numI);
            if (ignoreDielectrics) {
                outputName += "_no_dielectrics";
            }
            VisItDataCollection dC{ outputName, s.getMesh() };
            s.writeVisItFields(dC);
        }
    }
    return C;
}

mfem::DenseMatrix solveLMatrix(
    const Model& model, const DriverOptions& opts)
{
    // Inductance matrix can be computed from the 
    // capacitance obtained ignoring dielectrics as
    //          L = mu0 * eps0 * C^{-1}
    auto res{ solveCMatrix(model, opts, true) };
    res.Invert();
    res *= MU0_NATURAL * EPSILON0_NATURAL;
    return res;
}


Parameters Driver::getMTLPUL() const
{
    Parameters res;

    // Computes matrices in natural units.
    res.C = solveCMatrix(model_, opts_);
    res.L = solveLMatrix(model_, opts_);
    
    // Converts to SI units.
    res.C *= EPSILON0_SI;
    res.L *= MU0_SI;

    if (opts_.exportMatrices) {
        res.saveToJSONFile(opts_.exportFolder + "/matrices.pulmtln.out.json");
    }

    return res;
}




}