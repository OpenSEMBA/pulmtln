#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

using namespace mfem;

namespace pulmtln {

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

void exportFieldSolutions(
    const DriverOptions& opts, 
    ElectrostaticSolver& s,
    const std::string name, 
    bool ignoreDielectrics)
{
    if (opts.exportParaViewSolution) {
        std::string outputName{ opts.exportFolder + "/" + "ParaView/" + name };
        if (ignoreDielectrics) {
            outputName += "_no_dielectrics";
        }
        ParaViewDataCollection pd{ outputName, s.getMesh() };
        s.writeParaViewFields(pd);
    }

    if (opts.exportVisItSolution) {
        std::string outputName{ opts.exportFolder + "/" + "VisIt/" + name };
        if (ignoreDielectrics) {
            outputName += "_no_dielectrics";
        }
        VisItDataCollection dC{ outputName, s.getMesh() };
        s.writeVisItFields(dC);
    }
}


SolverParameters buildSolverParameters(
    const Model& model, 
    bool ignoreDielectrics)
{
    const Materials& mats = model.getMaterials();

    SolverParameters parameters;

    auto domainToEpsr{
        buildAttrToValueMap(mats.buildNameToAttrMap<Dielectric>(), 1.0)
    };
    if (!ignoreDielectrics) {
        for (const auto& d : mats.dielectrics) {
            domainToEpsr.at(d.attribute) = d.relativePermittivity;
        }
    }
    parameters.domainPermittivities = domainToEpsr;
    parameters.openBoundaries = getAttributesInMap(mats.buildNameToAttrMap<OpenBoundary>());

    if (model.determineOpenness() == Model::OpennessType::open) {
        parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMap<PEC>(), -1.0);
    }
    else {
        parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMap<PEC>(), 0.0);
    }
 
    return parameters;
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
    
    mfem::DenseMatrix C((int)mats.pecs.size() - 1);
    // Solves a electrostatic problem for each conductor besides the
    // reference conductor (Conductor_0).
    const auto pecToBdrMap{ mats.buildNameToAttrMap<PEC>() };  
    for (const auto& [nameI, bdrAttI] : pecToBdrMap) {
        auto numI{ Materials::getNumberContainedInName(nameI) };
        if (numI == 0) {
            continue;
        }
        
        auto parameters{ buildSolverParameters(model, ignoreDielectrics) };
        parameters.dirichletBoundaries[bdrAttI] = 1.0;

        Mesh mesh{ *model.getMesh() };
        
        ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
        s.Solve();
        
        // Fills row
        for (const auto& [nameJ, bdrAttJ] : pecToBdrMap) {
            auto numJ{ Materials::getNumberContainedInName(nameJ) };
            if (numJ == 0) {
                continue;
            }
            auto charge{ s.chargeInBoundary(pecToBdrMap.at(nameJ)) };
            if (model.determineOpenness() == Model::OpennessType::open) {
                C(numI - 1, numJ - 1) = charge / 2.0;
            }
            else {
                C(numI - 1, numJ - 1) = charge;
            }
        }

        exportFieldSolutions(opts, s, nameI, ignoreDielectrics);
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

PULParameters buildPULParametersForModel(const Model& m, const DriverOptions& opts)
{
    PULParameters res;

    res.C = solveCMatrix(m, opts);
    res.C *= EPSILON0_SI;

    res.L = solveLMatrix(m, opts);
    res.L *= MU0_SI;
    
    return res;
}

PULParameters Driver::getMTLPUL() const
{
    auto res{ buildPULParametersForModel(model_, opts_)};

    if (opts_.exportMatrices) {
        res.saveToJSONFile(opts_.exportFolder + "/matrices.pulmtln.out.json");
    }

    return res;
}

Model buildModelForDomain(Model& model, const Domain& domain)
{
    auto& globalMesh{ *model.getMesh() };
    
    // We must modify attributes in the global mesh to identify elements in subdomain.
    auto globalMeshBackup{ globalMesh };
    
    // --
    for (auto e{ 0 }; e < globalMesh.GetNE(); ++e) {
        if (domain.elems.count(e)) {
            globalMesh.GetElement(e)->SetAttribute(1);
        }
        else {
            globalMesh.GetElement(e)->SetAttribute(0);
        }
    }
    globalMesh.SetAttributes();

    // --
    Array<int> subdomainAttrs(1);
    subdomainAttrs[0] = 1;
    auto domainMesh{ SubMesh::CreateFromDomain(globalMesh, subdomainAttrs)};

    // Restores original attributes in global mesh.
    for (auto e{ 0 }; e < globalMesh.GetNE(); ++e) {
        globalMesh.GetElement(e)->SetAttribute(
            globalMeshBackup.GetElement(e)->GetAttribute()
        );
    }
    globalMesh.SetAttributes();


    return Model{domainMesh, model.getMaterials()};
}

PULParametersByDomain Driver::getMTLPULByDomain() const
{
    PULParametersByDomain res;

    auto idToDomain{ Domain::buildDomains(model_) };
    
    for (const auto& [id, domain] : idToDomain) {
        res.domainToPUL[id] = buildPULParametersForModel(
            buildModelForDomain(model_, domain),
            opts_
        );
    }
    
    res.domainTree = DomainTree{idToDomain};

    return res;
}


}