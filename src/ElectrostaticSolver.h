#pragma once

#include "SolverOptions.h"
#include "constants.h"
#include "FES.h"
#include "AttrToValueMap.h"

namespace pulmtln {

using namespace mfem;

struct SolverParameters {
    AttrToValueMap dirichletBoundaries; 
    AttrToValueMap neumannBoundaries;
    std::vector<int> openBoundaries;
    AttrToValueMap domainPermittivities;
};

class ElectrostaticSolver {
public:
    ElectrostaticSolver(
        Mesh& mesh,
        const SolverParameters&,
        const SolverOptions = SolverOptions{});
    ~ElectrostaticSolver();

    void Solve();

    void writeParaViewFields(ParaViewDataCollection&) const;
    void writeVisItFields(VisItDataCollection&) const;

    const GridFunction& GetPotential() { return *phi_; }
    const GridFunction& GetElectricField() { return *e_; }

    double totalChargeFromRho() const;
    double totalCharge() const;
    double chargeInBoundary(int bdrAttribute) const;
    double totalEnergy() const;

    Mesh* getMesh() { return mesh_; }

private:
    SolverOptions opts_;
    
    Mesh* mesh_;

    SolverParameters parameters_;

    H1_FESpace* H1FESpace_;    // Continuous space for phi
    ND_FESpace* HCurlFESpace_; // Tangentially continuous space for E
    RT_FESpace* HDivFESpace_;  // Normally continuous space for D
    L2_FESpace* L2FESpace_;    // Discontinuous space for rho

    BilinearForm* divEpsGrad_; // Laplacian operator
    BilinearForm* hDivMass_;   // For Computing D from E
    BilinearForm* h1SurfMass_; // For Surface Charge Density Source 

    MixedBilinearForm* hCurlHDivEps_; // For computing D from E
    
    LinearForm* rhod_; // Dual of Volumetric Charge Density Source

    DiscreteGradOperator* grad_; // For Computing E from phi
    DiscreteDivOperator* div_;  // For Computing rho from D

    GridFunction* phi_;       // Electric Scalar Potential
    GridFunction* rho_;       // Volumetric Charge Density (Div(D))
    GridFunction* e_;         // Electric Field
    GridFunction* d_;         // Electric Flux Density (aka Dielectric Flux)
    GridFunction* sigma_src_; // Surface Charge Density Source

    ConstantCoefficient oneCoef_;   // Coefficient equal to 1
    Coefficient* epsCoef_;   // Dielectric Permittivity Coefficient

    Array<int> ess_bdr_, ess_bdr_tdofs_, open_bdr_; // Essential Boundary Condition DoFs

    void Assemble();
    void applyBoundaryValuesToGridFunction(
        const AttrToValueMap& bdrValues, 
        GridFunction& gf) const;
};

}