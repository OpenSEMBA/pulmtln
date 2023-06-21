#pragma once

#include "SolverOptions.h"
#include "constants.h"
#include "FES.h"
#include "Model.h"

namespace pulmtln {

using namespace mfem;

class ElectrostaticSolver {
public:
    ElectrostaticSolver(
        Model&, 
        const SolverOptions&);
    ~ElectrostaticSolver();

    void Assemble();
    void Solve();

    void RegisterParaViewFields(ParaViewDataCollection&);
    void WriteParaViewFields(ParaViewDataCollection&);

    const GridFunction& GetVectorPotential() { return *phi_; }
    double computeTotalChargeFromRho() const;
    double computeTotalChargeFromP() const;
private:
    SolverOptions opts_;
    
    Mesh* mesh_;

    Array<int>* dbcs_; // Dirichlet BC Surface Attribute IDs
    Vector* dbcv_;     // Corresponding Dirichlet Values

    H1_FESpace* H1FESpace_;    // Continuous space for phi
    ND_FESpace* HCurlFESpace_; // Tangentially continuous space for E
    RT_FESpace* HDivFESpace_;  // Normally continuous space for D
    L2_FESpace* L2FESpace_;    // Discontinuous space for rho

    BilinearForm* divEpsGrad_; // Laplacian operator
    BilinearForm* h1Mass_;     // For Volumetric Charge Density Source
    BilinearForm* h1SurfMass_; // For Surface Charge Density Source
    BilinearForm* hDivMass_;   // For Computing D from E

    MixedBilinearForm* hCurlHDivEps_; // For computing D from E
    MixedBilinearForm* hCurlHDiv_;    // For computing D from E and P
    MixedBilinearForm* weakDiv_;      // For computing the source term from P

    LinearForm* rhod_; // Dual of Volumetric Charge Density Source

    LinearForm* l2_vol_int_;  // Integral of L2 field
    LinearForm* rt_surf_int_; // Integral of H(Div) field over boundary

    DiscreteGradOperator* grad_; // For Computing E from phi
    DiscreteDivOperator* div_;  // For Computing rho from D

    GridFunction* phi_;       // Electric Scalar Potential
    GridFunction* rho_;       // Volumetric Charge Density (Div(D))
    GridFunction* e_;         // Electric Field
    GridFunction* d_;         // Electric Flux Density (aka Dielectric Flux)

    ConstantCoefficient oneCoef_;   // Coefficient equal to 1
    ConstantCoefficient epsCoef_;   // Dielectric Permittivity Coefficient

    Array<int> ess_bdr_, ess_bdr_tdofs_; // Essential Boundary Condition DoFs
};

}