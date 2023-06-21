#pragma once

#include "SolverOptions.h"
#include "constants.h"
#include "FES.h"

namespace mfem {
namespace pulmtln {

class Solver {
public:
    Solver(Mesh& pmesh, int order,
        Array<int>& dbcs, Vector& dbcv,
        Array<int>& nbcs, Vector& nbcv,
        Coefficient& epsCoef,
        double (*phi_bc)(const Vector&),
        double (*rho_src)(const Vector&),
        void   (*p_src)(const Vector&, Vector&),
        Vector& point_charges);
    ~Solver();

    void PrintSizes();
    void Assemble();
    void Update();
    void Solve();
    
    const GridFunction& GetVectorPotential() { return *phi_; }

private:
    SolverOptions opts_;
    
    Mesh* mesh_;

    Array<int>* dbcs_; // Dirichlet BC Surface Attribute IDs
    Vector* dbcv_;     // Corresponding Dirichlet Values
    Array<int>* nbcs_; // Neumann BC Surface Attribute IDs
    Vector* nbcv_;     // Corresponding Neumann Values

    ParaViewDataCollection* paraview_dc_; // To prepare fields for Paraview viewing

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
    GridFunction* rho_src_;   // Volumetric Charge Density Source
    GridFunction* rho_;       // Volumetric Charge Density (Div(D))
    GridFunction* sigma_src_; // Surface Charge Density Source
    GridFunction* e_;         // Electric Field
    GridFunction* d_;         // Electric Flux Density (aka Dielectric Flux)

    ConstantCoefficient oneCoef_;   // Coefficient equal to 1
    Coefficient* epsCoef_;   // Dielectric Permittivity Coefficient
    Coefficient* phiBCCoef_; // Scalar Potential Boundary Condition
    Coefficient* rhoCoef_;   // Charge Density Coefficient

    // Source functions
    double (*phi_bc_func_)(const Vector&);          // Scalar Potential BC
    double (*rho_src_func_)(const Vector&);          // Volumetric Charge Density

    Array<int> ess_bdr_, ess_bdr_tdofs_; // Essential Boundary Condition DoFs
};

}
} // namespace mfem