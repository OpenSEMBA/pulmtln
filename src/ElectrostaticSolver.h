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

using multipolarCoefficient = std::pair<double, double>;
using multipolarCoefficients = std::vector<multipolarCoefficient>;

static double multipolarExpansion(
    const Vector& position, 
    const multipolarCoefficients& ab,
    const mfem::Vector expansionCenter)
{
    // 2D multiplar exapnsion from:
    // TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE. 
    // Lecture notes. 
    
    mfem::Vector rVec{ position - expansionCenter };
    double r{ rVec.Norml2() };
    double phi{ std::atan(rVec(1) / rVec(0)) };

    double res{ 0.0 };
    for (int n{ 0 }; n < ab.size(); ++n) {
        const auto& an = ab[n].first;
        const auto& bn = ab[n].second;
        if (n == 0) {
            res -= an * std::log(r);
            assert(bn == 0.0); // b0 should always be zero.
        }
        else {
            res +=  (an * std::cos(n * phi) + bn * std::sin(n * phi)) / std::pow(r, n);
        }
    }

    return res;
}

class ElectrostaticSolver {
public:
    ElectrostaticSolver(
        Mesh& mesh,
        const SolverParameters&,
        const SolverOptions = SolverOptions{});
    ~ElectrostaticSolver();

    void Solve();
    void setDirichletConditions(const AttrToValueMap& dbcs);
    void setNeumannCondition(const int bdrAttribute, Coefficient& chargeDensity);

    void writeParaViewFields(ParaViewDataCollection&) const;
    void writeVisItFields(VisItDataCollection&) const;

    const GridFunction& GetPotential() const { return *phi_; }
    const GridFunction& GetElectricField() const { return *e_; }

    double totalChargeFromRho() const;
    double totalCharge() const;
    
    double getCenterOfCharge() const;
    multipolarCoefficients getMultipolarCoefficients(std::size_t order) const;

    double chargeInBoundary(int bdrAttribute) const;
    double averagePotentialInDomain(int attr) const;
    double averagePotentialInBoundary(int bdrAttribute) const;
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
    void applyBoundaryConstantValuesToGridFunction(
        const AttrToValueMap& bdrValues, 
        GridFunction& gf) const;
};

}