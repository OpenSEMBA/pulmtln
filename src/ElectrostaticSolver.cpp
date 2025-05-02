#include "ElectrostaticSolver.h"

#include "constants.h"
#include "AttrToValueMap.h"

namespace pulmtln {

Array<int> AttrToMarker(const Mesh& mesh, const Array<int>& attrs)
{
    if (attrs.Size() == 0) {
        return Array<int>();
    }

    int maxAttrInConditions{*std::max_element(attrs.begin(), attrs.end())};
    int maxAttrInMesh{ mesh.bdr_attributes.Max() };
    int maxAttr{ std::max({maxAttrInConditions, maxAttrInMesh}) };

    Array<int> marker(maxAttr);
    marker = 0;

    
    MFEM_ASSERT(attrs.Max() <= maxAttr, "Invalid attribute number present.");

    if (attrs.Size() == 1 && attrs[0] == -1)
    {
        marker = 1;
    }
    else
    {
        marker = 0;
        for (int j = 0; j < attrs.Size(); j++)
        {
            int attr = attrs[j];
            MFEM_VERIFY(attr > 0, "Attribute number less than one!");
            marker[attr - 1] = 1;
        }
    }

    return marker;
}

Array<int> toArray(const std::vector<int>& v)
{
    Array<int> r{ (int)v.size() };
    for (int i{ 0 }; i < r.Size(); i++) {
        r[i] = v[i];
    }
    return r;
}

double firstOrderABC(const Vector& rVec)
{
    double r{ rVec.Norml2() };

    return - EPSILON0_NATURAL / (r * std::log(r));
}


ElectrostaticSolver::ElectrostaticSolver(
    Mesh& mesh,
    const SolverParameters& parameters,
    const SolverOptions opts) : 
    opts_(opts),
    mesh_(&mesh),
    parameters_(parameters),
    H1FESpace_(NULL),
    HCurlFESpace_(NULL),
    HDivFESpace_(NULL),
    L2FESpace_(NULL),
    divEpsGrad_(NULL),
    hDivMass_(NULL),
    hCurlHDivEps_(NULL),
    h1SurfMass_(NULL),
    rhod_(NULL),
    grad_(NULL),
    phi_(NULL),
    rho_(NULL),
    sigma_src_(NULL),
    e_(NULL),
    d_(NULL),
    oneCoef_(1.0)
{
    
    // Define compatible parallel finite element spaces on the parallel
    // mesh. Here we use arbitrary order H1, Nedelec, and Raviart-Thomas finite
    // elements.
    auto order{ opts_.order };
    H1FESpace_ = new H1_FESpace(mesh_, order, mesh_->Dimension());
    HCurlFESpace_ = new ND_FESpace(mesh_, order, mesh_->Dimension());
    HDivFESpace_ = new RT_FESpace(mesh_, order, mesh_->Dimension());
    L2FESpace_ = new L2_FESpace(mesh_, order - 1, mesh_->Dimension());

    // Select surface attributes for Dirichlet BCs
    ess_bdr_ = AttrToMarker(
        *mesh_,
        parameters.dirichletBoundaries.getAttributesAsArray()
    );

    // Setup domain permittivity coefficients.
    if (parameters_.domainPermittivities.empty()) {
        epsCoef_ = new ConstantCoefficient(EPSILON0_NATURAL);
    } else {
        mfem::Vector eps(mesh_->attributes.Max());
        eps = EPSILON0_NATURAL;
        for (const auto& [attr, epsr] : parameters_.domainPermittivities) {
            const auto epsSize{ eps.Size() };
            assert(attr <= epsSize);
            eps[attr-1] *= epsr;
        }
        epsCoef_ = new PWConstCoefficient(eps);
    }
    
    // Bilinear Forms
    divEpsGrad_ = new BilinearForm(H1FESpace_);
    divEpsGrad_->AddDomainIntegrator(new DiffusionIntegrator(*epsCoef_));

    // Setup open regions.
    std::unique_ptr<Coefficient> openRegionCoeff;
    if (!parameters.openBoundaries.empty()) {
        open_bdr_ = AttrToMarker(*mesh_, toArray(parameters.openBoundaries));
        openRegionCoeff.reset(new FunctionCoefficient(firstOrderABC));
        divEpsGrad_->AddBoundaryIntegrator(
            new BoundaryMassIntegrator(*openRegionCoeff),
            open_bdr_
        );
    }

    hDivMass_ = new BilinearForm(HDivFESpace_);
    hDivMass_->AddDomainIntegrator(new VectorFEMassIntegrator);
        
    h1SurfMass_ = new BilinearForm(H1FESpace_);
    h1SurfMass_->AddBoundaryIntegrator(new MassIntegrator);
    
    hCurlHDivEps_ = new MixedBilinearForm(HCurlFESpace_, HDivFESpace_);
    hCurlHDivEps_->AddDomainIntegrator(new VectorFEMassIntegrator(*epsCoef_));

    rhod_ = new LinearForm(H1FESpace_);

    // Discrete derivative operator
    grad_ = new DiscreteGradOperator(H1FESpace_, HCurlFESpace_);
    div_ = new DiscreteDivOperator(HDivFESpace_, L2FESpace_);

    // Build grid functions
    phi_ = new GridFunction(H1FESpace_);
    d_ = new GridFunction(HDivFESpace_);
    e_ = new GridFunction(HCurlFESpace_);
    rho_ = new GridFunction(L2FESpace_);
    sigma_src_ = new GridFunction(H1FESpace_);

    Assemble();
}

ElectrostaticSolver::~ElectrostaticSolver()
{   
    delete phi_;
    delete rho_;
    delete rhod_;
    delete d_;
    delete e_;
    delete sigma_src_;
    
    delete grad_;
    delete div_;

    delete divEpsGrad_;
    delete hDivMass_;
    delete hCurlHDivEps_;
    delete h1SurfMass_;
    
    delete H1FESpace_;
    delete HCurlFESpace_;
    delete HDivFESpace_;
    delete L2FESpace_;

    delete epsCoef_;
}

void ElectrostaticSolver::Assemble()
{
    divEpsGrad_->Assemble();
    divEpsGrad_->Finalize();

    hDivMass_->Assemble();
    hDivMass_->Finalize();

    hCurlHDivEps_->Assemble();
    hCurlHDivEps_->Finalize();

    h1SurfMass_->Assemble();
    h1SurfMass_->Finalize();

    *rhod_ = 0.0;
    rhod_->Assemble();

    grad_->Assemble();
    grad_->Finalize();

    div_->Assemble();
    div_->Finalize();
}

void ElectrostaticSolver::applyBoundaryValuesToGridFunction(
    const AttrToValueMap& bdrValues,
    GridFunction& gf) const
{
    auto attributes{ bdrValues.getAttributesAsArray() };
    auto values{ bdrValues.getValuesAsArray() };
    if (attributes.Size() == 0) {
        return;
    }

    // Apply piecewise constant boundary condition
    Array<int> bdr_attr(mesh_->bdr_attributes.Max());
    for (int i = 0; i < attributes.Size(); i++)
    {
        ConstantCoefficient voltage(values[i]);
        bdr_attr = 0;
        if (attributes[i] <= bdr_attr.Size())
        {
            bdr_attr[attributes[i] - 1] = 1;
        }
        gf.ProjectBdrCoefficient(voltage, bdr_attr);
    }
}

void ElectrostaticSolver::Solve()
{
    // Initialize the surface charge density (Neumann boundaries).
    *sigma_src_ = 0.0;
    applyBoundaryValuesToGridFunction(parameters_.neumannBoundaries, *sigma_src_);
    h1SurfMass_->AddMult(*sigma_src_, *rhod_);
    
    // Solves phi (electrostatic potential).
    *phi_ = 0.0; 
    {
        auto dbcs{ parameters_.dirichletBoundaries.getAttributesAsArray() };
        applyBoundaryValuesToGridFunction(parameters_.dirichletBoundaries, *phi_);

        // Determine the essential BC degrees of freedom
        if (dbcs.Size() > 0) {
            H1FESpace_->GetEssentialTrueDofs(ess_bdr_, ess_bdr_tdofs_);
        }
        else {
            ess_bdr_tdofs_.SetSize(1);
            ess_bdr_tdofs_[0] = 0;
        }


        // Apply essential BC and form linear system
        SparseMatrix DivEpsGrad;
        Vector Phi;
        Vector RHS;

        divEpsGrad_->FormLinearSystem(
            ess_bdr_tdofs_, *phi_, *rhod_, 
            DivEpsGrad, Phi, RHS);

        GSSmoother M(DivEpsGrad);
        PCG(DivEpsGrad, M, RHS, Phi, opts_.printIterations, 500, 1e-12, 0.0);

        divEpsGrad_->RecoverFEMSolution(Phi, *rhod_, *phi_);
    }

    // Computes E.
    grad_->Mult(*phi_, *e_); *e_ *= -1.0;
    
    // Computes D.
    {
        GridFunction ed(HDivFESpace_);
        hCurlHDivEps_->Mult(*e_, ed);
    
        SparseMatrix MassHDiv;
        Vector ED, D;

        Array<int> dbc_dofs_d;
        hDivMass_->FormLinearSystem(dbc_dofs_d, *d_, ed, MassHDiv, D, ED);

        GSSmoother M(MassHDiv);
        PCG(MassHDiv, M, ED, D, opts_.printIterations, 200, 1e-12);
        
        hDivMass_->RecoverFEMSolution(D, ed, *d_);
    }
    
    // Computes rho.
    div_->Mult(*d_, *rho_);
    
}

double ElectrostaticSolver::totalChargeFromRho() const
{
    LinearForm l2_vol_int{L2FESpace_};
    ConstantCoefficient oneCoef{1.0};
    l2_vol_int.AddDomainIntegrator(new DomainLFIntegrator(oneCoef));
    l2_vol_int.Assemble();
    return l2_vol_int(*rho_);
}

double ElectrostaticSolver::totalCharge() const
{
    LinearForm rt_surf_int{HDivFESpace_};
    rt_surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator);
    rt_surf_int.Assemble();
    return rt_surf_int(*d_);
}

std::unique_ptr<LinearForm> buildSurfaceIntegratorForBoundary(RT_FESpace* fes, int bdrAttribute)
{
    mfem::Array<int> attr(1);
    attr[0] = bdrAttribute;

    Array<Coefficient*> coefs(attr.Size());
    ConstantCoefficient one{ -1.0 };
    for (int i = 0; i < coefs.Size(); i++) {
        coefs[i] = &one;
    }
    PWCoefficient pwcoeff{ attr, coefs };

    auto surf_int =  std::make_unique<LinearForm>(fes);
    surf_int->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
    surf_int->Assemble();

    return surf_int;
}

double ElectrostaticSolver::chargeInBoundary(int bdrAttribute) const
{
    auto surf_int{ buildSurfaceIntegratorForBoundary(HDivFESpace_, bdrAttribute) };
    return (*surf_int)(*d_);
}

double ElectrostaticSolver::totalEnergy() const
{
    BilinearForm mass{ HCurlFESpace_ };
    ConstantCoefficient eps{ EPSILON0_NATURAL };
    mass.AddDomainIntegrator(new VectorFEMassIntegrator(eps));
    mass.Assemble();
    
    GridFunction aux(e_->FESpace());
    mass.Mult(*e_, aux);

    double energy{ 0.0 };
    energy = 0.5*InnerProduct(*e_, aux);

    return energy;
}

void ElectrostaticSolver::writeParaViewFields(ParaViewDataCollection& pv) const
{
    pv.SetHighOrderOutput(true);
    pv.SetLevelsOfDetail(3);
    pv.RegisterField("Phi", phi_);
    pv.RegisterField("D", d_);
    pv.RegisterField("E", e_);
    //pv.RegisterField("Rho", rho_);

    pv.Save();
}

void ElectrostaticSolver::writeVisItFields(
    VisItDataCollection& pv) const
{
    pv.SetLevelsOfDetail(3);
    pv.RegisterField("Phi", phi_);
    pv.RegisterField("D", d_);
    pv.RegisterField("E", e_);
    //pv.RegisterField("Rho", rho_);

    pv.Save();
}

}