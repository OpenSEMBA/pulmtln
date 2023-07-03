#include "ElectrostaticSolver.h"

#include "constants.h"
#include "BoundaryConditions.h"

namespace pulmtln {

void AttrToMarker(int max_attr, const Array<int>& attrs, Array<int>& marker)
{
    MFEM_ASSERT(attrs.Max() <= max_attr, "Invalid attribute number present.");

    marker.SetSize(max_attr);
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
}

ElectrostaticSolver::ElectrostaticSolver(
    Mesh& mesh,
    const BoundaryConditions& dbc,
    const std::map<int, double>& domainToEpsr,
    const SolverOptions& opts) : 
    opts_(opts),
    mesh_(&mesh),
    dbc_(dbc),
    domainToEpsr_(domainToEpsr),
    H1FESpace_(NULL),
    HCurlFESpace_(NULL),
    HDivFESpace_(NULL),
    L2FESpace_(NULL),
    divEpsGrad_(NULL),
    h1Mass_(NULL),
    h1SurfMass_(NULL),
    hDivMass_(NULL),
    hCurlHDivEps_(NULL),
    hCurlHDiv_(NULL),
    weakDiv_(NULL),
    rhod_(NULL),
    l2_vol_int_(NULL),
    rt_surf_int_(NULL),
    grad_(NULL),
    phi_(NULL),
    rho_(NULL),
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
    AttrToMarker(mesh_->bdr_attributes.Max(), dbc.getAttributes(), ess_bdr_);

    // Setup various coefficients
    if (domainToEpsr_.size() == 0) {
        epsCoef_ = new ConstantCoefficient(epsilon0_);
    }
    else {
        mfem::Vector eps(mesh_->attributes.Max());
        eps = epsilon0_;
        for (const auto& [attr, epsr] : domainToEpsr_) {
            assert(attr <= eps.Size());
            eps[attr-1] *= epsr;
        }
        epsCoef_ = new PWConstCoefficient(eps);
    }
    
    // Bilinear Forms
    divEpsGrad_ = new BilinearForm(H1FESpace_);
    divEpsGrad_->AddDomainIntegrator(new DiffusionIntegrator(*epsCoef_));

    hDivMass_ = new BilinearForm(HDivFESpace_);
    hDivMass_->AddDomainIntegrator(new VectorFEMassIntegrator);

    hCurlHDivEps_ = new MixedBilinearForm(HCurlFESpace_, HDivFESpace_);
    hCurlHDivEps_->AddDomainIntegrator(new VectorFEMassIntegrator(*epsCoef_));

    rhod_ = new LinearForm(H1FESpace_);

    
    l2_vol_int_ = new LinearForm(L2FESpace_);
    l2_vol_int_->AddDomainIntegrator(new DomainLFIntegrator(oneCoef_));

    rt_surf_int_ = new LinearForm(HDivFESpace_);
    rt_surf_int_->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator);

    // Discrete derivative operator
    grad_ = new DiscreteGradOperator(H1FESpace_, HCurlFESpace_);
    div_ = new DiscreteDivOperator(HDivFESpace_, L2FESpace_);

    // Build grid functions
    phi_ = new GridFunction(H1FESpace_);
    d_ = new GridFunction(HDivFESpace_);
    e_ = new GridFunction(HCurlFESpace_);
    rho_ = new GridFunction(L2FESpace_);

    Assemble();
}

ElectrostaticSolver::~ElectrostaticSolver()
{   
    delete phi_;
    delete rho_;
    delete rhod_;
    delete l2_vol_int_;
    delete rt_surf_int_;
    delete d_;
    delete e_;
    
    delete grad_;
    delete div_;

    delete divEpsGrad_;
    delete h1Mass_;
    delete h1SurfMass_;
    delete hDivMass_;
    delete hCurlHDivEps_;
    delete hCurlHDiv_;
    delete weakDiv_;

    delete H1FESpace_;
    delete HCurlFESpace_;
    delete HDivFESpace_;
    delete L2FESpace_;

    delete epsCoef_;
}

void ElectrostaticSolver::Assemble()
{
    std::cout << "Assembling ... " << std::flush; 

    divEpsGrad_->Assemble();
    divEpsGrad_->Finalize();

    hDivMass_->Assemble();
    hDivMass_->Finalize();

    hCurlHDivEps_->Assemble();
    hCurlHDivEps_->Finalize();

    *rhod_ = 0.0;
    rhod_->Assemble();

    l2_vol_int_->Assemble();
    rt_surf_int_->Assemble();

    grad_->Assemble();
    grad_->Finalize();

    div_->Assemble();
    div_->Finalize();

    if (h1Mass_)
    {
        h1Mass_->Assemble();
        h1Mass_->Finalize();
    }
    if (h1SurfMass_)
    {
        h1SurfMass_->Assemble();
        h1SurfMass_->Finalize();
    }
    if (hCurlHDiv_)
    {
        hCurlHDiv_->Assemble();
        hCurlHDiv_->Finalize();
    }
    if (weakDiv_)
    {
        weakDiv_->Assemble();
        weakDiv_->Finalize();
    }
    std::cout << "done." << std::endl; 
}

void ElectrostaticSolver::Solve()
{
    std::cout << "Running solver ... " << std::endl;

    std::cout << "Computing phi ..." << std::flush;
    *phi_ = 0.0;
    {
        auto dbcs{ dbc_.getAttributes() };
        auto dbcv{ dbc_.getValues() };
        if (dbcs.Size() > 0) {
            // Apply piecewise constant boundary condition
            Array<int> dbc_bdr_attr(mesh_->bdr_attributes.Max());
            for (int i = 0; i < dbcs.Size(); i++)
            {
                ConstantCoefficient voltage(dbcv[i]);
                dbc_bdr_attr = 0;
                if (dbcs[i] <= dbc_bdr_attr.Size())
                {
                    dbc_bdr_attr[dbcs[i] - 1] = 1;
                }
                phi_->ProjectBdrCoefficient(voltage, dbc_bdr_attr);
            }
        }

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
        PCG(DivEpsGrad, M, RHS, Phi, 1, 200, 1e-12, 0.0);

        divEpsGrad_->RecoverFEMSolution(Phi, *rhod_, *phi_);
    }

    std::cout << "Computing E ..." << std::flush;
    grad_->Mult(*phi_, *e_); *e_ *= -1.0;
    std::cout << "done." << std::endl;

    std::cout  << "Computing D ..." << std::flush;
    {
        GridFunction ed(HDivFESpace_);
        hCurlHDivEps_->Mult(*e_, ed);
    
        SparseMatrix MassHDiv;
        Vector ED, D;

        Array<int> dbc_dofs_d;
        hDivMass_->FormLinearSystem(dbc_dofs_d, *d_, ed, MassHDiv, D, ED);

        GSSmoother M(MassHDiv);
        PCG(MassHDiv, M, ED, D, 1, 200, 1e-12);
        
        hDivMass_->RecoverFEMSolution(D, ed, *d_);
    }
    std::cout << "done." << std::endl;

    std::cout << "Computing rho ..." << std::flush;
    div_->Mult(*d_, *rho_);
    std::cout  << "done." << std::endl;
    
    std::cout << "Solver done. " << std::endl; 
}

double ElectrostaticSolver::totalChargeFromRho() const
{
    return (*l2_vol_int_)(*rho_);
}

double ElectrostaticSolver::totalCharge() const
{
    return (*rt_surf_int_)(*d_);
}

double ElectrostaticSolver::chargeInBoundary(const Array<int>& attr) const
{
    Array<Coefficient*> coefs(attr.Size());
    ConstantCoefficient one{ -1.0 };
    for (int i = 0; i < coefs.Size(); i++) {
        coefs[i] = &one;
    }
    PWCoefficient pwcoeff{ attr, coefs };

    LinearForm surf_int(HDivFESpace_);
    surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
    surf_int.Assemble();
    auto charge{ surf_int(*d_) };

    return charge;
}

void ElectrostaticSolver::writeParaViewFields(
    ParaViewDataCollection& pv) const
{
    pv.SetHighOrderOutput(true);
    pv.SetLevelsOfDetail(3);
    pv.RegisterField("Phi", phi_);
    pv.RegisterField("D", d_);
    pv.RegisterField("E", e_);
    pv.RegisterField("Rho", rho_);

    pv.Save();
}

}