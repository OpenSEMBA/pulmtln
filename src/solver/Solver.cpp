#include "Solver.h"

#include "constants.h"

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

Solver::Solver(
    Model& model, 
    const SolverOptions& opts) : 
    opts_(opts),
    mesh_(&model.mesh),
    dbcs_(&model.dbcs),
    dbcv_(&model.dbcv),
    nbcs_(&model.nbcs),
    nbcv_(&model.nbcv),
    paraview_dc_(NULL),
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
    AttrToMarker(mesh_->bdr_attributes.Max(), *dbcs_, ess_bdr_);

    // Setup various coefficients
    epsCoef_ = ConstantCoefficient(epsilon0_);
    
    // Bilinear Forms
    divEpsGrad_ = new BilinearForm(H1FESpace_);
    divEpsGrad_->AddDomainIntegrator(new DiffusionIntegrator(epsCoef_));

    hDivMass_ = new BilinearForm(HDivFESpace_);
    hDivMass_->AddDomainIntegrator(new VectorFEMassIntegrator);

    hCurlHDivEps_ = new MixedBilinearForm(HCurlFESpace_, HDivFESpace_);
    hCurlHDivEps_->AddDomainIntegrator(new VectorFEMassIntegrator(epsCoef_));

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

    if (nbcs_->Size() > 0)
    {
        h1SurfMass_ = new BilinearForm(H1FESpace_);
        h1SurfMass_->AddBoundaryIntegrator(new MassIntegrator);
    }
}

Solver::~Solver()
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

}

void Solver::Assemble()
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

    std::cout << "done." << std::endl << std::flush; 
}

void Solver::Update()
{
    std::cout << "Updating ..." << std::endl;

    H1FESpace_->Update(false);
    HCurlFESpace_->Update(false);
    HDivFESpace_->Update(false);
    L2FESpace_->Update(false);

    // Inform the grid functions that the space has changed.
    phi_->Update();
    rhod_->Update();
    l2_vol_int_->Update();
    rt_surf_int_->Update();
    d_->Update();
    e_->Update();
    rho_->Update();

    // Inform the bilinear forms that the space has changed.
    divEpsGrad_->Update();
    hDivMass_->Update();
    hCurlHDivEps_->Update();

    if (h1Mass_) { h1Mass_->Update(); }
    if (h1SurfMass_) { h1SurfMass_->Update(); }
    if (hCurlHDiv_) { hCurlHDiv_->Update(); }
    if (weakDiv_) { weakDiv_->Update(); }

    // Inform the other objects that the space has changed.
    grad_->Update();
    div_->Update();
}


void Solver::Solve()
{
    std::cout << "Running solver ... " << std::endl;

    // Initialize the electric potential with its boundary conditions
    *phi_ = 0.0;
    {
        if (dbcs_->Size() > 0) {
            // Apply piecewise constant boundary condition
            Array<int> dbc_bdr_attr(mesh_->bdr_attributes.Max());
            for (int i = 0; i < dbcs_->Size(); i++)
            {
                ConstantCoefficient voltage((*dbcv_)[i]);
                dbc_bdr_attr = 0;
                if ((*dbcs_)[i] <= dbc_bdr_attr.Size())
                {
                    dbc_bdr_attr[(*dbcs_)[i] - 1] = 1;
                }
                phi_->ProjectBdrCoefficient(voltage, dbc_bdr_attr);
            }
        }

        // Determine the essential BC degrees of freedom
        if (dbcs_->Size() > 0)
        {
            // From user supplied boundary attributes
            H1FESpace_->GetEssentialTrueDofs(ess_bdr_, ess_bdr_tdofs_);
        }
        else
        {
            // Use the first DoF on processor zero by default
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

        // Extract the parallel grid function corresponding to the finite
        // element approximation Phi. This is the local solution on each
        // processor.
        divEpsGrad_->RecoverFEMSolution(Phi, *rhod_, *phi_);
    }
    // Compute the negative Gradient of the solution vector.  This is
    // the magnetic field corresponding to the scalar potential
    // represented by phi.
    grad_->Mult(*phi_, *e_); *e_ *= -1.0;
    
    // Compute electric displacement (D) from E and P (if present)
    std::cout  << "Computing D ..." << std::flush;

    GridFunction ed(HDivFESpace_);
    hCurlHDivEps_->Mult(*e_, ed);
    
    SparseMatrix MassHDiv;
    Vector ED, D;

    Array<int> dbc_dofs_d;
    hDivMass_->FormLinearSystem(dbc_dofs_d, *d_, ed, MassHDiv, D, ED);

    GSSmoother diagM(MassHDiv);
    PCG(MassHDiv, diagM, ED, D, 0, 200, 1e-12);
    diagM.Mult(ED, D);
    
    hDivMass_->RecoverFEMSolution(D, ed, *d_);
    div_->Mult(*d_, *rho_);

    std::cout  << "done." << std::flush;

    {
        // Compute total charge as volume integral of rho
        double charge_rho = (*l2_vol_int_)(*rho_);

        // Compute total charge as surface integral of D
        double charge_D = (*rt_surf_int_)(*d_);

        std::cout << std::endl << "Total charge: \n"
            << "   Volume integral of charge density:   " << charge_rho
            << "\n   Surface integral of dielectric flux: " << charge_D
            << std::endl << std::flush;
    }

    std::cout << "Solver done. " << std::endl; 
}

void Solver::RegisterParaViewFields(ParaViewDataCollection& pv)
{
    paraview_dc_ = &pv;

    pv.RegisterField("Phi", phi_);
    pv.RegisterField("D", d_);
    pv.RegisterField("E", e_);
    pv.RegisterField("Rho", rho_);
}

void Solver::WriteParaViewFields()
{
    if (!paraview_dc_) {
        throw std::runtime_error("Paraview register not initialized.");
    }
    std::cout << "Writing ParaView files ..." << std::flush;
    paraview_dc_->Save();

    std::cout << " done." << std::endl;
}

}