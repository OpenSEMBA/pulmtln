#include "Driver.h"

#include "ElectrostaticSolver.h"

#include <mfem.hpp>

namespace pulmtln {

using namespace mfem;

Driver::Driver(const Model& model, const SolverOptions& opts) :
    model_{model},
    opts_{opts}
{   
    int sdim = model.mesh.SpaceDimension();
    if (sdim != 2) {
        throw std::runtime_error("Solver can only run with 2D meshes.");
    }
    
    electrostaticSolver_ = 
        std::make_unique<ElectrostaticSolver>( 
            model_.mesh, model.dbc, opts_ );
    
    electrostaticSolver_->Assemble();
    electrostaticSolver_->Solve();

    ParaViewDataCollection paraview_dc{ "PULMTLN", &model_.mesh };
    electrostaticSolver_->writeParaViewFields(paraview_dc); 
}

const ElectrostaticSolver& Driver::getElectrostaticSolver() const
{
    return *electrostaticSolver_;
}

}