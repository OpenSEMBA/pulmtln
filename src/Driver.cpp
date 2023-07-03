#include "Driver.h"

#include "ElectrostaticSolver.h"

namespace pulmtln {

using namespace mfem;

Driver::Driver(const nlohmann::json& input)
{   
    //Model model;
    //int sdim = model.mesh.SpaceDimension();
    //if (sdim != 2) {
    //    throw std::runtime_error("Solver can only run with 2D meshes.");
    //}
    //
    //std::map<int, double> domainToEpsr;
    //electrostaticSolver_ = 
    //    std::make_unique<ElectrostaticSolver>( 
    //        model_.mesh, model.dbc, domainToEpsr, opts_ );
    //
    //electrostaticSolver_->Solve();

    //ParaViewDataCollection paraview_dc{ "PULMTLN", &model_.mesh };
    //electrostaticSolver_->writeParaViewFields(paraview_dc); 
}

MTLPULParameters Driver::getMTLPUL() const
{
    MTLPULParameters res;

    // TODO

    return res;
}




}