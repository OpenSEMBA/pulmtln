#include "Driver.h"

#include "Solver.h"
#include "constants.h"

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
    
    std::cout  << "Starting initialization." << std::endl;
    Solver pulmtln{ model_, opts_ };
    
    ParaViewDataCollection paraview_dc{ "PULMTLN", &model_.mesh };
    pulmtln.RegisterParaViewFields(paraview_dc);
    std::cout  << "Initialization done." << std::endl;

    // Assemble all forms
    pulmtln.Assemble();

    // Solve the system and compute any auxiliary fields
    pulmtln.Solve();

    pulmtln.WriteParaViewFields();   
}

}