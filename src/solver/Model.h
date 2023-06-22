#pragma once

#include "BoundaryConditions.h"

namespace pulmtln {

class Model {
public:
	Model() = default;
	Model(
		mfem::Mesh& m,
		const BoundaryConditions& dirichletBCs) :
		mesh{m},
		dbc{dirichletBCs}
	{}
	
	mfem::Mesh mesh;
	BoundaryConditions dbc;
};

}