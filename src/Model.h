#pragma once

#include "BdrConditionValues.h"
#include "Materials.h"

namespace pulmtln {

class Model {
public:
	Model() = default;
	Model(
		mfem::Mesh& mesh,
		const Materials& materials) :
		mesh_{mesh},
		materials_{materials}
	{}

	std::map<std::string, int> getMaterialsOfType(const MaterialType&);

	mfem::Mesh* getMesh() { return &mesh_; }
private:	
	Materials materials_;
	mfem::Mesh mesh_;
};

}