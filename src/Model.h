#pragma once

#include "AttrToValueMap.h"
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

	mfem::Mesh* getMesh() { return &mesh_; }
	const mfem::Mesh* getMesh() const { return &mesh_; }

	const Materials& getMaterials() const { return materials_;  }
private:	
	Materials materials_;
	mfem::Mesh mesh_;
};

}