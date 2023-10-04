#pragma once

#include "AttrToValueMap.h"
#include "Materials.h"

namespace pulmtln {

class Model {
public:
	enum class OpennessType {
		open,
		semiopen,
		closed
	};

	Model() = default;
	Model(
		mfem::Mesh& mesh,     // Model gets ownership.
		const Materials& materials) :
		mesh_{std::make_unique<mfem::Mesh>(std::move(mesh))},
		materials_{materials}
	{}

	mfem::Mesh* getMesh() { return mesh_.get(); }
	const mfem::Mesh* getMesh() const { return mesh_.get(); }

	const Materials& getMaterials() const { return materials_;  }
	
	OpennessType determineOpenness() const;
	
private:	
	Materials materials_;
	std::unique_ptr<mfem::Mesh> mesh_;
};

}