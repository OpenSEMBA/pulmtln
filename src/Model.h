#pragma once

#include "AttrToValueMap.h"
#include "Materials.h"

namespace pulmtln {

class Model {
public:
	enum class OpennessType {
		open,     // The most external boundary is open.
		semiopen, // The most external boundaries are open and a conductor.
		closed    // The most external boundary is a conductor.
	};

	Model() = default;
	Model(
		mfem::Mesh& mesh,          // Model gets ownership of mesh.
		const Materials& materials // Stores only materials present in mesh.
	);

	mfem::Mesh* getMesh() { return mesh_.get(); }
	const mfem::Mesh* getMesh() const { return mesh_.get(); }

	const Materials& getMaterials() const { return materials_; }
	
	void setGroundConductorId(MaterialId id) { groundConductorId_ = id; }
	MaterialId getGroundConductorId() const { return groundConductorId_; }

	OpennessType determineOpenness() const;

private:
	Materials materials_;
	std::unique_ptr<mfem::Mesh> mesh_;
	MaterialId groundConductorId_{ 0 };
};

}