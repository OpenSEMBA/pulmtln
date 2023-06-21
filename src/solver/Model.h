#pragma once

#include "Material.h"

#include <mfem.hpp>
#include <map>

namespace pulmtln {

enum class BdrCond {
	PEC,
	PMC
};

using Attribute = int;
using AttributeToMaterial = std::map<Attribute, Material>;
using AttributeToBoundary = std::map<Attribute, BdrCond>;

using BoundaryMarker = mfem::Array<int>;
using BoundaryToMarker = std::multimap<BdrCond, BoundaryMarker>;

class Model {
public:
	using Mesh = mfem::Mesh;
	Model() = default;
	Model(
		Mesh&, 
		const AttributeToMaterial& = AttributeToMaterial{},
		const AttributeToBoundary& = AttributeToBoundary{}
	);

	Mesh& getMesh() { return mesh_; };
	const Mesh& getConstMesh() const { return mesh_; }
	
	BoundaryToMarker& getBoundaryToMarker() { return bdrToMarkerMap_; }
	const BoundaryToMarker& getBoundaryToMarker() const { return bdrToMarkerMap_; }
	
	std::size_t numberOfMaterials() const;
	std::size_t numberOfBoundaryMaterials() const;
private:
	Mesh mesh_;
	
	AttributeToMaterial attToMatMap_;
	AttributeToBoundary attToBdrMap_;
	BoundaryToMarker bdrToMarkerMap_;
	
	void assembleAttToTypeMap(
		std::map<Attribute, BdrCond>& attToCond, 
		std::multimap<BdrCond, BoundaryMarker>& attToMarker);
};

}