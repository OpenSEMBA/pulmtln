#pragma once

#include "core/geometry/mesh/Unstructured.h"
#include "core/physicalModel/Group.h"

namespace semba {
namespace model {

using namespace geometry;

template<typename M = mesh::Unstructured>
class Model {
public:
	Model() = default;
	Model(const M&, const PMGroup&);
	Model(const Model& rhs);
	Model& operator=(const Model& rhs);

	M mesh;
	PMGroup physicalModels;
};

template<typename M>
Model<M>::Model(
	const M& mesh,
	const PMGroup& physicalModels
) : mesh(mesh),
physicalModels(physicalModels)
{
	this->mesh.reassignPointers(this->physicalModels);
}

template<typename M>
Model<M>::Model(const Model& rhs) {
	mesh = rhs.mesh;
	physicalModels = rhs.physicalModels;

	mesh.reassignPointers(physicalModels);
}

template<typename M>
Model<M>& Model<M>::operator=(const Model& rhs) {
	mesh = rhs.mesh;
	physicalModels = rhs.physicalModels;

	mesh.reassignPointers(physicalModels);

	return *this;
}

}
typedef model::Model<geometry::mesh::Unstructured> UnstructuredModel;
typedef model::Model<geometry::mesh::Structured> StructuredModel;
}