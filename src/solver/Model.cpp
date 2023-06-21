#include "Model.h"

namespace pulmtln {

Model::Model(
	Mesh& mesh, 
	const AttributeToMaterial& matMap, 
	const AttributeToBoundary& bdrMap) :
	mesh_(mesh)
{
	if (matMap.size() == 0) {
		attToMatMap_.emplace(1, Material(1.0, 1.0));
	}
	else {
		attToMatMap_ = matMap;
	}

	if (bdrMap.size() == 0) {
		for (int i = 1; i <= mesh.bdr_attributes.Size(); i++) {
			attToBdrMap_.emplace(i, BdrCond::PEC);
		}
	}
	else {
		attToBdrMap_ = bdrMap;
	}

	assembleAttToTypeMap(attToBdrMap_, bdrToMarkerMap_);
}

std::size_t Model::numberOfMaterials() const
{
	return attToMatMap_.size();
}

std::size_t Model::numberOfBoundaryMaterials() const
{
	return attToBdrMap_.size();
}

void Model::assembleAttToTypeMap(
	std::map<Attribute, BdrCond>& attToCond,
	std::multimap<BdrCond, BoundaryMarker>& attToMarker)
{
	for (const auto& kv : attToCond) {
		const auto& att{ kv.first };
		const auto& bdr{ kv.second };
		assert(att > 0);

		BoundaryMarker bdrMarker{ mesh_.bdr_attributes.Max() };
		bdrMarker = 0;
		bdrMarker[att - 1] = 1;

		attToMarker.emplace(bdr, bdrMarker);
	}
}


}
