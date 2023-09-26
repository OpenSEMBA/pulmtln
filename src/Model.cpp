#include "Model.h"

#include <set>
#include <assert.h>

namespace pulmtln {

using namespace mfem;

std::vector<const Element*> getElementsWithAttribute(const Mesh& mesh, int attr)
{	
	std::vector<const Element*> res;
	for (auto e{ 0 }; e < mesh.GetNBE(); ++e) {
		const auto elemPtr{ mesh.GetBdrElement(e) };
		if (elemPtr->GetAttribute() == attr) {
			res.push_back(elemPtr);
		}
	}
	return res;
}

std::multimap<int, const Element*> determineClosedLoops(const std::vector<const Element*>&elems)
{
	assert(std::all_of(
		elems.begin(), 
		elems.end(),
		[](const Element* e) {return e->GetType() == Element::Type::SEGMENT}
	));

	std::map<int, const Element*> vToE;
	for (const auto& e : elems) {
		std::set<int> segmentVertices;
		for (auto i{ 0 }; i < e->GetNVertices(); ++i) {
			segmentVertices.insert(e->GetVertices()[i]);
		}
		vToE.emplace(*segmentVertices.begin(), e);
	}

	for (auto it1{ vToE.begin() }; it1 != vToE.end(); ++it1) {
		auto it2{ std::next(it1) };
		if (it2 == vToE.end()) {
			it2 == vToE.begin();
		}
		
	}
}

bool elementsFormOpenLoops(const std::vector<const Element*>& elems)
{
	return determineClosedLoops(elems).size() != elems.size();
}

bool Model::isFullyOpen() const
{
	assert(materials_.openBoundaries.size() <= 1);
	
	for (const auto& m : materials_.openBoundaries) {
		if (!elementsFormOpenLoops(getElementsWithAttribute(mesh_, m.tag))) {
			return true;
		}
	}
	return false;
}


}