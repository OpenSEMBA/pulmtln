#include "Model.h"

#include <set>
#include <assert.h>

#include "DirectedGraph.h"

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
	DirectedGraph g;
	std::map<int, const Element*> vToE;
	for (const auto& e : elems) {
		assert(e->GetType() == Element::Type::SEGMENT);
		int v0{ e->GetVertices()[0] };
		int v1{ e->GetVertices()[1] };
		g.addEdge(v0, v1);
		vToE[v0] = e;
	}

	std::multimap<int, const Element*> res;
	int cycleCount{ 0 };
	for (const auto& cycle : g.findCycles()) {
		for (const auto& v : cycle) {
			res.emplace(cycleCount, vToE.at(v));
		}
		cycleCount++;
	}
	
	return res;
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