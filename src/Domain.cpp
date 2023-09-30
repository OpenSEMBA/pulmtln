#include "Domain.h"

namespace pulmtln {

DirectedGraph buildMeshGraph(const mfem::Mesh& mesh)
{
	DirectedGraph meshGraph;
	for (auto e{ 0 }; e < mesh.GetNE(); ++e) {
		const mfem::Element* elem{ mesh.GetElement(e) };
		DirectedGraph::Path vs(elem->GetNVertices());
		for (auto i{ 0 }; i < vs.size(); ++i) {
			vs[i] = elem->GetVertices()[i];
		}
		meshGraph.addClosedPath(vs);
	}

	return meshGraph;
}

bool isPECInDomain(const PEC& pec, const IdSet& verticesInDomain, const mfem::Mesh& mesh)
{
	IdSet verticesWithAttribute;
	for (auto e{ 0 }; e < mesh.GetNBE(); ++e) {
		const mfem::Element* elem{ mesh.GetBdrElement(e) };
		if (elem->GetAttribute() != pec.attribute) {
			continue;
		}
		verticesWithAttribute.insert(
			elem->GetVertices(),
			elem->GetVertices() + elem->GetNVertices());
	}

	IdSet common;
	std::set_intersection(
		verticesInDomain.begin(), verticesInDomain.end(),
		verticesWithAttribute.begin(), verticesWithAttribute.end(),
		std::inserter(common, common.begin())
	);

	return common.size() > 0;
}
	
DomainTree::DomainTree(const Domain::IdToDomain& domains)
{
	// Check conductor 0 must be in a single domain.
	bool foundOnce{ false };
	for (const auto& [id, dom] : domains) {
		if (dom.conductorIds.count(0) && !foundOnce) {
			foundOnce = true;
			continue;
		}
		if (dom.conductorIds.count(0) && foundOnce) {
			throw std::runtime_error(
				"Conductor 0 must be only in a single domain."
			);
		}
	}
	if (!foundOnce) {
		throw std::runtime_error("Conductor 0 is not present");
	}
	
	std::multimap<Domain::ConductorId, Domain::Id> condIdToDomId;
	for (const auto& [domId, dom] : domains) {
		addVertex(domId);
		for (const auto& condId : dom.conductorIds) {
			condIdToDomId.emplace(condId, domId);
		}
	}

	int prevCondId{ -1 };
	int prevDomId{ -1 };
	for (const auto& [cId, dId] : condIdToDomId) {
		if (cId == prevCondId) {
			addEdge(prevDomId, dId);
		}
		else {
			prevCondId = cId;
			prevDomId = dId;
			addVertex(dId);
		}
	}

	// Post-conditions
	assert(findCycles().size() == 0);
}

Domain::IdToDomain Domain::buildDomains(const Model& model)
{
	if (model.determineOpenness() != Model::OpennessType::closed) {
		throw std::runtime_error("Domains can only be determined for closed problems.");
	}

	const mfem::Mesh& mesh{ *model.getMesh() };

	Domain::IdToDomain res;
	Domain::Id id{ 0 };
	for (const auto& domainMeshGraph : buildMeshGraph(mesh).split()) {

		const auto vsInDomain{ domainMeshGraph.getVertices() };

		Domain domain;
		// Determine which elements belong to domain.
		for (auto e{ 0 }; e < mesh.GetNE(); ++e) {
			const auto& elem{ *mesh.GetElement(e) };
			if (std::all_of(
				elem.GetVertices(), 
				elem.GetVertices() + elem.GetNVertices(),
				[&vsInDomain](int i) { return vsInDomain.count(i) == 1; }
			)) {
				domain.elements.insert(e);
			}
		}

		// Determine conductors in domain.
		for (const auto& pec : model.getMaterials().pecs) {
			if (isPECInDomain(pec, vsInDomain, mesh)) {
				domain.conductorIds.insert(
					Materials::getNumberContainedInName(pec.name) 
				);
			}
		}

		//
		res[id++] = domain;
	}

	// Sets grounds.
	for (auto& [domId, dom] : res) {
		if (dom.conductorIds.count(0) == 1) {
			dom.ground = 0;
		}
	}
	for (const auto& edge : DomainTree{ res }.getEdgesAsPairs()) {
		const auto& c1{ res[edge.first].conductorIds };
		const auto& c2{ res[edge.second].conductorIds };
		std::set<ConductorId> common;
		std::set_intersection(
			c1.begin(), c1.end(),
			c2.begin(), c2.end(),
			std::inserter(common, common.begin())
		);
		assert(common.size() == 1);
		res[edge.second].ground = *common.begin();
	}
	
	// Post conditions.
	// All elements belong to a single domain.
	// Conductor 0 is always ground and is in a single domain.

	return res;
}


}