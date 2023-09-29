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
}

bool isPECInDomain(const PEC& pec, const IdSet& verticesInDomain, const mfem::Mesh& mesh)
{
	// TODO
	return false;
}
	
DirectedGraph buildGraph(const Domain::IdToDomain& domains)
{
	int mostExternalDomain{ -1 };
	for (auto& [id, dom] : domains) {
		if (dom.conductorIds.count(0) == 1) {
			mostExternalDomain = id;
			break;
		}
	}

	// TODO

	DirectedGraph res;

	return res;
}

Domain::IdToDomain Domain::buildDomains(const Model& model)
{
	if (model.determineOpenness() != Model::OpennessType::closed) {
		throw std::runtime_error("Domains can only be determined for closed problems.");
	}

	const mfem::Mesh& mesh{ *model.getMesh() };

	Domain::IdToDomain res;
	Domain::Id id{ 0 };
	for (const auto& domainGraph : buildMeshGraph(mesh).split()) {
		Domain domain;
		auto vsInDomain{ domainGraph.getVertices() };

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
	auto domainGraph{ Domain::buildGraph(res) };


	


	// Post conditions.
	// All elements belong to a single domain.
	// If domain contains conductor 0. That one is always ground.

	return res;
}

DirectedGraph Domain::buildGraph(const IdToDomain&)
{
	// Builds graph in which vertices represent domains with edges pointing to subdomains.
	// Starts with domain 0.
		
	DirectedGraph res;

	// TODO

	return res;
}

}