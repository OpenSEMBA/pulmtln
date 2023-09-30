#pragma once

#include "DirectedGraph.h"
#include "Model.h"

namespace pulmtln {

struct Domain {
    using Id = int;
    using ConductorId = int;
    using IdToDomain = std::map<Id, Domain>;

    ConductorId ground{-1};
    std::set<ConductorId> conductorIds;
    std::set<int> elements;

    static IdToDomain buildDomains(const Model&);

};


// DomainTree is a tree graph with a single root having
// a vertex for each domain.
// The root contains ConductorId 0.
class DomainTree : private DirectedGraph {
public:
    DomainTree() = default;
    DomainTree(const Domain::IdToDomain&);

    using DirectedGraph::getEdgesAsPairs;
};

}