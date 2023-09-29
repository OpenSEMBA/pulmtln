#pragma once

#include "DirectedGraph.h"
#include "Model.h"

namespace pulmtln {

struct Domain {
    using Id = int;
    using ConductorId = int;
    using IdToDomain = std::map<Id, Domain>;

    ConductorId ground;
    std::set<ConductorId> conductorIds;
    std::set<int> elements;

    static IdToDomain buildDomains(const Model&);
    static DirectedGraph buildGraph(const IdToDomain&);
};

}