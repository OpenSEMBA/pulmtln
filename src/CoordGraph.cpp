#include "CoordGraph.h"

#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <set>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>

namespace boost { template <typename G> inline void renumber_vertex_indices(G& g) { throw; } }

namespace pulmtln {

IdSet CoordGraph::getVertices() const 
{
    IdSet res;
    auto verts = vertices(graph_);
    for (auto& vp = verts; vp.first != vp.second; ++vp.first) {
        res.insert(graph_[*vp.first].id);
    }
    return res;
}

void CoordGraph::addVertex(const CoordinateId& id) 
{
    if (vertexMap_.find(id) == vertexMap_.end()) {
        vertexMap_.emplace(id, add_vertex(VertexId({ id }), graph_));
    }
}

void CoordGraph::addEdge(const CoordinateId& id1, const CoordinateId& id2) 
{
    if (id1 == id2) {
        throw std::runtime_error("Edges starting and finishing in same vertex are not allowed.");
    }
    
    if (vertexMap_.count(id1) == 0) {
        addVertex(id1);
    }
    if (vertexMap_.count(id2) == 0) {
        addVertex(id2);
    }
    
    add_edge(vertexMap_[id1], vertexMap_[id2], graph_);
}

struct cycle_recorder
{
    template <typename Path, typename graph_t>
    inline void cycle(const Path& p, const graph_t& g)
    { 
        cycles->push_back({});
        for (auto v : p) {
            cycles->back().push_back(g[v].id);
        }
    }
    std::vector<CoordGraph::Path>* cycles;
};

std::vector<Path> CoordGraph::findCycles() const
{
    std::vector<Path> res;
    cycle_recorder vis{&res};
    tiernan_all_cycles(graph_, vis);
    return res;
}

}