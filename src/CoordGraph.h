#pragma once

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/function.hpp>

#include <iostream>
#include <vector>
#include <set>
#include <map>

namespace pulmtln {

using namespace boost;

using CoordinateId = int;
using Path = std::vector<CoordinateId>;
using IdSet = std::set<CoordinateId>;

struct VertexId { 
    CoordinateId id; 
    bool operator==(const VertexId& rhs) const { return id == rhs.id; }
    bool operator<(const VertexId& rhs) const { return id < rhs.id; }
};

typedef adjacency_list<vecS, vecS, bidirectionalS, VertexId, property<edge_index_t, int>> graph_t;

typedef filtered_graph<
    graph_t, 
    function<bool(graph_t::edge_descriptor)>, 
    function<bool(graph_t::vertex_descriptor)> > ComponentGraph;

typedef property_map<graph_t, vertex_index_t>::type IndexMap;
typedef VertexId vertex_t;

typedef graph_traits<graph_t>::vertex_iterator vertex_iter;
/// @brief Implements an oriented graph class based on <a href="https://www.boost.org/doc/libs/1_78_0/libs/graph/doc/index.html">the BOOST Graph Library</a> 
class CoordGraph {
public:
    typedef std::pair<CoordinateId, CoordinateId> EdgeIds;
    typedef std::map<CoordinateId, graph_t::vertex_descriptor> VertexMap;
    typedef std::vector<CoordinateId> Path;
    typedef std::vector<Path> Paths;

    CoordGraph() = default;
    
    /// @brief Add vertex identified by CoordinateId id to the graph.
    void addVertex(const CoordinateId& id);
    /// @brief Add oriented edge from vertex id1 to vertex id2 to the graph. 
    void addEdge(const CoordinateId& id1, const CoordinateId& id2);
    
    /// @brief Returns the number of vertices in the graph. 
    std::size_t verticesSize() const { return num_vertices(graph_); }
    /// @brief return the number of edges in the graph.
    std::size_t edgesSize() const { return num_edges(graph_); }

    /// @brief Returns the CoordinateId corresponding to the graph vertices. 
    IdSet getVertices() const;
    
    /// @brief Returns all the Path in the graph that form a cycle.
    Paths findCycles() const;
    
private:
    VertexMap vertexMap_;
    graph_t graph_;
};


}
