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
    /// @brief CoordGraph class constructor. 
    
    /// @brief Add vertex identified by CoordinateId id to the graph.
    void addVertex(const CoordinateId& id);
    /// @brief Add oriented edge from vertex id1 to vertex id2 to the graph. 
    void addEdge(const CoordinateId& id1, const CoordinateId& id2);
    /// @brief Remove vertex identified by CoordinateId id from the graph.
    void removeVertex(const CoordinateId& id);
    /// @brief Remove oriented edge from vertex id1 to vertex id2 to the graph. 
    void removeEdge(const CoordinateId& id1, const CoordinateId& id2);
    /// @brief Returns the number of vertices in the graph. 
    std::size_t verticesSize() const { return num_vertices(graph_); }
    /// @brief return the number of edges in the graph.
    std::size_t edgesSize() const { return num_edges(graph_); }

    /// @brief Return a graph for each disconnected subgraph within the graph. 
    std::vector<CoordGraph> split() const;

    /// @brief Returns the CoordinateId corresponding to the graph vertices. 
    IdSet getVertices() const;
    std::vector<CoordinateId> getOrderedVertices() const;
    /// @brief Returns the CoordinateId of the vertices adjacent to CoordianteId id.
    IdSet getAdjacentVertices(const CoordinateId id) const;
    /// @brief Returns the CoordinateId of the vertices in the interior of a graph 
    /// with a single interior region. 
    IdSet getInterior() const;
    /// @brief Returns the CoordinateId of the vertices in the extremes of a graph, 
    /// defined as those vertices with 0 or 1 neighbours. 
    IdSet getExterior() const;
    /// @brief Returns the CoordinateId of the vertices on the boundary of the graph, 
    /// and the vertices not on the boundary of the graph. 
    std::pair<IdSet, IdSet> getBoundAndInteriorVertices() const;

    IdSet getClosestVerticesInSet(const CoordinateId& id, const IdSet& coordSet) const;

    /// @brief Returns all the Path in the graph that form a cycle.
    Paths findCycles() const;
    /// @brief Returns true if the Path forms a cycle and has a definite orientation. 
    bool isOrientableAndCyclic(const Path&) const;

    /// @brief Returns all the Path in the graph that do not form a cycle, and whose vertices
    /// do not belong to any cycle.
    Paths findAcyclicPaths() const;
    /// @brief Returns all the edges in the graph that are not part of a cycle.
    std::set<EdgeIds> getAcyclicEdges() const;
    /// @brief Transforms the input Path into an oriented Path according to the graph edge orientation.
    Path orderByOrientation(const Path&) const;
    /// @brief Return the shortest Path in the graph connecting Coordinateid ini and CoordinateId end,
    Path findShortestPath(const CoordinateId& ini, const CoordinateId& end) const;

    /// @brief Returns the boundary of the graph, formed by all the external edges. 
    CoordGraph getBoundaryGraph() const;
    /// @brief Returns the graph formed by all the internal edges. 
    CoordGraph getInternalGraph() const;
    /// @brief Returns the intersection of the graph with another graph.
    CoordGraph intersect(const CoordGraph& rhs) const;
    /// @brief Returns the intersection of the graph with a set of CoordinateId.
    template<typename Container> CoordGraph intersect(const Container& rhs) const;
    /// @brief Returns the difference of the graph with another graph. 
    CoordGraph difference(const CoordGraph& rhs) const;
    
    bool canBeSplit() const;

private:
    VertexMap vertexMap_;
    graph_t graph_;

    std::vector<std::vector<VertexId>> findCycles_() const;
    IdSet getExtremes() const;

    static Paths removeRepeated(const Paths&);
    bool isForwardOriented(const Path&) const;
};


}