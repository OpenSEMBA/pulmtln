#pragma once

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <type_traits>

#include "core/geometry/mesh/Geometric.h"
#include "core/ProblemDescription.h"
#include "parsers/Parser.h"

namespace semba::parsers::JSON {

const std::string VERSION{ "0.16" };

using namespace geometry;
using namespace math;

using json = nlohmann::json;
using PM = physicalModel::PhysicalModel;

class Parser : public semba::parsers::Parser {
public:  
    Parser(const std::string& filename);
    UnstructuredProblemDescription read() const;
};

json readAnalysis(const json&);
std::unique_ptr<mesh::Unstructured> readUnstructuredMesh(const PMGroup&, const json&, const std::string& folder);
Grid3 readGrids(const json&);
PMGroup readMaterials(const json&);
SourceGroup readSources(mesh::Unstructured& mesh, const json&);
OutputRequestGroup readProbes(mesh::Unstructured& mesh, const json&);
void readBoundary(mesh::Unstructured& mesh, const json& j, PMGroup& physicalModelGroup, const Grid3& grid);


} 